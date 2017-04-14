#!/usr/bin/python
# -*- coding: utf8 -*-

#~ import psycopg2 as sql
#~ import getpass
import rec_to_db
import numpy as np
import sys, getopt, os, shutil
import re, copy, random
import tree2
import blockevents
import time
import shelve
# call remider
callline = 'python getblockvents.py ...'

# default flags
writeResults=True
saveObjectDicts=True
silent=True
# default parameters
buildancs = False
dbclue = 'phylarianne'
dbpwd = None
queue = None
dirxmltrees = None
minblocksize = 2
bsthreshold = 0.90
maxgapsize = {'transfer':4, 'speciation':0, 'duplication':0, 'gain':1}
skipeventtypes = {'speciation':'all'}
commonsetoperation = 'intersection'
repliconchecktable = 'public.buildingblocks'
lreplicon = None
randomseed = None
splitdiscordant=True	
starttime = time.time()
tl = 0
safetl = 0
prevelapsedtime = 0


# default fields for SQL queries
queryFields = ['gene_id', 'chromosome', 'genomic_beg', 'genomic_end', 'hogenom_gene_id', 'name', 'locus_tag', 'family_accession']
statFields = ['nblock', 'repliconid', 'compatible', 'not_against', 'against']
		
		
def shuffleReplicon(tprot, randomseed):
	"""create a false random order of genes on the replicon"""
	if randomseed:
		# saves original prtein order
		nonshuffledtprot = copy.deepcopy(tprot)
		# randomizes the protein order
		random.seed(randomseed)
		lprot = list(tprot)
		random.shuffle(lprot)
		tprot = tuple(lprot)
	else:
		nonshuffledtprot = None
	return (tprot, nonshuffledtprot)

def checkProtAroundReplicon(tprot, nonshuffledtprot, replicon, nblock, reftree, dirtreepickles, cn, maxgapsize=maxgapsize, skipeventtypes=skipeventtypes, silent=silent, randomseed=randomseed):
	"""look systematically around the replicon for proteins with an event and build blocks from those seeds
	
	code very similar to blockevents.LeafBlock.searchNgbrEvent(), BOTH SHOULD BE FACTORIZED INTO ONE FONCTION.
	"""
	dprot_events = {}
	#~ i = 0
	#~ while i < len(tprot):
	for i in range(len(tprot)):
		if not silent:
			#~ sys.stderr.write("\r%d"%i)
			print 'checkProtAroundReplicon(), i =', i, "; next prot:",
		apid = tprot[i][0]
		# loads the gene tree
		dgene = rec_to_db.getGeneInfo(apid, queryFields, cn)
		# check for existence of the gene tree
		apfam = dgene['family_accession']
		if not silent: print "%s (%s)"%(apid, apfam)
		nfgenetreepickle = "%s/%s.pickle"%(dirtreepickles, apfam)
		if not os.access(nfgenetreepickle, os.F_OK):
			# may be a small (n<4) family with no tree computed
			tdfamily = rec_to_db.getFamilyInfo(apfam, '*', cn, returnDict=False)
			if len(tdfamily) < 4:
				#~ i += 1
				if not silent: print "\tskip prot from small family", apfam
				continue
			else:
				raise IndexError, "family %s has no tree computed according to database record:\n%s"%(apfam, str(tdfamily))
		genetree = tree2.load_pickle(nfgenetreepickle)
		# get the first events of each type located at nodes above the gene in the reconciled gene tree
		#~ ddevents = genetree.getEvents(lineage=apid, returnDict=True)
		#~ lendindex = []
		lineage = genetree[apid].lineage(value='id')
		recgiid = rec_to_db.getGeneTreeInfo(apfam, ['rec_gi_id'], cn, returnDict=False)[0]
		if not silent: print apid, "lineage:", lineage
		for startnodeid in lineage:
			eventnode = genetree.idgetnode(startnodeid)
			tdevents = rec_to_db.getEventInfo(recgiid, startnodeid, ['event_id', 'event_type', 'rec_gi_end_node'], cn, returnDict=True)
			# for eventnodelab in ddevents:
				# eventnode = genetree[eventnodelab]
				# devent = ddevents[eventnodelab]		# only one event per node, ignores possible conflicting transfer inferences by Prunier at the same node
			for devent in tdevents:
				eventtype = devent['eventtype']
				eventid = devent['event_id']
				# check if event was not already recorded in another block
				if apid in dprot_events:
					if eventid in dprot_events[apid]:
						# already incorpored in a block as a compatible protein
						if not silent: print "\tstartnodeid: %d, eventid: %d ; already incorpored in a block"%(startnodeid, eventid)
						continue
				devent['eventlocation'] = rec_to_db.getEventLocation(eventid, eventtype, cn)
				if eventtype in skipeventtypes:
					# builds block only from chosen event types
					if 'all' in skipeventtypes[eventtype]:
						# all events of this type are skipped
						continue
					if (set(devent['eventlocation'][1]) & set(skipeventtypes[eventtype])):
						# the location of the event intersects with those for which this type of events are skipped
						continue
					if 'old' in skipeventtypes[eventtype]:
						# all events of this type older than or as old as the gene tree ancestor are skipped
						treeanc = reftree.map_to_node(genetree.dictLeafLabelsToSpecies().values())
						oldercoord = reftree.coalesce(devent['eventlocation'][1])
						if oldercoord==reftree or (oldercoord==treeanc or treeanc.is_child(oldercoord)):
							continue
				#~ elif eventtype=='transfer':
					#~ if rec_to_db.getTransferObservation(eventid, cn)<1:
						#~ # ignore transfers not supported by a Prunier replicate for bloc initiation
						#~ continue
					#~ elif eventnode==genetree:
						#~ # ignore artefactual transfers annotated at the root of the gene tree
						#~ continue
				if not silent: print "\tstartnodeid: %d, devent: %s"%(startnodeid, str(devent))
				seed = blockevents.ProtInBlock( dfields=dgene, eventnode=eventnode, dicevent=devent)
				if randomseed:
					# changes protein coordinates of the shuffled replicon into those of the ordered replicon (for compatibility with position dependent functions)
					dref = rec_to_db.getGeneInfo(nonshuffledtprot[i][0], queryFields, cn)
					seed.editCoords(dref['begin'], dref['end'])
	#						print seed
				# recognition of the side of the event
				#~ side = seed.getEventSide()
				side = None
				if eventtype == 'transfer':
					trchildid = devent['eventlocation'][4]
					trchild = genetree.idgetnode(trchildid)
					if trchild:
						if genetree[seed.getid()]==trchild or genetree[seed.getid()].is_child(trchild):
							side = 'rec'
						else:
							side = 'don'
				seed.setStatus('seed')
				## initiation of list values in block dictionaries
				nblock += 1
				block = blockevents.LeafBlock(nblock, seed, side=side)
				replicon.addBlock(block)
				## Searching on replicon for neighbour genes with similar transfer patterns
				n = block.searchNgbrEvent(tprot, i, +1, queryFields, maxgapsize, reftree, cn, dirtreepickles, nonshuffledtprot, silent=silent)
				if not silent: print "new block", block, '\ni' , i, 'n', n
				# trim the non-transfered gaps not closed by transfered proteins at extremities
				(ntrimleft, ntrimright) = block.trimNonTransferedEnds(maxgapsize=maxgapsize, silent=silent)
				if not silent: print "trimmed block", block, '\nntrimright', ntrimright,'\n'
				#~ # block ended at n (end of block exploration) - ntrimright (number of trimmed proteins)
				#~ lendindex.append(n - ntrimright)
				## sets the definitive block coordinates in replicon record
				replicon.addBlockCoords(block)
				## record the explored events for the proteins in the block
				dblockprotevents = block.getDictProtToEventids()
				for protid in dblockprotevents:
					dprot_events.setdefault(protid, []).append(dblockprotevents[protid])
		#~ if lendindex:
			#~ # continues the loop after the last protein added to the shorter block
			#~ i = min(lendindex)
		#~ else:
			#~ i += 1
	return nblock

def checkGapsAndSplitBlocks(replicon, nblock, reftree, dirtreepickles, ct, foutprestat, maxgapsize=maxgapsize, bsthreshold=bsthreshold, splitdiscordant=splitdiscordant, writeResults=writeResults, silent=silent):
	def checkGap(block, nblock):
		lnewid = []
		eventtype = block.eventtype()
		if eventtype != 'transfer':
			# gap cannot bear signal against a duplication or speciation (can't they ???), skip this check
			return lnewid, nblock
#			print block
		lnewblocks = []
		# searches for proteins with signal against block transfer history
		lagainst = block.checkGapProtSignal(reftree, bsthreshold, dirtreepickles, ct)
		if lagainst and not silent: print "lagainst:", [agp.getid() for agp in lagainst]
		# makes statistics before spliting blocks and withdrawing against protein
		if writeResults:				
			# counts protein status in block
			(ncompat, nnotagainst, nagainst) = block.countStatus()
			# writing block statistics
			foutprestat.write('%s\t%s\t%d\t%d\t%d\n'%(block.getid(), replicon.getid(), ncompat, nnotagainst, nagainst))							
		# splits the blocks according to list of 'against' proteins
		if splitdiscordant:
			maxsplit = 1
			# if found at least one 'against' protein
			if lagainst:
				# deletes block coordinates from replicon record before modifying them
				replicon.delBlockCoords(block, silent=silent)
				while lagainst:
					# splits the block at this point
					lnb, nblock = block.splitDiscordantBlock(lagainst, maxgapsize, nblock, maxsplit)
					lnewblocks += lnb
					if not silent: print "lnewblocks:", [lnb.getid() for lnb in lnewblocks]
					# and refreshes the block rec/don sets and status of the block proteins
					block.commonRecDonSet()
					lagainst = block.checkGapProtSignal(reftree, bsthreshold, dirtreepickles, ct)
					if not silent: print "lagainst:", [agp.getid() for agp in lagainst]
				else:
					# reloads block coordinates in replicon records
					replicon.addBlockCoords(block, silent=silent)
					# appends the new blocks yielded by block splits to the end of replicon list ; they will be treated similarly at the end of the loop.
					for newblock in lnewblocks:
						replicon.addBlock(newblock)
						replicon.addBlockCoords(newblock, silent=silent)
						lnewid.append(newblock.getid())
		# when the block is clean, uses information from 'not_against' proteins to refine rec/don sets
		block.commonRecDonSet(refine=True)
		return lnewid, nblock
		
	lnewblockid = []
	for block in replicon:
		lni, nblock = checkGap(block, nblock)
		lnewblockid += lni
	# test to see if misses to check new blocks
	while lnewblockid:
		newblockid = lnewblockid.pop(0)
		lni, nblock = checkGap(replicon[newblockid], nblock)
		lnewblockid += lni
	if splitdiscordant:	
		# after having split blocks, might find new redundancy on replicons
		purge = replicon.purgeBlockList(reftree, silent=True)
		#~ if purge:
			#~ raise IndexError, "Blocks should not be redundant" 
	return nblock

def buildAncestralBlocks(replicon, nances, nblock, dancblocks, dfam_blocks, dfam_ancs, skipeventtypes=skipeventtypes, commonsetoperation=commonsetoperation, silent=silent):
	newleafblocks = []	# list of new leaf blocks created when disassembling uncoherent blocks
	for block in replicon:
		eventtype = block.eventtype()
		if eventtype in skipeventtypes:
			# builds block only from chosen event types
			if 'all' in skipeventtypes[eventtype]:
				# all events of this type are skipped
				continue
			if block.getRecSet() & set(skipeventtypes[eventtype]):
				# the location of the event intersects with those for which this type of events are skipped
				continue
		if not silent: print "\n", block.description()
		nances, nblock = block.assignBlock(nances, nblock, dancblocks, dfam_blocks, dfam_ancs, commonsetoperation=commonsetoperation, silent=silent)
	return nances, nblock
	
def parseReplicon(repliconid, nblock, dirtreepickles, reftree, dirout, foutstat, foutprestat, db, jobid, skipeventtypes=skipeventtypes, minblocksize=minblocksize, repliconchecktable=repliconchecktable, queue=queue, tl=tl, safetl=safetl, prevelapsedtime=prevelapsedtime, starttime=starttime, randomseed=randomseed, silent=silent, callline=callline):	
	# create session cursors
	cr=db.cursor()
	cp=db.cursor()
	cn=db.cursor()
	ct=db.cursor()
	
	# time check
	currtime = time.time()
	elapsedtime = currtime - starttime
	print "elapsed time", elapsedtime
	if queue and (elapsedtime > safetl) :
		# better stop at 90% of time limit
		lastrepli = repliconid
		return 1, nblock
	if queue and (elapsedtime - prevelapsedtime) > (tl - elapsedtime):
		# time to parse one replicon is higher than remaining time, must stop
		lastrepli = repliconid
		return 1, nblock
	prevelapsedtime = elapsedtime
	
	cr.execute( "UPDATE %s SET job_id=%s, current_status=1, date=now() WHERE chromosome='%s';"%(repliconchecktable, str(jobid), str(repliconid)))
	db.commit()
	
	print 'Scanning %s'%repliconid
	replicon = blockevents.RepliconMap(repliconid)
	# queries all the proteins
	cp.execute( "SELECT DISTINCT hogenom_gene_id, genomic_beg FROM genome.gene WHERE chromosome=%s ORDER BY genomic_beg;", (repliconid,) )
	tprot = cp.fetchall()	# tuple of all protein information tuples matching the query
	(tprot, nonshuffledtprot) = shuffleReplicon(tprot, randomseed)
	
	# looping over proteins (seeds for blocks)
	print 'Building blocks'
	nblock = checkProtAroundReplicon(tprot, nonshuffledtprot, replicon, nblock, reftree, dirtreepickles, cn, maxgapsize=maxgapsize, skipeventtypes=skipeventtypes, randomseed=randomseed, silent=silent)
	
	print 'Purging redundant blocks'
	purge = replicon.purgeBlockList(reftree, silent=silent)
	#~ if purge:
		#~ raise IndexError, "Blocks should not be redundant" 	
	
	print "Checking gap protein for signal against / not against block transfer history" 		
	nblock = checkGapsAndSplitBlocks(replicon, nblock, reftree, dirtreepickles, ct, foutprestat, maxgapsize=maxgapsize, bsthreshold=bsthreshold, splitdiscordant=splitdiscordant, silent=silent)
				
	if writeResults:	
		dirblockout = '%s/blocks'%dirout
		fout = open('%s/blocks_%s'%(dirblockout, repliconid), 'w')
		fout.write('### call : %s\n### %d blocks ; detail of blocks of size >= %d:\n\n'%(callline, len(replicon.getBlockList()),minblocksize))	
		for block in replicon:
			# counts protein status in block
			(ncompat, nnotagainst, nagainst) = block.countStatus()
			# writing block statistics
			foutstat.write('%s\t%s\t%d\t%d\t%d\n'%(block.getid(), repliconid, ncompat, nnotagainst, nagainst))		
			if block.getBlockSize() >= minblocksize:
				# filter donor-side transfer blocks
				if block.eventtype()=='transfer' and block.eventside()!='rec':
					continue
				# writing block description
				block.writeBlockToFile(fout)
		fout.close()
#			print '\tLeaf blocks description file written at \n\t%s/blocks_%s'%(dirblockout, repliconid)		
	
	#~ # fills leaf block dictionary with blocks from this replicon	
	#~ for blockid in replicon.getBlockList():
		#~ if not blockid in dblocks:
			#~ dblocks[blockid] = replicon[blockid]
		#~ else:
			#~ print "old entry\n", dblocks[blockid].description()
			#~ print "new entry\n", replicon[blockid].description()
			#~ raise IndexError, "block id %s should be unique"%str(blockid)

	if saveObjectDicts:
		print 'Saving Pyhton objects in shelve...'
		repliconshelve = shelve.open("%s/replicons.shelve"%dirout, flag='c')
		repliconshelve[repliconid] = replicon
		repliconshelve.close()
		print "RepliconMap object saved in %s/replicons.shelve at entry '%s'"%(dirout, str(repliconid))
	else:
		# fills replicon dictionary
		dreplicons[repliconid] = replicon

	# state end of the job in control table
	cr.execute( "UPDATE %s SET  job_id=%s, current_status=2, date=now() WHERE chromosome='%s';"%(repliconchecktable, str(jobid), str(repliconid)))
	db.commit()
	print "\n- - - - - - \n"
	
	return 0, nblock
	# end of replicon loop 1 (construction of leaf blocks) ; individual replicon iteration can be split in several jobs

	
def usage():
	s =  "Usage : python getblockevents.py genetreeobjects reftree outputdir jobid [options]\n"
	s += " Arguments:\n"
	s += "  genetreeobjects (directory path)\n\tdirectory containing the pickled reconciled gene tree objects\n"
	s += "  reftree (file path)\n\treference phylogeny in Newick format.\n"
	s += "  outputdir (directory path)\n\tdirectory where to write ouput files.\n"
	s += "  jobid (integer)\n\tnumber multiplied by one billion to initiate the blocks identifiers built in this job (to avoid id redundancy between jobs). For instance, '5' will initiate the counter at 5000000000.\n"
	s += " Options:\n"
	s += "  -a \n\twill build ancestral blocks once leaf block contruction is complete.\n"
	s += "  -q integer\n\t maximum duration time for farming the job (in hours), execution will stop when reaching > 85% of it.\n"
	s += "  -t directory_path\n\tdirectory containing all trees of the database in phyloXML format.\n"
	s += "  -R file_path | replicon1[,replicon2[,...]]\n\tlist of replicons to be treated\n"
	s += "  -m integer\t(default=2)\n\tthe minimal size (in number of genes involved in an event) to display the block in output\n"
	s += "  -b float\t(default=90)\n\tthe bootstrap threshold for significancy of a branch when testing the signal against a transfer event in a gene tree\n"
	s += "  -g dict\t(default={transfer:2, speciation:0, duplication:0, gain:0})\n\tthe maximal gap size for transfered gene block extension\n\tas a (bracketted) comma-separated list of event types detailing type-wise maximum gaps\n\t\te.g. \"8\" or \"transfer:2,speciation:0\" or \"{transfer:5,speciation:0,duplication:0}\"\n"
	s += "  -e dict\t(default=empty)\n\ta (bracketted) comma-separated list of the set of event types to be skipped for block construction\n\ta dot-separated list of node labels in species tree where the event should be overlooked can be appended with a colon separator\n\t\te.g. \"speciation:N1.N2.N3,duplication:N1\"\n"
	#~ s += "  -o boolean\n\tthe type of operation used for determining common sets of receptor or donor clades in a block : 'union' (0) or 'intersection' (1, default)\n"
	s += "  -s boolean\n\tspecifies if leaf blocks are to be split at proteins with \"against\" status (1, default) or not (0)\n"
	s += "  -r integer\n\tif provided, any non-null integer argument will be used as random seed to shuffle the gene order\n"	
	s += "  -c str\t(default='public.buildingblocks')\n\tthe table in agrogenom PostgreSQL database where to check for replicon to screen\n"
	s += "  -v\n\tverbose mode\n"
	s += "  -h\n\tprint this help message and exit\n"
	s += " Phylogenomic database connection options:\n"
	s += "  -d str\n\tPostgreSQL database guess string\n"
	s += "  -D str\n\tPostgreSQL database name\n"
	s += "  -H str\n\tPostgreSQL host server\n"
	s += "  -U str\n\tPostgreSQL user\n"
	s += "  -p str\n\tPostgreSQL user password\n"
	return s

def main(buildancs=buildancs, dbpwd=dbpwd, dbclue=dbclue, queue=queue, dirxmltrees=dirxmltrees, minblocksize=minblocksize, bsthreshold=bsthreshold, maxgapsize=maxgapsize, skipeventtypes=skipeventtypes, repliconchecktable=repliconchecktable, commonsetoperation=commonsetoperation, lreplicon=lreplicon, randomseed=randomseed, splitdiscordant=splitdiscordant, starttime=starttime, tl=tl, safetl=safetl, silent=silent):
#~ def main(writeResults=writeResults, saveObjectDicts=saveObjectDicts, maxgapsize=maxgapsize, skipeventtypes=skipeventtypes, commonsetoperation=commonsetoperation):

	if len(sys.argv) < 5:
		print "Missing arguments\n", usage()
		sys.exit(2)	
		
	dirtreepickles = sys.argv[1]
	if not os.access(dirtreepickles, os.F_OK):
		print 'wrong path to directory containing the pickled reconciled gene tree objects'
		sys.exit(2)

	nfreftree = sys.argv[2]
	if not os.access(nfreftree	, os.F_OK):
		print 'wrong path to refererence tree file'
		sys.exit(2)

	outdir = sys.argv[3]
		
	jobid = int(sys.argv[4])

	callline = 'python '+' '.join(sys.argv)

	try:
		options, args = getopt.getopt(sys.argv[5:], 'avp:q:t:R:m:b:g:e:o:r:c:s:hd:D:H:U:')
	except getopt.GetoptError, err:
		print str(err)
		print usage()
		sys.exit(2)	
	dopt = dict(options)

	if '-h' in dopt:
		print usage()
		sys.exit(0)
	buildancs = ('-a' in dopt)
	silent = not ('-v' in dopt)
	dbname = dopt.get('-D')
	dbuser = dopt.get('-U')
	dbhost = dopt.get('-H')
	dbpwd = dopt.get('-p')
	dbclue = dopt.get('-d', 'phylarianne')
	if dbpwd: callline.replace(dbpwd, '*****')
	queue = int(dopt.get('-q', 0))
	dirxmltrees = dopt.get('-t')
	minblocksize = int(dopt.get('-m', 2))
	bsthreshold = float(dopt.get('-b', 0.9))
	if '-g' in dopt: 
		spdmgs = dopt['-g'].strip(' ').split(',')
		for dmgs in spdmgs:
			spdmgs = dmgs.strip(' ').split(':')
			if len(spdmgs) == 2:
				maxgapsize[spdmgs[0]] = spdmgs[1]
			else:
				print "unadequate syntax for specifying maximum gap sizes"
				print usage()
				sys.exit(2)	
	print 'maxgapsize:', maxgapsize
	if '-e' in dopt:
		lskipeventtypes = dopt['-e'].split(',')
		for skipe in lskipeventtypes:
			spskipe = skipe.strip(' ').split(':')
			if len(spskipe)==1:
				skipeventtypes[skipe] = ['all']
			elif len(spskipe)==2:
				skipeventtypes[spskipe[0]] = spskipe[1].split('.')
			else:
				print "unadequate syntax for specifying events to be skipped"
				print usage()
				sys.exit(2)
	print 'skipeventtypes:', skipeventtypes
	repliconchecktable = dopt.get('-c', 'public.buildingblocks')
	#~ if '-o' in dopt: 
		#~ if dopt['-o']=='0': 
			#~ commonsetoperation = 'union'
		#~ elif dopt['-o']=='1': 
			#~ commonsetoperation = 'intersection'
		#~ else: 
			#~ print "Wrong set operation definition: %s"%dopt['-o']
			#~ print usage()
			#~ sys.exit(2)	
	print 'commonsetoperation:', commonsetoperation
	if '-R' in dopt: 
		repliconlist = dopt['-R']
		if os.access( repliconlist, os.F_OK ):
			frepliconlist = open(repliconlist, 'r')
			for line in frepliconlist:
				lreplicon.append(line.strip('\n '))
			frepliconlist.close()
		else:
			lreplicon = repliconlist.split(',')
	if '-r' in dopt: 
		randomseed=int(dopt['-r'])
		print "Protein order will be randomized on replicons. Random seed = %d"%randomseed
	splitdiscordant = bool(int(dopt.get('-s', True)))
		
	starttime = time.time()
	print "Start job %d at %s"%(jobid, time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
	if queue:
		print "queue = %dh"%queue
		tl = float(queue * 3600) 
		safetl = tl * 0.90
		print "time limit = %fs"%tl
		print "safe time limit = %fs"%safetl
		

	# sets connection objects to AGROGENOMDB	
	db, cj = rec_to_db.dbconnect(dbname=dbname, dbuser=dbuser, dbhost=dbhost, dbpwd=dbpwd, dbclue=dbclue)
		
	if buildancs:
		dirout = "%s/getblockevent_out.a.%d"%(outdir,jobid)
	else:
		dirout = "%s/getblockevent_out.l.%d"%(outdir,jobid)

	dirconcattrees = '%s/block_trees'%dirout
	dirblockout = '%s/blocks'%dirout
	dirsqlout = "%s/blocks_db_dump"%dirout

	nfstat = '%s/transfered_blocks.stats'%dirout
	nfprestat = '%s/before_split_blocks.stats'%dirout
	nfancblocks = "%s/ancestral_blocks"%dirout
	nfancstats = '%s/ancestral_blocks.stats'%dirout
	nfmatrix = "%s/blocks_events.mat"%dirout

	# loads reference tree (species tree)
	reftree = tree2.ReferenceTree(fic=nfreftree)
	reftree.complete_node_ids()

	if os.access(dirout, os.F_OK):
		print "output directory %s already exists. Risk to overwrite existing files. Exit."%dirout
		sys.exit(2)
	for d in [dirout, dirblockout, dirconcattrees, dirsqlout]:
		os.mkdir(d)

	if writeResults:
		foutstat = open(nfstat, 'w')
		foutstat.write('\t'.join(statFields)+'\n')
		foutprestat = open(nfprestat, 'w')
		foutprestat.write('\t'.join(statFields)+'\n')

	nblock = jobid * 1000000000
	nances = jobid * 1000000000
				
	dreplicons = {}

	lblocks = []
	dblocks = {}
	dfam_blocks = {'transfer':{}, 'speciation':{}, 'duplication':{}, 'gain':{}}

	dancblocks = {}
	dfam_ancs= {'transfer':{}, 'speciation':{}, 'duplication':{}, 'gain':{}}

	if not buildancs:
			
		## construction of leaf blocks
		print "Building leaf blocks"
		# initiate looping over replicons	
		if not lreplicon:
			cj.execute("LOCK %s ;"%repliconchecktable)
			cj.execute( "SELECT chromosome FROM %s WHERE job_id=0 AND current_status=0 ORDER by nb_genes desc;"%repliconchecktable)
			treplicon = cj.fetchone()
		else:
			treplicon = (lreplicon.pop(),)

		# initiate clock
		currtime = time.time()
		lastrepli = None
		prevelapsedtime = 0
			
		# looping over replicons	
		while treplicon:
			repliconid = treplicon[0]
			rcode, nblock = parseReplicon(repliconid, nblock, dirtreepickles, reftree, dirout, foutstat, foutprestat, db, jobid, skipeventtypes=skipeventtypes, minblocksize=minblocksize, repliconchecktable=repliconchecktable, queue=queue, tl=tl, safetl=safetl, prevelapsedtime=prevelapsedtime, starttime=starttime, silent=silent, randomseed=randomseed)
			if rcode > 0:
				break
			if not lreplicon:
				if lreplicon==None:
					cj.execute("LOCK %s ;"%repliconchecktable)
					cj.execute( "SELECT chromosome FROM %s WHERE job_id=0 and current_status=0 ORDER by nb_genes desc;"%repliconchecktable)
					treplicon = cj.fetchone()
				else:
					if len(lreplicon) > 0:
						treplicon = (lreplicon.pop(),)
					else:
						break

		print 'Leaf block description files written at %s'%dirblockout	
		if writeResults:	
			foutstat.close()
			foutprestat.close()
			print 'Block statistics written at \n%s'%nfstat

		if lastrepli:
			print "Not enough remaining time for engaging a block construction campaign on a new replicon.\n Stopped after %s"%lastrepli
			print "End of job %d at %s"%(jobid, time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))	
			
	# elif buildancs:
	else:
		## construction of ancestral blocks
		print "Building ancestral blocks"
		#~ if saveObjectDicts:
		ljobdir = os.listdir(outdir)
		jobdirids = []
		for jobdir in ljobdir:
			# check existing job ids
			jobdirids.append(int(jobdir.split('.')[-1]))
		# reset first leaf block identifier if redundant with existing ones (needed as in case of multiple association of a leaf block with ancestral blocks, a new leaf block is created)
		if not jobid > max(jobdirids):
			nblock = (max(jobdirids)+1) * 1000000000
		print "nblock initiated at", nblock
		for jobdir in ljobdir:
			# loads saved Python objects
			nfshelve = "%s/%s/replicons.shelve"%(outdir, jobdir)
			print nfshelve
			if os.access(nfshelve, os.F_OK):
				print "Open RepliconMap object shelve : %s"%nfshelve
				repliconshelve = shelve.open(nfshelve, flag='r')
				for repliconid in repliconshelve:
					replicon = repliconshelve[repliconid]
					print "\tadds blocks from %s"%repliconid
					nances, nblock = buildAncestralBlocks(replicon, nances, nblock, dancblocks, dfam_blocks, dfam_ancs, skipeventtypes, commonsetoperation, silent=silent)
				repliconshelve.close()	
		#~ else:
			#~ # use objects from current memory
			#~ for repliconid in dreplicons:
				#~ replicon = dreplicons[repliconid]
				#~ print "\tadds blocks from %s"%repliconid
				#~ nances = buildAncestralBlocks(replicon, nances, nblock, dancblocks, dfam_blocks, dfam_ancs, skipeventtypes, commonsetoperation, silent=silent)
		# end of replicon loop 2 (construction of ancestral blocks) ; individual replicon iteration increment on the same objects, every iteration must be executed on the same job/memory
		print ""
		
		if saveObjectDicts:
			print 'Saving Python objects in shelve...'
			ancblockshelve = shelve.open("%s/ancblocks.shelve"%dirout)
			for ancblockid in dancblocks:
				ancblockshelve[str(ancblockid)] = dancblocks[ancblockid]
	#			print "\tAncestralBlock objects saved in %s/ancblocks.shelve at entry '%s'"%(dirout, str(ancblockid))
			ancblockshelve.close()
			print "Saved all AncestralBlock objects"
		if writeResults:	
			print 'Writing ancestral blocks descriptions...'
			blockevents.writeAncestralBlockEvents(dancblocks, reftree, nfancblocks, nfmatrix, nfancstats)	
			print 'Ancestral blocks description file written at \n%s'%nfancblocks
			print 'Block transfers events written in clade matrix \n%s'%nfmatrix
			print 'Writing PostgreSQL database dump...'
			blockevents.writeSQLBlockTables(dancblocks, dirsqlout, reftree, cj)
			print 'agrogenom PostgreSQL database dump written at \n%s'%dirsqlout
			if dirxmltrees:
				blockevents.writeXMLBlockTreeCollections(dancblocks, dirxmltrees, dirconcattrees)
				print 'Concatenate trees of block families written in directory \n%s'%dirconcattrees
		

	print "End of job %d at %s"%(jobid, time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))	
		
if __name__ == '__main__':
	main()
