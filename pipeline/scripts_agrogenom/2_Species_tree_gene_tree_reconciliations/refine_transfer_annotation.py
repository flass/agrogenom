#!/usr/bin/python
# -*- coding: utf8 -*-

import sys, getopt
import os, shutil
import tree2
import rec_to_db
from numpy import mean
from blockevents import compatibleRecDonSet
import time
cwd = os.getcwd()

#~ def chooseEvent():

def main(dirgenetrees, nfreftree, dirout, reccolid, testsession=True, appendmode=False, **kw):
	
	silent = kw.get('silent', True)
	useTempTables = kw.get('useTempTables', 2)
	minbs = float(kw.get('minbs', 0.9))	
	tables = ['phylogeny.event', 'phylogeny.event_possible_species_node', "phylogeny.represented_event"]
	
	if not appendmode:
		if os.access(dirout, os.F_OK):
			#~ shutil.rmtree(dirout)
			if not testsession: raise IOError, "the output directory already exists:\n%s"%dirout
			else: shutil.rmtree(dirout)
		os.mkdir(dirout)
	nconflictingnodes = 0
	ncoveragehelpdecide = 0
	nblocksizehelpdecide = 0
	ncohabitexclusive = 0
	
	reftree = tree2.ReferenceTree(fic=nfreftree) #, branch_lengths=False)
	reftree.complete_internal_labels(prefix = 'N')
	reftree.complete_node_ids()
	
	dbcon, dbcur = rec_to_db.dbconnect(**kw)
	if useTempTables:
		# must create temporary tables like ones in the current db (only to store the new speciation events that may be turned into duplication events)
		for table in tables:
			rec_to_db.createTempTable(dbcur, table, truncateOnCommit=True, silent=silent)

	lnfgenetrees = os.listdir(dirgenetrees)
	treeext = lnfgenetrees[0].rsplit('.', 1)[1]
	lfam = [nfgenetree.rsplit('.', 1)[0] for nfgenetree in lnfgenetrees]
	if appendmode:
		lnfalreadycomputedtrees = os.listdir("%s/normed_phyloxml_trees"%dirout)
		lalreadycomputedfams = [nfalreadycomputedtree.rsplit('.', 1)[0] for nfalreadycomputedtree in lnfalreadycomputedtrees]
		lfam = list(set(lfam) - set(lalreadycomputedfams))
		lenlfam = len(lfam)
		print lenlfam, "families remain to be computed:"
		if lenlfam > 25: print ', '.join(lfam[:25]), "... [+ %d others]"%(lenlfam-25)
		else: print ', '.join(lfam)
	
	nfam = 0
	if not silent: print "- - - - - - -\n"
	for fam in lfam:
		nfam += 1
		# characterize gene tree
		nfgenetree = '%s.%s'%(fam, treeext)
		genetree = tree2.GeneTree(fic="%s/%s"%(dirgenetrees, nfgenetree))
		genetree.complete_internal_labels()
		genetree.complete_node_ids()
		dgt = rec_to_db.getGeneTreeInfo(fam, '*', dbcur)
		recgiid = int(dgt['rec_gi_id'])
		
		deventid_devent = {}
		dconflictreplicates_nodes = {}
		
		if not silent: print "# %d family:"%nfam, fam, '\n'
		else: print nfam, fam
		for node in genetree:
			# explore the events at the node
			startnodeid = node.nodeid()
			nrnodeid = rec_to_db.getNRNode(fam, startnodeid, dbcur)
			if not silent: print node.label(), "nodeid:", startnodeid, "nrnodeid:", nrnodeid
			tdevent = rec_to_db.getTreeEventsInfo(fam, startnodeid, ['event_id', 'event_type'], dbcur)
			
			deventid_blocksize = {}
			deventid_coverage = {}
			# evaluate the total coverage of the node by Prunier replicates
			# to spot nodes not determined as 'transfer' in all replicates -> the rest of replicates tell 'speciation'
			totalnodecoverage = set(rec_to_db.getNodeObservingReplicates(nrnodeid, dbcur))
			if not silent: print "totalnodecoverage (%d):"%len(totalnodecoverage), totalnodecoverage
			totaltransfercoverage = set([])
			possibleeventtypes = set([])
			
			for devent in tdevent:
				if not silent: print "\tevent:", devent
				eventid = int(devent['event_id'])
				eventtype = devent['event_type']
				possibleeventtypes |= set([eventtype])
				deventid_devent[eventid] = devent
				# characterize linked block event size
				tblocksizecounts = rec_to_db.getEventLinkedGeneCounts(eventid, dbcur)
				deventid_blocksize[eventid] = (mean(tblocksizecounts), max(tblocksizecounts))
				if not silent: print "\t\tblocksizes: mean:", mean(tblocksizecounts), "max:", max(tblocksizecounts)
				# characterize event coverage by Prunier replicates
				if eventtype=='T':
					trsupreps = rec_to_db.getTransferSupportingReplicates(eventid, dbcur)
					totaltransfercoverage |= set(trsupreps)
				else:
					trsupreps = []
				if not silent: print "\t\ttrsupreps (%d):"%len(trsupreps), trsupreps
				deventid_coverage[eventid] = trsupreps
						
			if len(totaltransfercoverage) > 0:
				# node has transfer event(s) recorded
				speciationcoverage = totalnodecoverage - totaltransfercoverage
				if (len(speciationcoverage)>0) and (possibleeventtypes==set(['T'])):
					# the node is not considered a transfer in every Prunier test, but no alternative event entry exists in the database for the moment
					# create an alternative event of speciation (temporary, may be dropped to replace by a duplication in call to integrateReconciliations())
					eventid = rec_to_db.insertNewSerialVal('phylogeny.event', [recgiid, 'S', None, startnodeid, None], dbcur, useTempTables=useTempTables)
					node.set_speciation(reftree=reftree)
					refloc = node.speciation()[0]
					for loc in node.speciation()[1]:
						epsn = reftree[loc]
						epsnid = epsn.nodeid()
						refnodebool = (loc==refloc)
						rec_to_db.executeSQLinsert(dbcur, 'phylogeny.event_possible_species_node', [eventid, epsnid, 'location', refnodebool], tmptable=useTempTables, silent=silent)
					devent = {'event_id':eventid, 'event_type':'S'}
					tdevent += (devent,)
					deventid_coverage[eventid] = speciationcoverage
					deventid_blocksize[eventid] = (1, 1)
					deventid_devent[eventid] = devent
					if not silent:
						print "\tevent:", devent
						print "\t\tblocksizes: mean:", 1, "max:", 1
						print "\t\tspeciationcoverage (%d):"%len(speciationcoverage), speciationcoverage
				
					
			## choose the event
			refdevent = None
			if len(tdevent) > 1:
				# must choose one event
				nconflictingnodes += 1
				meanblocksize = 0
				maxblocksize = 0
				maxcoverage = 0
				refeventtype = None
				decisions = (False, False, False)
				for devent in tdevent:
					eventid = devent['event_id']
					eventtype = devent['event_type']
					blocksize = deventid_blocksize[eventid]
					coverage = len(deventid_coverage[eventid])
					
					maxblocksizedecide = (blocksize[1] > maxblocksize)
					maxblocksizeequal = (blocksize[1] == maxblocksize)
					meanblocksizedecide = (blocksize[0] > meanblocksize)
					meanblocksizeequal = (blocksize[0] == meanblocksize)
					coveragedecide = ((coverage >= maxcoverage) and (coverage != 0))
					coverageequalnull = ((coverage == maxcoverage) and (coverage == 0))
					speciationovertransfer = ((eventtype in ['S', 'D']) and (refeventtype == 'T'))
					
					if maxblocksizedecide or (maxblocksizeequal and (meanblocksizedecide or (meanblocksizeequal and (coveragedecide or (coverageequalnull and speciationovertransfer))))):
						refdevent = devent
						refeventtype = eventtype
						maxblocksize = blocksize[1]
						meanblocksize = blocksize[0]
						maxcoverage = coverage
						decisions = (maxblocksizedecide, meanblocksizedecide, coveragedecide)
					
				# on what basis was last choice made?
				maxblocksizedecide, meanblocksizedecide, coveragedecide = decisions
				if maxblocksizedecide or meanblocksizedecide:
					nblocksizehelpdecide += 1
				if coveragedecide:
					ncoveragehelpdecide += 1
			else:
				# trivial choice, 
				refdevent = devent	
						
						
			refeventtype = rec_to_db.dabrev_evt[refdevent['event_type']]
			refeventid = int(refdevent['event_id'])
			if not silent: print "\t# refdevent", refdevent
			refeventloc = rec_to_db.getEventLocation(refeventid, refeventtype, dbcur)
			node.set_anyevent(refeventtype, refeventloc, refeventid)
			if not silent: print "\t\t", node.getdicevent()
			
			## map conflicting Prunier replicates 
			refcov = deventid_coverage[refeventid]
			if not silent: print "\t\trefcov", refcov
			for alteventid in deventid_coverage:
				if alteventid!=refeventid:
					# events that are not the reference event
					altdevent = deventid_devent[alteventid]
					if not silent: print "\tx altdevent", altdevent
					alteventtype = rec_to_db.dabrev_evt[altdevent['event_type']]
					if refeventtype == alteventtype:
						# event of the same type may be different but not in conflict
						alteventloc = rec_to_db.getEventLocation(alteventid, alteventtype, dbcur)
						if refeventtype=='transfer':
							refrecdonset = (refeventloc[1], refeventloc[3])
							altrecdonset = (alteventloc[1], alteventloc[3])
						else:
							refrecdonset = (refeventloc[1], [])
							altrecdonset = (alteventloc[1], [])
						if compatibleRecDonSet(refrecdonset, altrecdonset):
							# events are different but compatible = not in conflict
							continue # the for alteventid loop
					
					altcov = deventid_coverage[alteventid]
					if not silent: print "\t\taltcov", altcov
					for refreplicate in refcov:
						for altreplicate in altcov:
							tconf = (refreplicate, altreplicate)	# ordered tuple: (majority replicate, minority replicate)
							dconflictreplicates_nodes.setdefault(tconf, []).append(startnodeid)
		
		# search for exclusive replicate pairs that may have the reference event chosen diferently in different nodes
		nodescohabitexclusive = set([])
		spotted = []
		nprint = 0
		for reftconf in dconflictreplicates_nodes:
			if reftconf in spotted:
				continue
			refreplicate, altreplicate = reftconf
			revtconf = (altreplicate, refreplicate)
			if revtconf in dconflictreplicates_nodes:
				if nprint==0: print "Exclusive reconciliations by Prunier replicates cohabit in %s tree:"%(fam)
				if nprint<10: print "replicates %s and %s:"%(altreplicate, refreplicate)
				for tconf in (reftconf, revtconf):
					nodeconf = set(dconflictreplicates_nodes[tconf])
					nodescohabitexclusive |= nodeconf
					if nprint<10: print "\tchosen events supported by replicate %s:\n%s"%(tconf[0], '\n'.join([ "\t%s: %s"%(genetree.idgetlab(nodeid), str(genetree.idgetnode(nodeid).event())) for nodeid in nodeconf ]) )
				spotted.append(revtconf)
				nprint += 1
		ncen = len(nodescohabitexclusive)
		if nprint>=10: print "[%d more conflicting replicate pairs]"%(nprint-10)
		if nprint>0: print "%d nodes bear events that are exclusive with other node events"%ncen
		ncohabitexclusive += ncen
			
		if not silent: print "\nCompletion of reconciliation of gene tree knowing transfer events\n"
		rec_to_db.integrateReconciliations(fam, genetree, dirout, dbcur, reftree, minbs=minbs, reccolid=reccolid, recgiid=recgiid, checkEvents=True, inferTPMStransfers=True, silent=silent, useTempTables=useTempTables*2)	# useTempTables=2 : check non-temporary and temporary tables

					
		if useTempTables:
			if not silent: print "\nExport of SQL database to dump files"
			nftmpdumpfile = "%s/tmp_dump.tab"%cwd
			for table in tables:
				temptable = table.split('.')[1]
				rec_to_db.dumpTableToFile(dirout, temptable, dbcur, delim='\t', append=True, nftmpdumpfile=nftmpdumpfile, silent=silent)
			
		if testsession:
			dbcon.rollback()
		else:
			dbcon.commit()
			
		if not silent: print "- - - - - - -\n"
			
	print "-- -- -- -- -- -- --"
	print "\n\nfound", nconflictingnodes, "nodes where conflicting events where mapped"
	print "coverage helped to decide in", ncoveragehelpdecide, "cases"
	print "block size helped to decide in", nblocksizehelpdecide, "cases"
	print ncohabitexclusive, "nodes bear events that are exclusive with other node events"
	

	
def usage():
	s = "python refine_transfer_annotation.py [options] path/dirgenetrees path/reftree path/dirout\n"
	s+= "-d str\n\tPostgreSQL database guess string\n"
	s+= "-D str\n\tPostgreSQL database name\n"
	s+= "-H str\n\tPostgreSQL host server\n"
	s+= "-U str\n\tPostgreSQL user\n"
	s+= "-p str\n\tPostgreSQL user password\n"
	s+= "-b float\n\tminimum bootstrap threshold for event consideration\n\n"
	s+= "-r date\n\tidentifier of recociliation collection in the database, refrence to the coherent set of event recorded in output trees\n"
	s+= "-a\tappend mode: do not erase result directory previous to running the program\nand do not compute reconicliation when a file is already present in the output tree directory\n"
	s+= "-v\tverbose mode\n"
	s+= "-t\ttest mode\n"
	return s		

if __name__ == '__main__':
	options, args = getopt.getopt(sys.argv[1:], 'D:H:U:p:d:b:r:vta')
	if len(args) < 3:
		print usage()
		sys.exit(2)		
	dopt = dict(options)
	dbname = dopt.get('-D')
	dbuser = dopt.get('-U')
	dbhost = dopt.get('-H')
	dbpwd = dopt.get('-p')
	dbclue = dopt.get('-d')
	minbs = float(dopt.get('-b', 0.9))
	reccolid = dopt.get("-r", time.strftime("%Y-%m-%d"))
	silent = ('-v' not in dopt)
	testsession = ('-t' in dopt)
	appendmode = ('-a' in dopt)
	if None in [dbhost, dbuser, dbname] and not dbclue: dbclue = 'pbil-sgbd'
	dirgenetrees = args[0]
	nfreftree = args[1]
	dirout = args[2]
	if len(args)>3:
		print "unused arguments:", args[3:]		
	print "dirgenetrees = '%s'"%dirgenetrees
	print "dirout = '%s'"%dirout
	print "minbs =", minbs
	print "testsession =", testsession
	
	main(dirgenetrees, nfreftree, dirout, reccolid, minbs=minbs, dbhost=dbhost, dbuser=dbuser, dbname=dbname, dbpwd=dbpwd, dbclue=dbclue, silent=silent, testsession=testsession, appendmode=appendmode)

