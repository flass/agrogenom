#!/usr/bin/python
# -*- coding: utf8 -*-

import sys, getopt
import os, shutil
import copy
import tree2
import rec_to_db
import cPickle
cwd = os.getcwd()

silent=True
excludeddatasetspe = ['LIBAP', 'LIBSC']
# key specificity contrasting clades for defining "closely related" clades: 
constrastingclades = ['N15', 'N10','N5','N4']

annotfields = ["subfamily", "organism", "family", "hogenom_gene_id", "gene_name", "locus_tag", "new_locus_tag", "chromosome", "begin", "end", "description", \
"gene_tree_event", "block_event", "nb_related_strains_possessing_the_subfamily", "nb_remote_strains_possessing_the_subfamily", "related_clades_possessing_the_subfamily_(codes)", "remote_clades_possessing_the_subfamily_(codes)", "related_clades_possessing_the_subfamily_(names)", "remote_clades_possessing_the_subfamily_(names)"]

def write_event_leaves(subfam, fout, lglab, dgene_annot, dcode_names, dbcur, nbrelatedstrains='', nbotherstrains='', relatedclades=[], otherclades=[], event=None):
	if event:
		if event['eventtype']=='transfer': eventstr = "%s "%event['eventlocation'][2]
		else: eventstr = ""
		eventstr += "(%s)-> %s"%(event['eventtype'], event['eventlocation'][0])
		nbevtinblock = 0
		ancblockid = None
		if event['eventid']:
			dbcur.execute("SELECT anc_block_id FROM phylogeny.event WHERE event_id=%d;"%event['eventid'])
			ancblockid = dbcur.fetchone()[0]
		levtstr = [eventstr, str(ancblockid)]
	else:
		levtstr = [""]*2
		
	nbstr = [str(nbrelatedstrains), str(nbotherstrains)]
	codecla = [', '.join(list(relatedclades)), ', '.join(list(otherclades))]
	namecla = [', '.join([dcode_names[cla] for cla in relatedclades]), ', '.join([dcode_names[cla] for cla in otherclades])]
	s = '\n'.join(['\t'.join([subfam]+[str(annot) for annot in dgene_annot[glab]]+levtstr+nbstr+codecla+namecla) for glab in lglab])
	fout.write(s+'\n')

def main(nflnffamfasta, dirtreepickles, nfreftree, diroutancestral, infer="Dollo", appendmode=True, \
    nfinpresabsmat=None, dirfamhystorypickles=None, dirsubfamostortpickles=None, \
    outcladegainlosstable=True, outcladespecificgainlosstable=True, outsubfamtrees=True, outsubfamtables=True, \
    outfamphyloprofilemat=True, outsubfamphyloprofilemat=True, outgenomesynthesis=True, outfamnhx=True, outfampickle=True, **kw):
	
	silent = kw.get('silent', True)
	useTempTables = kw.get('useTempTables', 2)
	CountOptions = kw.get('CountOptions', "")
	
	flnffamfasta = open(nflnffamfasta,'r')
	lnffamfasta = flnffamfasta.readlines()
	flnffamfasta.close()
	reftree = tree2.ReferenceTree(fic=nfreftree)
	reftree.complete_node_ids()

	if infer:
		if not CountOptions: diroutscenarii = '%s/%s'%(diroutancestral, infer)
		else: diroutscenarii = '%s/%s.%s'%(diroutancestral, infer, CountOptions.strip('-').replace(' ', '_').replace('.', ''))
	else:
		diroutscenarii = diroutancestral
	print "output directory:", diroutscenarii
	diroutstate = '%s/family_states'%diroutscenarii
	diroutevents = '%s/family_events'%diroutscenarii
	diroutpickles = '%s/family_history_pickles'%diroutscenarii
	dirouttables = '%s/annottables'%diroutscenarii
	diroutsubtrees = '%s/ortho_subtrees'%diroutscenarii
	#~ diroutsubleaves = '%s/ortho_subfams'%diroutscenarii
	diroutphyloprofiles = '%s/phylogenetic_profiles'%diroutscenarii
	diroutmodtrees = '%s/modified_trees'%diroutscenarii
	
	rootdirouts = [diroutancestral, diroutscenarii]
	dirouts = rootdirouts
	if outfamnhx: dirouts += [diroutstate, diroutevents] 
	if outfampickle: dirouts.append(diroutpickles)
	if outcladegainlosstable or outcladespecificgainlosstable: dirouts.append(dirouttables)
	if outsubfamtrees: dirouts.append(diroutsubtrees)
	#~ if outsubfamtables: dirouts.append(diroutsubleaves)
	if outfamphyloprofilemat or outsubfamphyloprofilemat: dirouts.append(diroutphyloprofiles)
	if infer: dirouts.append(diroutmodtrees)	# ancestral content inferences may imply adding transfers to reconciliations
	
	for dirout in dirouts:
		if os.access(dirout, os.F_OK):
			if not appendmode and (not dirout in rootdirouts):
				# never erase root directories
				shutil.rmtree(dirout)
				os.mkdir(dirout)
		else:
			os.mkdir(dirout)
	
	# set AGROGENOM database connection
	dbcon, dbcur = rec_to_db.dbconnect(**kw)
	
	if infer:
		# possibility of new gene tree reconciliations 
		tables = ['phylogeny.event', 'phylogeny.event_possible_species_node', "phylogeny.represented_event"]
		if useTempTables:
			# must create temporary tables like ones in the current db (only to store the new speciation events that may be turned into duplication events)
			for table in tables:
				rec_to_db.createTempTable(dbcur, table, truncateOnCommit=True, silent=silent)

	refnodelabs = reftree.get_children_labels()
	refleaflabs = reftree.get_leaf_labels()
	if outfamphyloprofilemat: 
		# prepare family phylogenetic profiles/scenarii output matrices
		foutmatfamall = open('%s/fam_node_counts.mat'%diroutphyloprofiles, 'w')
		foutmatfamleaf = open('%s/fam_leaf_counts.mat'%diroutphyloprofiles, 'w')
		foutmatfamall.write('\t'.join(['family']+refnodelabs)+'\n')
		foutmatfamleaf.write('\t'.join(['family']+refleaflabs)+'\n')
	if outsubfamphyloprofilemat: 
		# prepare subfamily phylogenetic profiles/scenarii output matrices
		foutmatfamall = open('%s/fam_node_counts.mat'%diroutphyloprofiles, 'w')
		foutmatsubfamall = open('%s/subfam_node_counts.mat'%diroutphyloprofiles, 'w')
		foutmatfamleaf = open('%s/fam_leaf_counts.mat'%diroutphyloprofiles, 'w')
		foutmatsubfamleaf = open('%s/subfam_leaf_counts.mat'%diroutphyloprofiles, 'w')
		foutmatsubfamall.write('\t'.join(['family']+refnodelabs)+'\n')
		foutmatsubfamleaf.write('\t'.join(['family']+refleaflabs)+'\n')

	if outcladegainlosstable or outcladespecificgainlosstable:
	# prepare detrailed gained/loss gene annoation tabular output
		drefnodesout = {}
		for refnode in reftree:
			dnodefiles = {}
			levt = []
			if outcladegainlosstable: levt += ['gain', 'loss']
			if outcladespecificgainlosstable: levt += ['specific_presence', 'specific_absence']
			for evt in levt:
				if appendmode: m = 'a'
				else: m = 'w'
				dnodefiles[evt] = open("%s/%s.%s"%(dirouttables, refnode.label(), evt), m)
				dnodefiles[evt].write( '\t'.join( annotfields )+'\n')
			drefnodesout[refnode.label()] = dnodefiles
		# fetch annotations for all genes in db
		dgene_annot = rec_to_db.get_all_gene_annotation(dbcur)	
		# fetch scientific names of taxa
		dcode_names = rec_to_db.get_scientific_name_dict(reftree, dbcur)

	if outgenomesynthesis:
		genomeHistory = copy.deepcopy(reftree)
		
	if outsubfamtables:
		foutsubfamleaves = open('%s/ortho_subfam_leaves.tab'%diroutscenarii, 'w')
		foutsubfamevents = open('%s/ortho_subfam_events.tab'%diroutscenarii, 'w')


	nfam = 0
	if dirfamhystorypickles:
		print 'load pre-computed family histories from:\t%s\n'%dirfamhystorypickles
	else:
		print 'load families inventory from:\t%s\n'%nflnffamfasta
		print 'load (if exist) families trees from:\t%s\n'%dirtreepickles
		
	for line in lnffamfasta:
		nffasta = line.rstrip('\n')
		nfam += 1
		fam = nffasta.rsplit('/',1)[1].rsplit('.', 1)[0]
		
		nfgenetree = "%s/%s.pickle"%(dirtreepickles, fam)
		if not os.access(nfgenetree, os.F_OK):
			indata = 'fasta'
			print 'fasta #%d:\t%s'%(nfam,fam)
			ffasta = open(nffasta, 'r')
			dspe_leaves = {}
			lleaves = []
			for line in ffasta:
				if line.startswith('>'):
					leaf = line.strip('>\n').split()[0]
					spe = leaf.split('_')[0]
					if not spe in excludeddatasetspe:
						dspe_leaves.setdefault(spe, []).append(leaf)
						lleaves.append(leaf)
			genetree = tree2.Node(lleaves=lleaves)
		else:
			indata = 'tree'
			print 'tree #%d:\t%s'%(nfam,fam)
			genetree = tree2.load_pickle(nfgenetree)			
			
		if dirfamhystorypickles:
			# load previously computed history
			# NB: cannot write family-specific informations such as detailled gain-loss tables
			nfhistorypickle = '%s/%s.pickle'%(dirfamhystorypickles, fam)
			if os.access(nfhistorypickle, os.F_OK):
				print 'pickle #%d:\t%s'%(nfam,fam)
				famHistory = tree2.load_pickle(nfhistorypickle)
			else:
				continue
			if dirsubfamostortpickles:
				nfsubfamostortpickles = '%s/%s.ostort.pickle'%(dirsubfamostortpickles, fam)
				if os.access(nfsubfamostortpickles, os.F_OK):
					fostort = open(nfsubfamostortpickles,'r')
					[orthoSubtrees, orthoReftrees] = cPickle.load(fostort)
					fostort.close()
		else:			
				
			if indata =='fasta':
				orthoSubtrees = {1:genetree}
				famHistory = copy.deepcopy(reftree)
				famHistory.addCountInTree(dspe_leaves)
				if len(lleaves)==0:
					print "family specific to species in excluded dataset", excludeddatasetspe, "; skip."
					continue
				#~ elif len(lleaves) in range(1,4):
				else:
					# just want families with at least 1 sequence (that are not avoided species) and with at most 3 sequences, because trees are not computed for them ; 
					# can reconstrcut states with DolloPars or with Count's AssymetricWagner
					if infer=="Dollo":
						famHistory.DolloPars()
					elif infer.startswith("Count"):
						famHistory.Count(prog=infer.split('.')[1], options=CountOptions, **kw)
				# else:
					# tree was missed for another reason that a small fam (PhyML or Prunier SegFaults...); cannot reconstruct states with DolloPars without reconciliatiopn information: ignore this family	
				orthoReftrees = {1:famHistory}
				
			else:
				lleaves = genetree.get_leaf_labels()
				lspe = []
				for leaf in lleaves:
					spe = leaf.split('_')[0]
					if spe not in lspe:
						lspe.append(spe)	
						
				# Infers presence/absence history and orthologous groups linked to events
				if infer:
					famHistory, orthoReftrees, orthoSubtrees, modifiedgt = reftree.ancestralStatesFromGeneTree(genetree, method=infer, **kw)
					if modifiedgt:
						if not silent: print "gene tree  was modified"
						rec_to_db.integrateReconciliations(fam, modifiedgt, diroutmodtrees, dbcur, reftree, complete=True, inferTPMStransfers=False, duplications=True, speciations=True, checkEvents=True, **kw)
						# update events in ost/ort dictionaries
						leventnodes = orthoReftrees.keys()
						for evtnode in leventnodes:
							modevtnode = modifiedgt[evtnode.label()]
							if evtnode.eventid() != modevtnode.eventid():
								if not silent: print "update event at %s in ost/ort dictionaries:\n%s -> %s"%(evtnode.label(), evtnode.event(), modevtnode.event())
								ost = orthoSubtrees[evtnode]
								ort = orthoReftrees[evtnode]
								if modevtnode not in orthoSubtrees: m = modevtnode
								else: m = copy.copy(modevtnode)
								orthoSubtrees[m] = ost
								orthoReftrees[m] = ort
								del orthoSubtrees[evtnode]
								del orthoReftrees[evtnode]
				else:
					raise ValueError, "Reftree.ancestralStatesFromMatrix() not implemented yet, need to infer ancestral states"
					
				
				devt_leaves = {}
				for evtnode in orthoReftrees:
					ost = orthoSubtrees[evtnode]
					devt_leaves[evtnode] = ost.get_leaf_labels()
			
			if not silent: print "\n# # # # # # # # # # # # # # #\n"

		famHistory.checkGainLossPresenceCoherence()

		if outfampickle:
			# Dump gene family history
			tree2.dump_pickle(famHistory, '%s/%s.pickle'%(diroutpickles, fam))
		if outfamnhx:
			# Summarize gene family history
			famHistory.write_newick('%s/%s.nhx'%(diroutstate, fam), comment='presence', mode='write')
			famHistory.write_newick('%s/%s.nhx'%(diroutevents, fam), comment='events', mode='write')
		if outgenomesynthesis:
			# Sum it to the whole history of genomes
			genomeHistory += famHistory
			genomeHistory.checkGainLossPresenceCoherence()				
		
		# get detailled events from orthologous groups and fetch annotations of concerned genes
		if (not dirfamhystorypickles) or (dirfamhystorypickles and dirsubfamostortpickles):
			sostleaves = set()
			leventnodes = orthoSubtrees.keys()
			if isinstance(leventnodes[0], tree2.Node):
				leventnodes = genetree.sort(leventnodes, order=2)
			for e, evtnode in enumerate(leventnodes):				
				if isinstance(evtnode, tree2.Node):
					subfamid = str(e)
				else:
					subfamid = 's'
				subfam = "%s.%s"%(fam, subfamid)
				ort = orthoReftrees[evtnode]
				ost = orthoSubtrees[evtnode]
				lglab = ost.get_leaf_labels()
				# check that leaves are not assigned to several ost
				if sostleaves & set(lglab):
					print subfam
					print ost
					print "sostleaves", sostleaves
					raise IndexError, "leaves alredy assigned to a previous subtree"
				else:
					#~ print sostleaves, "|=", lglab
					sostleaves |= set(lglab)
				
				ortpresents = ort.presents(returnLabels=True)
				mrcapresents = ort.map_to_node(ortpresents)
				# clades where it is present/absent everywhere
				specificpresclades = ort.whichClade(ortpresents)
				specificabsclades = []
				# clades where it is present almost everywhere
				aproxclades = ort.whichClade(ortpresents, holes=0.2)
				#~ print "presclades", presclades
				#~ print "aproxclades", aproxclades
				#~ print "\nspecificpresclades", specificpresclades
				lnodevt = []
				
				if not ((subfamid == 's') and (ost.nb_leaves() > 3)):
					# filter big families with no trees
					#~ print subfam
					for refnode in ort:
						#~ print "refnode", refnode.label()
						restrost = ost.restrictToLeaves(refnode.get_leaf_labels(), useSpeDict=True, force=True)
						if restrost: glabs = restrost.get_leaf_labels()
						else: glabs = []
						devts = refnode.getEvents()
						# the other clades in which the gene is present
						nbotherstrains = len(set(ortpresents) - set(refnode.get_leaf_labels()))
						nbrelatedstrains = 0
						# the other clades in which the gene is present
						otherclades = set(aproxclades) - set([refnode.label()])
						# the closely related clades in which the gene is present
						relatedclades = set([])
						# the ancestor under which other clades are considered close relatives
						contrastclade = None
						for contcla in constrastingclades:
							if refnode.is_child(ort[contcla]):
								contrastclade = contcla
								break
						if contrastclade:
							for othcla in otherclades:
								if ort[othcla].is_child(ort[contrastclade]):
									relatedclades.add(othcla)
									nbothclastr = len(set(ort[othcla].get_leaf_labels()) & set(ortpresents))
									nbrelatedstrains += nbothclastr
						otherclades -= relatedclades
						nbotherstrains -= nbrelatedstrains
							
						#~ print "otherclades", otherclades
						if outcladegainlosstable or outcladespecificgainlosstable:
						# write gene gain/loss (when specifically present/absent) by clade with detailed annotation in tables
							if devts['gain']:
								# associate a gene tree event to each gain
								for ng in range(devts['gain']):
									#~ print refnode.label(), "gain"
									ev = None
									#~ print "here 0", refnode.label()
									if isinstance(evtnode, tree2.Node):
										ev = evtnode.getdicevent()
										if evtnode.eventloc()[1]==refnode.label() and (not (evtnode in lnodevt)):
											# referenece subtree event
											#~ glabs = ost.restrictToLeaves(refnode.get_leaf_labels()).get_leaf_labels()
											#~ ev = evtnode.getdicevent()
											lnodevt.append(evtnode)
										else:
											# look in events within the subtree that could match
											for nost in ost:
												#~ if nost.eventtype()=='transfer': print "nost.event()", nost.event()
												if (not (nost in lnodevt)) and nost.eventtype()=='transfer' and nost.eventloc()[1]==refnode.label():
													glabs = nost.restrictToLeaves(refnode.get_leaf_labels(), useSpeDict=True, force=True).get_leaf_labels()
													ev = nost.getdicevent()
													lnodevt.append(nost)
											
									if outcladegainlosstable:
										write_event_leaves(subfam, drefnodesout[refnode.label()]['gain'], glabs, dgene_annot, dcode_names, dbcur, nbrelatedstrains=nbrelatedstrains, nbotherstrains=nbotherstrains, relatedclades=relatedclades, otherclades=otherclades, event=ev)
									if outcladespecificgainlosstable and refnode.homogeneous_state_in_clade()==True:
										for apcla in aproxclades:
											# test that refnode is not included in a higher clade where it is present almost everywhere
											if refnode.is_child(ort[apcla]):
												#~ print "present in", refnode.label(), "and in other child clades of", apcla
												# present in the clade above but in the same copy number?
												for spprcl in specificpresclades:
													if ort[spprcl].is_child(ort[apcla]) and refnode.is_child(ort[spprcl]):
														# subsitute true homogeneously occuring clade for contrasting
															contrapcla = spprcl
															break # the sppprcl loop
												else:
													contrapcla=apcla
												foreback = refnode.clade_specific_contrast(contrast=contrapcla, counts=True, returnforeback=True)
												if (nost.nb_leaves() > 1) and foreback and (foreback[0] > foreback[1]) and (foreback[0] <= 2):
													# the subfamily exists in one more copy in this clade
													mulleaves = refnode.get_leaf_labels()
													mulancs = nost.whichMulticopySubtrees(multi=mulleaves, returnLabels=False)
													oneanc = nost.map_to_node(ort[contrapcla].get_leaf_labels(), useSpeDict=True)
													for mulanc in mulancs:
														if not mulanc.is_child(oneanc):
															# finds the gene surnumerary copy subtree which is not included in the larger subtree of original copy
															glabs = list(set(glabs) & set(mulanc.get_leaf_labels()))
															if glabs: break # the for mulanc loop
													else:
														for mulanc in mulancs:
															# can be included in the larger subtree if duplicative transfer
															if refnode.label() in [emul['eventlocation'][0] for emul in oneanc.getEvents(eventtype='transfer', lineage=mulanc.label())]:
																glabs = list(set(glabs) & set(mulanc.get_leaf_labels()))
																if glabs: break # the for mulanc loop
														else:
															break # the for apcla loop
												else:
													break # the for apcla loop
										else:
											#~ print "ev", ev
											#~ print "here 1", refnode.label()
											write_event_leaves(subfam, drefnodesout[refnode.label()]['specific_presence'], glabs, dgene_annot, dcode_names, dbcur, nbrelatedstrains=nbrelatedstrains, nbotherstrains=nbotherstrains, relatedclades=relatedclades, otherclades=otherclades, event=ev)
											# to not write twice the gain
											if refnode.label() in specificpresclades:
												specificpresclades.remove(refnode.label())
												#~ print "specificpresclades", specificpresclades
								
							elif devts['loss']:
								if outcladegainlosstable:
									write_event_leaves(subfam, drefnodesout[refnode.label()]['loss'], lglab, dgene_annot, dcode_names, dbcur)
								if outcladespecificgainlosstable and refnode.homogeneous_state_in_clade()==False:
									write_event_leaves(subfam, drefnodesout[refnode.label()]['specific_absence'], lglab, dgene_annot, dcode_names, dbcur, nbrelatedstrains=nbrelatedstrains, nbotherstrains=nbotherstrains, relatedclades=relatedclades, otherclades=otherclades)
									specificabsclades.append(refnode.label())
									#~ print "here 3", refnode.label()
									
						if outcladespecificgainlosstable: 
							if refnode.label() in specificpresclades:
								for apcla in aproxclades:
									# test that refnode is not included in a higher clade where it is present almost everywhere
									if refnode.is_child(ort[apcla]):
										#~ print "present in", refnode.label(), "and in other child clades of", apcla
										break
								else:
									#~ print "ev", ev
									#~ print "here 2", refnode.label()
									write_event_leaves(subfam, drefnodesout[refnode.label()]['specific_presence'], glabs, dgene_annot, dcode_names, dbcur, nbrelatedstrains=nbrelatedstrains, nbotherstrains=nbotherstrains, relatedclades=relatedclades, otherclades=otherclades)
									#~ print refnode.label(), "specific pres"
									specificpresclades.remove(refnode.label())
									#~ print "specificpresclades", specificpresclades
							elif (not (refnode.label() in specificabsclades)) and (refnode.clade_specific_contrast(contrast='allpres', force=True)==False):
								write_event_leaves(subfam, drefnodesout[refnode.label()]['specific_absence'], lglab, dgene_annot, dcode_names, dbcur, nbrelatedstrains=nbrelatedstrains, nbotherstrains=nbotherstrains, relatedclades=relatedclades, otherclades=otherclades)
								specificabsclades.append(refnode.label())
								#~ print "here 4", refnode.label()
							
				if outsubfamphyloprofilemat:			
					# write subfamily phylogenetic profile in matrices
					for fmatout, reflabs, leavesOnly in [(foutmatsubfamall, refnodelabs, False), (foutmatsubfamleaf, refleaflabs, True)]:
						p = ort.writePhyloProfileMatrix(fmatout, subfam, reflabs, leavesOnly)	
						if leavesOnly & (sum(p.values()) != len(ost.get_leaves())):
							#~ print 'p = ', p
							#~ print 'lleaves = ', lglab
							#~ print 'ostnwk = "%s"'%ost.newick()
							for leaf in lglab:
								spe = leaf.split('_')[0]
								p[spe] -= 1
							dspeleaves = ost.dictSpeciesToLeafLabels()
							for spe in p:
								if p[spe]!=0:
									print spe,'; diff count = ', p[spe], '; lleaves = ', dspeleaves[spe]
							raise IndexError, "subfamily %s (under node %s, event %s) phylogenetic profile is false"%(subfam, evtnode.label(), str(evtnode.event()))
							
				if outsubfamtables:
					foutsubfamleaves.writelines(['%s\t%s\n'%(glab, subfam) for glab in lglab])
					if isinstance(ost, tree2.GeneTree): foutsubfamevents.writelines(['%d\t%s\t%s\n'%(ostnode.eventid(), subfam, 'F') for ostnode in ost])
					if isinstance(evtnode, tree2.GeneTree) and evtnode.eventid(): foutsubfamevents.write('%d\t%s\t%s\n'%(evtnode.eventid(), subfam, 'T')) # origin=True

			if outsubfamtrees:
				fostort = open("%s/%s.ostort.pickle"%(diroutsubtrees, fam), 'w')
				cPickle.dump([orthoSubtrees, orthoReftrees], fostort)
				fostort.close()
				
				
		if outfamphyloprofilemat:
			# write family phylogenetic profile in matrices
			for fmatout, reflabs, leavesOnly in [(foutmatfamall, refnodelabs, False), (foutmatfamleaf, refleaflabs, True)]:
				p = famHistory.writePhyloProfileMatrix(fmatout, fam=fam, reflabs=reflabs, leavesOnly=leavesOnly)
				if not dirfamhystorypickles:
					if leavesOnly and (sum(p.values()) != len(genetree.get_leaves())):
						#~ print 'p = ', p
						lleaves = genetree.get_leaf_labels()
						#~ print 'lleaves = ', lleaves
						sostleaves = set()
						for evtnode in orthoReftrees:
							ost = orthoSubtrees[evtnode]
							sostleaves |= set(ost.get_leaf_labels())
						#~ print 'sostleaves = ', sostleaves
						for leaf in lleaves:
							spe = leaf.split('_')[0]
							p[spe] -= 1
						dspeleaves = genetree.dictSpeciesToLeafLabels()
						for spe in p:
							if p[spe]!=0:
								print spe,'; diff count = ', p[spe], '; lleaves = ', dspeleaves[spe]
						raise IndexError, "family phylogenetic profile is false"
					
	if outgenomesynthesis:
		# save the whole history of genomes
		genomeHistory.write_newick('%s/genome_event_synthesis.nhx'%(diroutscenarii), comment='events', mode='write')
		genomeHistory.write_newick('%s/genome_states_synthesis.nhx'%(diroutscenarii), comment='presence', mode='write')
		genomeHistory.write_newick('%s/genome_synthesis.nhx'%(diroutscenarii), comment='gainlosscount', mode='write')
		tree2.dump_pickle(genomeHistory, '%s/genome_synthesis.pickle'%(diroutscenarii))
	if outcladegainlosstable:
		# close tabular result files
		for refnodelab in drefnodesout:
			for evt in drefnodesout[refnodelab]:
				fout = drefnodesout[refnodelab][evt]
				fout.close()			
	if outsubfamphyloprofilemat:
		# close count matrix files
		foutmatfamall.close()
		foutmatsubfamall.close()
		foutmatfamleaf.close()
		foutmatsubfamleaf.close()
		#~ sys.stdout.write('\n')
		
	if outsubfamtables:
		foutsubfamleaves.close()

	if infer and infer.startswith("Count"):
		if useTempTables:
			# export of SQL database to dump files
			nftmpdumpfile = "%s/tmp_dump.tab"%cwd
			for table in tables:
				temptable = table.split('.')[1]
				rec_to_db.dumpTableToFile(diroutmodtrees, temptable, dbcur, delim='\t', append=appendmode, nftmpdumpfile=nftmpdumpfile, silent=silent)
		
	# close database connection
	dbcon.close()

def tf(val):
	if isinstance(val, bool):
		return val
	elif isinstance(val, str):
		if val.upper() in ['TRUE', 'T', '1']:
			return True
		elif val.upper() in ['FALSE', 'F', '0']:
			return False
		else:
			raise ValueError, val
	else:
		raise ValueError, val
	
def usage():
	s = "python ancestral_content.py lnffamfasta dirtreepickles reftree diroutscenarii\n"
	s+= "Phylogenomic database connection options:\n"
	s+= "-d str\n\tPostgreSQL database guess string\n"
	s+= "-D str\n\tPostgreSQL database name\n"
	s+= "-H str\n\tPostgreSQL host server\n"
	s+= "-U str\n\tPostgreSQL user\n"
	s+= "-p str\n\tPostgreSQL user password\n"
	s+= "-b float\tminimum bootstrap threshold for reconciliation event consideration\n"
	s+= "-o str\n\tCount options (surrounded by quotes, ex: -o \"-gain 3 -max_paralogs 100\")\n"
	s+= "-i\t\tpopulate database with INSERT commands directly on true tables;\n\t\tdefault behaviour uses INSERT and COPY commands on temporary tables to export dumps\n"
	s+= "\nAncestral Content inference options:\n"
	s+= "--infer=method\n\tinfer ancestral states and gain/loss scenarii using:\n"
	s+= "\tDollo parsimony (method='Dollo')\n"
	s+= "\tAssymetric Wagner parsimony (method='Count.AssymetricWagner')\n"
	s+= "\tPosterior probabilities of DTL model (method='Count.Posteriors', need --model-rates to be set)\n"
	s+= "--tmp-dir=directory path\n\tpath to directory where Count writes temporary files\n"
	s+= "--count-path=directory path\n\tpath to the Count.jar file\n"
	s+= "--model-rates=file path\n\tpath to the DTL model rates file\n"
	s+= "\nInput options:\n"	
	s+= "--input-presence-matrix=file path\n\tancestral presence/absence state matrix (as exported from Count [Csuros & Miklos, 2009])\nNOT IMPLEMENTED YET\n"
	s+= "--input-family-history-pickles=directory path\n\tuse previously computed family history trees to generate ancestral content data\n"
	s+= "--input-subfamily-trees=directory path\n\tuse previously computed subfamily histories and gene subtrees to generate subfamily-specific data\n"
	s+= "\nOutput options:\n"
	s+= "--output-all\n\tcombine all following output options (can be overriden for specifying particular output options if there provided as false)\n"
	s+= "--output-clade-gainloss-tables=bool\n\twrite gene gain/loss tables by clade with detailed annotation (default false)\n"
	s+= "--output-clade-specific-gainloss-tables=bool\n\twrite gene gain/loss tables by clade when specifically present/absent in all species descendent, with detailed annotation (default false)\n"
	s+= "--output-subfamily-trees=bool\n\twrite subfamily trees (default false)\n"
	s+= "--output-subfamily-tables=bool\n\twrite into tables HOGENOM gene ids and event ids included in a subfamily (default false)\n"
	s+= "--output-family-phyloprofile-matrices=bool\n\twrite full family phylogenetic profiles in matrices (default false)\n"
	s+= "--output-subfamily-phyloprofile-matrices=bool\n\twrite subfamily phylogenetic profiles in matrices (default false)\n"
	s+= "--output-NHX-formated-family-histories=bool\n\twrite full family history summaries in two NHX-formated tree files, one for events, another for homolog counts (default false)\n"
	s+= "--output-pickled-family-histories=bool\n\tdump full family history in tree2 object pickles (default true)\n"
	s+= "--output-genome-synthesis=bool\n\tsummarize events and presence informations of all gene families in a genome history (default false)\n"
	s+= "-l\tlocal computing mode for execution of Count, sets to --tmp-dir and --count-path to the current directory\n"
	s+= "-a\tappend mode: do not erase result directory previous to running the program and concatenate table dumps to already existing files\n"
	s+= "-v\tverbose mode\n"
	s+= "-t\ttest mode\n"
	return s		

if __name__ == '__main__':
	options, args = getopt.getopt(sys.argv[1:], 'D:H:U:p:d:o:b:vtla', \
	["infer=", "input-presence-matrix=", "input-family-history-pickles=", "input-subfamily-trees=", "model-rates=", "tmp-dir=", 'count-path=', \
	'output-all', 'output-clade-gainloss-tables=', 'output-clade-specific-gainloss-tables=', 'output-subfamily-trees=', 'output-subfamily-tables=', 'output-family-phyloprofile-matrices=', 'output-subfamily-phyloprofile-matrices=', \
	'output-NHX-formated-family-histories=', 'output-pickled-family-histories=', 'output-genome-synthesis='])
	if len(args) < 3:
		print usage()
		sys.exit(2)		
	dopt = dict(options)
	dbname = dopt.get('-D')
	dbuser = dopt.get('-U')
	dbhost = dopt.get('-H')
	dbpwd = dopt.get('-p')
	dbclue = dopt.get('-d')
	CountOptions = dopt.get('-o', "")
	minbs = dopt.get('-b')
	if minbs: minbs = float(minbs)
	nfinpresabsmat = dopt.get('--input-presence-matrix')
	dirfamhystorypickles = dopt.get('--input-family-history-pickles')
	dirsubfamostortpickles = dopt.get('--input-subfamily-trees')
	allout = ('--output-all' in dopt)
	outcladegainlosstable = tf(dopt.get('--output-clade-gainloss-tables', allout))
	outcladespecificgainlosstable = tf(dopt.get('--output-clade-specific-gainloss-tables', allout))
	outsubfamtrees = tf(dopt.get('--output-subfamily-trees', allout))
	outsubfamtables = tf(dopt.get('--output-subfamily-tables', allout))
	outfamphyloprofilemat = tf(dopt.get('--output-family-phyloprofile-matrices', allout))
	outsubfamphyloprofilemat = tf(dopt.get('--output-subfamily-phyloprofile-matrices', allout))
	outgenomesynthesis = tf(dopt.get('--output-genome-synthesis', allout))
	outfamnhx = tf(dopt.get('--output-NHX-formated-family-histories', allout))
	outfampickle = tf(dopt.get('--output-pickled-family-histories', True))
	nfrates  = dopt.get('--model-rates', "")
	infer = dopt.get('--infer')
	silent = ('-v' not in dopt)
	testsession = ('-t' in dopt)
	appendmode = ('-a' in dopt)
	useTempTables = (not ('-i' in dopt))*2
	if '-l' in dopt:
		tmpdir = countdir = os.getcwd()
	else:
		tmpdir = dopt.get('--tmp-dir', "/panhome/lassalle/tmp")
		countdir = dopt.get('--count-path', "/panhome/lassalle/count")
	if None in [dbhost, dbuser, dbname] and not dbclue: dbclue = 'pbil-sgbd'
	nflnffamfasta = args[0]
	dirtreepickles = args[1]
	nfreftree = args[2]
	diroutancestral = args[3]
	if len(args)>4:
		print "unused arguments:", args[3:]		
	if nfinpresabsmat:
		print "use %s matrix for ancestral states"%nfinpresabsmat
		if infer:
			print "(no ancestral state reconstruction, ignore '--infer' option value)"
			infer = False
	elif infer:
		print "infer ancestral states with", infer
		if CountOptions: print "with options", CountOptions
		if nfrates: print "with rates from", nfrates
	call = ' '.join(sys.argv)
	if dbpwd: call = call.replace(dbpwd, '*****')
	print "# call:\npython", call
	main(nflnffamfasta, dirtreepickles, nfreftree, diroutancestral, \
	nfinpresabsmat=nfinpresabsmat, dirfamhystorypickles=dirfamhystorypickles, dirsubfamostortpickles=dirsubfamostortpickles, \
	outcladegainlosstable=outcladegainlosstable, outcladespecificgainlosstable=outcladespecificgainlosstable, \
	outsubfamtrees=outsubfamtrees, outsubfamtables=outsubfamtables, outfampickle=outfampickle, \
	outfamphyloprofilemat=outfamphyloprofilemat, outsubfamphyloprofilemat=outsubfamphyloprofilemat, \
	outgenomesynthesis=outgenomesynthesis, outfamnhx=outfamnhx, \
	infer=infer, CountOptions=CountOptions, nfrates=nfrates, \
	dbhost=dbhost, dbuser=dbuser, dbname=dbname, dbpwd=dbpwd, dbclue=dbclue, \
	silent=silent, minbs=minbs, tmpdir=tmpdir, countdir=countdir, useTempTables=useTempTables, appendmode=appendmode)
