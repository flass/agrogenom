#!/usr/bin/python
# -*- coding: utf8 -*-

import sys
import os
import copy
import tree2
import rec_to_db
cwd = os.getcwd()

def setToStr(s):
	l = list(s)
	l.sort()
	return "|".join(l)

#~ dirgenetrees = sys.argv[1]
#~ dirsubfamostortpickles = sys.argv[2]
nfreftree = sys.argv[1]
topoccurclade = sys.argv[2]
dirout = sys.argv[3]


dbcon, dbcur = rec_to_db.dbconnect(dbclue='phylariane', dbpwd='********')

reftree = tree2.ReferenceTree(fic=nfreftree)
reflabs = reftree.sort(reftree[topoccurclade].get_children_labels(), labels=True)
dtaxid_code = {}
for node in reftree:
	code = node.label()
	dtaxid_code[rec_to_db.get_taxid(code, dbcur)] = code
			
nsubfam = 0
# subfam-to-species matrix of replicon location
fmatout = open('%s/subfam_location_profile.mat'%dirout, 'w')
# species-to-location subfamily trees
ortlocout = '%s/mappping_subfam_locations'%dirout
if not os.access(ortlocout, os.F_OK):
	os.mkdir(ortlocout)
# table of translocation events
ftranslocout = open('%s/subfam_translocations.tab'%dirout, 'w')
dbfields = ["subfam_id", "name_txt", "hogenom_gene_id", "name", "locus_tag", "new_locus_tag", "genomic_beg", "genomic_end", "description"]
headerfields = ["subfamily", "organism", "hogenom_gene_id", "gene_name", "locus_tag", "new_locus_tag", "begin", "end", "description"]
translocfields = ["ancestor.code", "ancestor.location", "descendent.code", "descendent.location"]
ftranslocout.write('\t'.join(headerfields+translocfields)+'\n')
# sum over the genome of species-to-location trees
stateloc = copy.deepcopy(reftree)
transloc = copy.deepcopy(reftree)

	
tmptabsubfam = 'subfam_%s'%topoccurclade
tsubfams = rec_to_db.getSubfamFromPhyloPattern(dbcur, specificity=(tuple(reftree[topoccurclade].get_children_labels()),True), tempTable=tmptabsubfam)
# fetch annotations for all genes in db
dgene_annot = rec_to_db.get_all_gene_annotation(dbcur, cols=dbfields, join_clause=[('genome.gene2subfam', 'hogenom_gene_id'), (tmptabsubfam, 'subfam_id')])	
#~ testfam = '49RHIZOB_5155.1'
#~ tsubfams = (testfam,)
for subfam in tsubfams:
	#~ print subfam
	
	ortloc = copy.deepcopy(reftree)
	ortloc.resetPresence(state='absent')
	ortloc.cleanEvents()
	profile = rec_to_db.getSubfamOccurence(dbcur, subfam, occurence='count', returnDict=True)
	#~ print subfam, profile
	# filters subfamilies with multiple occurences (might be families without trees and possibly very large gene counts/genome)
	if max(profile.values())>1: continue
	ttloc = rec_to_db.getRepliconLoc(dbcur, subfam=subfam)
	for tloc in ttloc:
		hogenomid, taxid, fam, subfam, replicon = tloc
		code = dtaxid_code.setdefault(taxid, rec_to_db.get_code(taxid, dbcur))
		if replicon==None: replicon='?'
		ortloc[code].presenceAtNode(state=replicon)
	absents = []
	for node in profile:
		if profile[node]==0: absents.append(ortloc[node])
	ortloc.FitchPars(excludedNodes=absents)
	before = ortloc.newick(comment='presence', ignoreBS=True)
	try:
		ortloc.refineFitchPars(excludedNodes=absents, silent=True)
	except ValueError, e:
		print subfam
		print before
		raise ValueError, e
	
	#~ if subfam==testfam:
		#~ print ortloc[topoccurclade].newick(comment='presence', ignoreBS=True)
		#~ break
	ortloc.writePhyloProfileMatrix(fmatout, fam=subfam, reflabs=reflabs, leavesOnly=False, header=(nsubfam==0), counts=False)
	
	# find translocation events
	for node in ortloc[topoccurclade]:
		f = node.go_father()
		if f:
			fstates = set(f.state().split('|'))
			nstates = set(node.state().split('|'))
			if (fstates <= set(['?', '-']) or (nstates <= set(['?', '-']))): continue
			difstates = nstates - fstates
			if (difstates and not (difstates <= set(['?', '-']))):
				# replicon location of the subfamily changed between the father and the node, with a informative change
				floc = setToStr(fstates)
				nloc = setToStr(nstates)
				# sum the replicon translocation counts over the whole genome
				dtrans = transloc[node.label()].misc()
				dtrans[(floc, nloc)] = dtrans.get((floc, nloc), 0) + 1
				# get the genes of the subfamily under the father
				lspe = f.get_leaf_labels()
				supannots = [f.label(), floc, node.label(), nloc]
				for tloc in ttloc:
					hogenomid = tloc[0]
					if hogenomid.split('_')[0] in lspe:
						rec_to_db.write_gene_annot(ftranslocout, hogenomid, dgene_annot, fields=dbfields, supvalues=supannots)				
			
	# transform annotation
	ortloc.factorStateToDict(alternativeAsFraction=True)
	
	# save the location tree
	tree2.dump_pickle(ortloc, '%s/%s.ortloc.pickle'%(ortlocout,subfam))
	
	# sum the replicon location counts over the whole genome
	stateloc += ortloc
		
	nsubfam += 1
	
	sys.stdout.write('\r%d\t\t'%nsubfam)

fmatout.close()

stateloc += transloc
tree2.dump_pickle(stateloc, '%s/genome_synthesis.replicon_location.pickle'%(dirout))
stateloc.write_newick('%s/genome_synthesis.replicon_location.nwk'%(dirout), comment='locationcount')
stateloc.write_newick('%s/genome_synthesis.replicon_translocation.nwk'%(dirout), comment='translocationcount')
