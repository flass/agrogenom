#!/usr/bin/python

import os, sys
import rec_to_db, tree2


flnfgt = open(sys.argv[1], 'r')
lnfgt = []
for line in flnfgt:
	lnfgt.append(line.rstrip('\n'))
flnfgt.close()
reftree = tree2.ReferenceTree(fic=sys.argv[2])
dirprout = sys.argv[3]
dirouttrees = sys.argv[4]

for nfgt in lnfgt:
	print nfgt
	fam = nfgt.rstrip('\n').split('/')[-1].split('.')[0]
	#~ if not fam == "49RHIZOB_10225": continue
	genetree = tree2.GeneTree(fic=nfgt)
	# ignore families to small to have an meaningful tree
	if genetree.nb_leaves() < 4: continue
	# insure that each sub-root node has an annotated branch support	(creating a node at root make one of the new sons to not have supports)
	lbs = []
	for srn in genetree.get_children():
		if srn.is_leaf():
			# a branch separating one leaf from the rest of the tree have trivial branch support, skip
			break
		if srn.bs()!=None:
			lbs.append(srn.bs())
	else:
		for srn in genetree.get_children():
			if not srn.bs():
				try:
					srn.set_bs(lbs[0])
				except IndexError, e:
					print genetree
					raise IndexError, e
	# loads dictionary of species in the unicopy subtree to the genetree leaves
	dspe_leaves = genetree.dictSpeciesToLeafLabels()
	dleaf_spe = genetree.dictLeafLabelsToSpecies()
	nfprout = "%s/%s.1.prout"%(dirprout, fam)
	# ignore families where Prunier calculation has not succeeded
	if not os.access(nfprout, os.F_OK): continue
	fprout = open(nfprout, 'r')
	prout, ntrans, minbs, ny = rec_to_db.readPrunierOut(fprout, 0, True)
	fprout.close()
	if ntrans > 0:
		### identifies transfered leaves in gene tree
		ltrans, linfsup, ltransSpe = rec_to_db.parsePrunierOut(prout, ntrans, ny)
		#~ t = 0
		#~ while t < len(ltrans):
			#~ recspe = ltrans[t]
			#~ # mapping of Prunier transfer receptor and donor on gene tree
			#~ for spe in recspe:
				#~ leaflab = dspe_leaves[spe]
				#~ genetree.pop(leaflab)
		for llab in genetree.get_leaf_labels():
			genetree[llab].edit_label(dleaf_spe[llab])
		genetree.rootGeneTreeAsRefTree(reftree, ltransSpe, silent=True, findARoot=True)
		for llab in genetree.get_leaf_labels():
			genetree[llab].edit_label(dspe_leaves[llab][0])
		#~ print genetree
	genetree.write_newick("%s/%s.phb"%(dirouttrees, fam), mode='write')
