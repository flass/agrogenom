#!/usr/bin/python
# -*- coding: utf-8 -*-

import tree2
import shutil, os, sys
import shelve, cPickle
import time, copy, subprocess

cwd = os.getcwd()
famprefix = '"'
alnfileext = ".codon-gb.phyl"
phymltreefileext = "%s_phyml_tree.txt"%alnfileext
#~ RAxMLpath = "/home/lassalle/Programmes/RAxML-7.2.8-ALPHA/raxmlHPC-PTHREADS -T 2"
PhyMLpath = "/usr/bin/phyml"

def finddups(genetree, fam, maxlossrate, fmodif, force_phyml=False):
	
	nleaves = genetree.nb_leaves()
	nrecur = max(1000, nleaves*3)
	sys.setrecursionlimit(nrecur)
	print "%d leaves, set max recursion to %d"%(nleaves, nrecur)
	gc = copy.deepcopy(genetree)
	try:
		## annotates duplications (no reference to species tree) and provide list of conflicting nodes (duplications or other events) 
		duplabels, changedTopo = genetree.annotateDuplications(reftree=reftree, minbs=minbs, maxlossrate=maxlossrate, modifyTopology=True)
	except RuntimeError:
		# the model underlying moreParsimoniousScenario() algorithm expects a minimum of ancestral duplications, so the algorithm is trying to merge numerous recent dulpication events into a few ancestral ones.
		# there is a problem dealing with families too far from the model, such as that of EF-Tu translation factors show what seem to be strain-specific duplications in each strain (only duplication at leaves  = "cheries"), probably arising from continuous intr-genomic gene conversion between paralogs.
		# as this is the opposite of the model, the algorithm has great difficulty finding out a solution and may remain stuck in recursion. this has to be avoided, and the model being unadequat, the topology modification algorithm is not applied.
		genetree = gc
		duplabels, changedTopo = genetree.annotateDuplications(reftree=reftree, minbs=minbs, maxlossrate=maxlossrate, modifyTopology=False)
	sys.stdout.flush()
	#~ ## optimizes the branch lengths under the new topology using RAxML 
	#~ tmptreefile = "%s/%s.phb.tmp"%(workdir, fam)
	#~ genetree.write_newick(tmptreefile, comment=None)
	#~ alnfile = "%s/%s%s"%(alndir, fam, alnfileext)
	#~ outtreefile = "%s/%s.phb"%(modtreedir, fam)
	#~ command = RAxMLpath + " -t %s -s %s -n %s -f e -m GTRGAMMAI"%(tmptreefile, alnfile, outtreefile)
	#~ larg = command.split()
	#~ status = subprocess.check_call(larg)	# raises CallProcessError if RAxML does not return with status 0.
	#~ # RAxML call went well, erase the temporary input tree and load the output tree
	#~ os.remove(tmptreefile)
	
	outtreefile = "%s/%s%s"%(modtreedir, fam, '.phb')
	if changedTopo or force_phyml:
		# for sake of homogeneity of proceeding, better apply this step to all trees, even those without topolgy change
		# notably cures the possible consequences of rooting an unrooted tree with TPMS, yielding some branches without branch supports
		# !!! potential problem: these new branch supports may contradict the preceding reconciliation 
		# (clean procedure would need branch support re-estimation once after TPMS and a second time after duplication reconciliation)
		if changedTopo: fmodif.write("%s\n"%fam)
		## optimizes the branch lengths under the new topology using PhyML and updates the branch lengths/supports on the current GeneTree object
		alnfile = "%s/%s%s"%(alndir, fam, alnfileext)	# make sure to use the same alignment that for tree inference (i.e. that trimmed with Gblocks)
		genetree.recomputeBranchLengthsAndSupports(alnfile, outtreefile, d='nt', q='', m='GTR', c='8', a='e', v='e', o='lr', b='-3', quiet=True)
	else:
		genetree.write_newick(outtreefile, comment=None, mode='write')

	# write annotated tree in normed xml format
	genetree.write_phyloXML(nfout="%s/%s.xml"%(normtreedir, fam), treename=fam, normparams='default')
	## preprares unicopy subtree collection : defines subfamilies of orthologs and finds unicopy subtrees to test for transfer
	# list of subtrees of orthologs
	orthoSubtrees = genetree.get_all_ortholog_groups(duplabels)
	# dictionaries of sub-families
	orthoSubfams = {}
	# dictionaries of unicopy subtrees
	unicopySubtrees = {}
	for ost in orthoSubtrees:
		lost = ost.label()
		osf = ost.get_leaf_labels()
		orthoSubfams[lost] = osf
		dust = ost.combinationsUnicopySubtrees(duplabels, silent=True)
		unicopySubtrees[lost] = dust		
	# verifies that all leaves in the genetree are covered by a subfam
	coverage = reduce(lambda x,y: x+y, orthoSubfams.values())
	diff = set(coverage) ^ set(genetree.get_leaf_labels())
	if diff: raise IndexError, "leaves were not recovered in sub-families of %s:\n%s"%(fam, str(diff))
	# write records of genes in subfamilies
	for lost in orthoSubfams:
		fortho = open("%s/%s.%s"%(subfamdir, fam, lost), 'w')
		fortho.write("\n".join(orthoSubfams[lost]))
		fortho.close()
	# reduce the subtree set to a non-redundant set of topologically significant trees (at least 4 leaves)
	subtreeDB = {}
	for lost in unicopySubtrees:
		dust = unicopySubtrees[lost]
		for co in dust:
			# 'co' describes the combination of orthologous groups together with 'ost' (first member of dot-separated chain)
			ust = dust[co]
			kost = (lost, co)
			if ust.nb_leaves() > 3:
				lleaves = ust.get_leaf_labels()
				sleaves = set(lleaves)
				lleaves.sort()
				tleaves = tuple(lleaves)
				for tl in subtreeDB:
					if sleaves <= set(tl):
						# subtree leaf set is equivalent or included in a previous one, append it at the end of list
						subtreeDB[tl].append(kost)
						break
					elif sleaves > set(tl):
						# subtree leaf set includes a previous one, append it at the begin of list
						subtreeDB[tleaves] = [kost] + subtreeDB[tl]
						del subtreeDB[tl]
						break
				else:
					# subtree leaf set is neither included nor includes a previous one, create a new one
					subtreeDB[tleaves] = [kost]
	# write the dictionary of leaf sets to sub-families
	fl2fdict = open("%s/%s.dict"%(leaf2subfamdir, fam), 'w')
	fl2fdict.write(str(subtreeDB))
	fl2fdict.close()
	# write the (first of equivalent) unicopy subtree(s) in newick format
	nst = 0
	ft2ldict = open("%s/%s.dict"%(subtree2leafdir, fam), 'w')
	for tleaves in subtreeDB:
		nst += 1
		ft2ldict.write("%s.%d\t%s\n"%(fam, nst, ",".join(list(tleaves))))
		kost = subtreeDB[tleaves][0]
		ust = unicopySubtrees[kost[0]][kost[1]]
		# transform leaf labels into species names, erase any other label
		dust = ust.dictLeavesToSpecies()
		for node in ust:
			if node in dust: node.edit_label(dust[node])
			else: node.edit_label('')
		# write the tree with the root clear from any information (for unrooted rendering in newick output)
		ust.write_newick("%s/%s.%d.phb"%(subtreedir, fam, nst), mode='write', comment=None, unrooted=True) #, comment='unicity')
	ft2ldict.close()


#~ def main():

if len(sys.argv) < 7:
	print "Usage: python find_ancestral_duplications.py (u|c) path/input_file path/alignments path/outdir path/reftree min_branch_support min_duplication_consistency_score" 
	print "input_file can be a unique tree file (u) or a collection of gene trees as ouputed by TPMS (c)."
	print "\tif 'u', reconciled gene tree objects will be serialized as an individual cPickle file: famname.pickle;"
	print "\tif 'c', they will be stored in a common shelve database file 'allfam.shelve'. ALSO, THE OUPUT DIRECTORY 'path/outdir' IS ERASED AND RESET."
	sys.exit(2)

if sys.argv[1] == 'u':
	collection = False
elif sys.argv[1] == 'c':
	collection = True
else:
	collection = False
	print "treat input_file as a unique tree file"
	
nftree = sys.argv[2]
if collection:
	fam = None
else:
	fam = nftree.split('/')[-1].split('.')[0]
alndir = sys.argv[3]
outdir = sys.argv[4]
reftree = tree2.ReferenceTree(fic=sys.argv[5])
minbs = float(sys.argv[6])
maxlossrate = float(sys.argv[7])
	
if collection: ftpms = open(nftree, 'r')
modtreedir = outdir+"/modified_trees"
normtreedir = outdir+"/normed_trees"
subtreedir = outdir+"/subtrees"
subfamdir = outdir+"/subfams"
leaf2subfamdir = outdir+"/leaf2subfams"
subtree2leafdir = outdir+"/subtree2leaves"
shelves = outdir+"/treeshelves"
nfmodif = "%s/changed_topologies"%outdir
print nfmodif
if collection:
	print "erase/reset output directory"
	if os.access(outdir, os.F_OK):
		shutil.rmtree(outdir)
for dirout in [outdir, modtreedir, normtreedir, subtreedir, subfamdir, leaf2subfamdir, subtree2leafdir, shelves]:
	if not os.access(dirout, os.F_OK):
		os.mkdir(dirout)

workdir = os.getcwd()
print workdir


if collection:
	ntrees = 0
	# inititate database for storage of gene tree objects in pickled format
	treeshelf = shelve.open("%s/allfam.shelve"%shelves)
	ftpms = open(nftree, 'r')
	fmodif = open("%s/changed_topologies"%outdir, 'w')
	print 
	print "process gene trees"
	for line in ftpms:
		if line.startswith('%s'%famprefix):
			fam = line.strip('"\n').split('.')[0]
			ntrees += 1
		else:
			sys.stdout.write("\r%d\t\t%s\t\t"%(ntrees,fam))
			sys.stdout.flush()
			genetree = tree2.GeneTree(newick=line.rstrip('\n'), keep_comments=True, branch_lengths=True)
			finddups(genetree, fam, maxlossrate, fmodif, force_phyml=True)
			treeshelf[fam] = genetree
	treeshelf.close()
	ftpms.close()	
else:
	genetree = tree2.GeneTree(fic=nftree, keep_comments=True, branch_lengths=True)
	fmodif = open(nfmodif, 'a')
	finddups(genetree, fam, maxlossrate, fmodif, force_phyml=True)
	# stores the reconciled gene tree objects in pickled format (one file / object)
	fpickle = open("%s/%s.pickle"%(shelves, fam), 'w')
	cPickle.dump(genetree, fpickle, protocol=2)
	fpickle.close()
	
fmodif.close()
sys.stdout.write("\n")				

#~ if __name__ == "__main__":
    #~ main()
