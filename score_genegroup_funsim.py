#!/usr/bin/env python

import os, sys
import psycopg2
import rec_to_db
import getopt, getpass
import tree2

from numpy import *
import random

from AIGO import logger

from AIGO.ReferenceSet import RefSet
from AIGO.FunctionalAnnotation import FuncAnnot
from AIGO.go.OBO import readGOoboXML


from AIGO.Similarity import GO_Similarity, GOSet_PWSimilarity


from AIGO.Analyse import AnalyseFA
from AIGO.Report import ReportFA

from AIGO.utils.Execute import batchExecute

#### default values for variables 

GOdir = "/pandata/lassalle/agrogenom/GeneOntology"
randsimdir = "%s/rand_gene_pair_sim"%GOdir
outdir = "%s/rand_gene_group_sim"%GOdir
nfGOtermgraph = "%s/go_daily-termdb.obo-xml"%GOdir

# agrogenomdb on phylarianneInstall
dbuser = 'agrogenomadmin'
dbhost = 'phylariane.univ-lyon1.fr'
dbname = 'agrogenomdb'
dbpwd = None

# example: windows of 7 genes in A. tumefaciens C58 genome
taxid = 176299
windowsizes = [7]
blockscores = False
reftree = None

metrics = ["funSimMax", "funSimAverage"]  
aspects = ['biological_process', 'molecular_function']

outcols = ['repli', 'windowsize', 'loc', 'leafblockid', 'nannot', 'labels', 'GOaspect', 'metric', 'meanPWSim', 'meanMaxSim']
analysisList=["coverage",  "richness", "numberAnnot", "redundancy", "specificity", "frequency", "informationContent"]
random.seed(1234567)

##### Functions

def repliconWindows(llab, windowsize, replitype):
	"""make gene windows around the circular or linear replicon"""
	windows = []
	if replitype=="circular":
		for i in range(len(llab)):
			j = i + windowsize
			if j < len(llab):
				windows.append(llab[i:j])
			else:
				windows.append(llab[i:len(llab)-1] + llab[:j-len(llab)])
	else:
		for i in range(len(llab)-windowsize):
			windows.append(llab[i:i+windowsize])
	return enumerate(windows)

def getBlockGenesToFile(nfbg, dbcon, dbcur, geneidcol='hogenom_gene_id', taxid=None, repli=None, eventtype=None, silent=True):
	q = "CREATE TEMPORARY TABLE blockgenes ON COMMIT DROP AS ("
	q += " SELECT chromosome, anc_block_id, leaf_block_id, %s, event_type, rec_sp_node_id, don_sp_node_id"%(geneidcol)
	q += " 	FROM genome.gene"
	q += " 	INNER JOIN blocks.gene_block_event_reflocs USING (gene_id)"
	if eventtype or repli or taxid:
		lw = []
		if taxid:
			lw.append("tax_id=%d"%(taxid))
		if repli:
			lw.append("chromosome='%s'"%(repli))
		if eventtype:
			lw.append("event_type='%s'"%(eventtype))
		q += " WHERE %s"%(' AND '.join(lw))
	q += " ORDER BY anc_block_id, leaf_block_id "
	q += ");"
	if not silent: print "retrieve block event genes from database\n", q
	dbcur.execute( q )
	fout = open(nfbg, 'w')
	dbcur.copy_to(fout, 'blockgenes', null='-1')
	fout.close()
	dbcon.commit()


def meanRandPWSim(G, FA, drepli_lab, randsimdir, GOaspects='all_aspects', compute=False, silent=True):
	"""evealuate the average functional similarity share by a random sample of pair of genes"""
	if GOaspects=='all_aspects':
		aspects = ['biological_process', 'molecular_function', 'cellular_component']
	else:
		aspects = list(GOaspects)
	dmeanrandsim = {}
	for repli in drepli_lab:
		llab = drepli_lab[repli]
		dmeanrandsim[repli] = {}
		#~ for GOaspect in FA.G.aspect:
		for GOaspect in aspects:
			dmeanrandsim[repli][GOaspect] = {}
			if not compute:
				for metric in metrics:
					nfrand = "%s/%s.%s.%s"%(randsimdir, repli, GOaspect, metric)
					if not os.access(nfrand, os.F_OK):
						compute = True
						break
			if compute:
				logger.info("Compute functional similarity on %s aspect between all random gene pair in %s"%(GOaspect, repli))
				dfout = {}
				for metric in metrics:
					nfrand = "%s/%s.%s.%s"%(randsimdir, repli, GOaspect, metric)
					dfout[metric] = open(nfrand, 'w')
				lsimmax = []
				lsimmean = []
				# list of genes covered by an annotation in this replicon
				lGP = list( set(FA.GPtoGO[GOaspect].keys()) & set(llab) )
				lGP.sort()
				for GP1 in lGP:
					lstrmax = []
					lsmax = []
					lstrmean = []
					lsmean = []
					for GP2 in lGP:
						if GP2 >= GP1:
							# only exlpore the lower triangular matrix
							continue
						else:
							GO1 = FA.GPtoGO[GOaspect][GP1]
							GO2 = FA.GPtoGO[GOaspect][GP2]
							maxsim, l = GOSet_PWSimilarity(G, GO1, GO2, FA=FA, metric="funSimMax")
							# profit from the fact that computing with "funSimMax" or "funSimAverage" metrics is the same (yield the same list l)
							lsmax.append(maxsim)
							lstrmax.append("%.3f"%maxsim)
							meansim = mean(l)
							lsmean.append(meansim)
							lstrmean.append("%.3f"%meansim)
					msmax = mean(lsmax)
					lsimmax.append(msmax)
					if lstrmax:
						dfout["funSimMax"].write(' '.join(lstrmax)+'\n')
						msmean = mean(lsmean)
						lsimmean.append(msmean)
						dfout["funSimAverage"].write(' '.join(lstrmean)+'\n')
						if not silent: print GP1, msmax, msmean
				for metric in dfout:
					dfout[metric].close()
			# read random similarity records
			for metric in metrics:
				lsim = []
				nfrand = "%s/%s.%s.%s"%(randsimdir, repli, GOaspect, metric)
				logger.info("Read in file %s for random gene pair similarities"%nfrand)
				foutrand = open(nfrand, 'r')
				for line in foutrand:
					ls = line.rstrip('\n').split()
					lfs = []
					for sim in ls:
						lfs.append(float(sim))
					lsim += lfs
				msim = mean(lsim)
				dmeanrandsim[repli][GOaspect][metric] = msim
				if not silent: print "on %s for aspect %s with %s metric: %f"%(repli, GOaspect, metric, msim)
	return dmeanrandsim

def genepair(GP1, GP2):
	gpair = [GP1, GP2]
	gpair.sort()
	return(tuple(gpair))
	
def groupPWSim(geneset, FA, G, GOaspect, default=None, silent=True):
	dgene_maxsim = dict(zip( metrics, [{}]*len(metrics))) 
	dpair_sim = dict(zip( metrics, [{}]*len(metrics)))
	nannot = 0
	geneset = list(geneset)
	geneset.sort()
	for GP1 in geneset:
		if not silent: print GP1
		lsim = dict(zip( metrics, [[]]*len(metrics)))
		if GP1 in FA.GPtoGO[GOaspect]:
			nannot += 1
			GO1 = FA.GPtoGO[GOaspect][GP1]
			# compute the matrix of gene pair similarities
			for GP2 in geneset:
				if GP2 <= GP1:
					# only compute the upper triangular matrix of similarities
					continue
				if not silent: print '\t', GP2,
				if GP2 in FA.GPtoGO[GOaspect]:
					GO2 = FA.GPtoGO[GOaspect][GP2]
					maxsim, l = GOSet_PWSimilarity(G, GO1, GO2, metric="funSimMax", FA=FA)
					meansim = mean(l)
					if not silent: print meansim, maxsim
				else:
					maxsim = default["funSimMax"]
					meansim = default["funSimAverage"]
					if not silent: print 'NA', 'NA'
				lsim["funSimMax"].append(maxsim)
				lsim["funSimAverage"].append(meansim)
				gpair = genepair(GP1, GP2)
				dpair_sim["funSimMax"][gpair] = maxsim
				dpair_sim["funSimAverage"][gpair] = meansim
			# find the maximum similarity of a gene with other genes of the group
			for metric in metrics:
				for GP2 in geneset:
				# fetch the previously computed similarities
					if GP2 < GP1:
						gpair = genepair(GP1, GP2)
						lsim[metric].append(dpair_sim[metric][gpair])
				dgene_maxsim[metric][GP1] = max(lsim[metric])
		else:
			# no similarity can be computed, take the default value
			for metric in metrics:
				dgene_maxsim[metric][GP1] = default[metric]
			for GP2 in geneset:
				if GP2 > GP1:
					# only fill the upper triangular matrix of similarities
					if not silent: print '\t', GP2, 'NA', 'NA'
					gpair = genepair(GP1, GP2)
					for metric in metrics:
						dpair_sim[metric][gpair] = default[metric]
	return dpair_sim, dgene_maxsim, nannot


def groupSummarySim(lnannot, window, FA, G, aspects, repli, fout, metrics=metrics, default=None, loc=0, leafblockid=0, taxid=None, silent=True):
	if taxid: geneset = ["%d.%s"%(taxid, u) for u in window]
	else : geneset = window
	labels = ','.join(window)
	windowsize = len(geneset)
	dmetric_sim = {}
	for GOaspect in aspects:
		if not default: d = None
		else: d = default[repli][GOaspect]
		if not silent: print GOaspect
		dpair_sim, dgene_maxsim, nannot = groupPWSim(geneset, FA, G, GOaspect, default=d, silent=silent)
		if not silent: print "default similarities (mean of random pairs ) replacing 'NA's:\n", d
		for metric in metrics:
			meanPWSim = mean(dpair_sim[metric].values())
			meanMaxSim = mean(dgene_maxsim[metric].values())
			dmetric_sim[metric] = [meanPWSim, meanMaxSim]
			strout = ('\t'.join([str(locals()[col]) for col in outcols])+'\n').replace('None', 'NA')
			fout.write(strout) 
	if not silent: print ' '.join([labels, ',', str(nannot), '/', str(windowsize), 'annotated'])
	lnannot.append(nannot)
	return dmetric_sim


def loadFA(G, norganism, dbcur, drepli, drepli_lab, taxid, aspects=aspects, metrics=metrics, analysisList=analysisList):
	inrefset = set([])
	for repli in drepli_lab:
		inrefset |= set(drepli_lab[repli])
	refSet = RefSet(organism=norganism, inSet=inrefset, refType="DB")
	FA = FuncAnnot(norganism, refSet, G, organism=norganism)
	FA.read_from_db(dbcur, replicons=drepli.keys())
	analyseFA = AnalyseFA()
	#print FA.GPtoGO['biological_process'].keys()
	analyseFA.largestSet([FA])
	logger.info("Largest sets of annotations:")
	logger.info("\t%d for %s" % (FA['largestSet']['All_aspects_of_GO'], FA.name))
	batchExecute(analysisList, analyseFA, [FA])
	#~ drepli_lab = {}
	#~ genelabeldir = "%s/genelabels/%s"%(outdir, norganism)
	#~ nflabels = "%s/%s_all_gene_labels"%(genelabeldir, norganism)
	#~ flab = open(nflabels, 'r')
	#~ for line in flab:
		#~ lsp = line.rstrip('\n').split('\t')
		#~ drepli_lab[lsp[0]] = drepli_lab.setdefault(lsp[0], []) + ["%s.%s"%(str(taxid), lsp[1])]	#[lsp[1]]
	#~ flab.close()
	return FA #, drepli_lab

def loadSubfamFA(G, norganism, dbcur, taxid, aspects=aspects, specificity=None, logicalop='AND', tempTable=None, analysisList=analysisList):
	inrefset = set(rec_to_db.getSubfamFromPhyloPattern(dbcur, specificity=specificity, logicalop=logicalop, tempTable=tempTable)
	for repli in drepli_lab:
		inrefset |= set(drepli_lab[repli])
	refSet = RefSet(organism=norganism, inSet=inrefset, refType="DB")
	FA = FuncAnnot(norganism, refSet, G, organism=norganism)
	FA.read_from_db(dbcur, getsubfams=True, specificity=specificity, logicalop=logicalop, fromTempTable=tempTable)
	analyseFA = AnalyseFA()
	analyseFA.largestSet([FA])
	logger.info("Largest sets of annotations:")
	logger.info("\t%d for %s" % (FA['largestSet']['All_aspects_of_GO'], FA.name))
	batchExecute(analysisList, analyseFA, [FA])
	return FA

def main(nfGOtermgraph, taxid, randsimdir, outdir, windowsizes, outcols, dbname, dbuser, dbhost, dbpwd, aspects=aspects, metrics=metrics, blockscores=blockscores, silent=True):
	
	dbcon, dbcur = rec_to_db.dbconnect(dbhost=dbhost, dbuser=dbuser, dbname=dbname, dbpwd=dbpwd)
	
	try:
		taxid = int(taxid)
		code = rec_to_db.get_code(taxid, dbcur)
	except ValueError:
		# try to fetch the (numerical) tax_id from the 5-letter code of the taxon
		code = taxid
		taxid = rec_to_db.get_taxid(code, dbcur)
	norganism = rec_to_db.get_organism_name(taxid, dbcur)
	drepli = rec_to_db.get_drepli_topology(taxid, dbcur)
	drepli_lab = rec_to_db.get_drepli_genes(taxid, dbcur, 'hogenom_gene_id')
	
	G = readGOoboXML(nfGOtermgraph)
	FA = loadFA(G, norganism, dbcur, drepli, drepli_lab, taxid, aspects=aspects, metrics=metrics)
	
	logger.info("Evaluate the average functional similarity share by a random sample of pair of genes") 
	dmeanrandsim = meanRandPWSim(G, FA, drepli_lab, randsimdir, GOaspects=aspects, silent=silent)
	
	logger.info("Evaluate the average functional similarity share by a random sample of groups of genes") 
	nfout = "%s/%s.genegroup_funsim"%(outdir, code)
	if not os.access(nfout, os.F_OK):
		fout = open(nfout, 'w')
		fout.write('\t'.join(outcols)+'\n')
	else:
		fout = open(nfout, 'a')

	default = dmeanrandsim
	for repli in drepli:
		logger.info("In replicon %s:"%repli)
		llab = drepli_lab[repli]
		for windowsize in windowsizes:
			logger.info("\tEvaluate the functional similarity shared by groups of gene (random gene sampling across the replicon)") 
			for i in range(len(llab)):
				lnannot = []
				window = random.sample(llab, windowsize)
				groupSummarySim(lnannot, window, FA, G, aspects, repli, fout, metrics=metrics, default=default, loc=-1, silent=silent)
			logger.info("\n\t%f / %d annotated in average for random groups"%(mean(lnannot), windowsize))
			logger.info("\tEvaluate the functional similarity clusters of gene share (systematic gene window scan)") 
			lnannot = []
			windows = repliconWindows(llab, windowsize, drepli[repli])
			for loc, window in windows:
				if not len(window)==windowsize:
					continue
				groupSummarySim(lnannot, window, FA, G, aspects, repli, fout, metrics=metrics, default=default, loc=loc, silent=silent)
			print '\n', mean(lnannot), '/', windowsize, 'annotated in average for contiguous groups of size %d in %s'%(windowsize, repli)
	fout.close()

	if blockscores:
		# load block event gene information and write it in a temporary file
		nfbg = "%s/%s.blockgenes"%(outdir, code)
		if not os.access(nfbg, os.F_OK):
			getBlockGenesToFile(nfbg, dbcon, dbcur, taxid=taxid, silent=silent)
		fbg = open(nfbg, 'r')
		
		lnannot = []
		geneset = []
		leafblockid = None
		bevt = None
		# initiate output files
		# will store informations combined by blocks in the 'blockinfos' output file
		nfbi = "%s/%s.blockinfos"%(outdir, code)
		fbi = open(nfbi, 'w')
		fbi.write('\t'.join(["chromosome", "anc_block_id", "leaf_block_id", "event_type", "rec_sp_node_id", "don_sp_node_id"])+'\n')
		# will store functional similarities among blocks in the 'eventblocks_funsim' output file
		nfbout = "%s/%s.eventblocks_funsim"%(outdir, code)
		fbout = open(nfbout, 'w')
		fbout.write('\t'.join(outcols)+'\n')
		# parse temporary file
		for line in fbg:
			[repli, abid, lbid, hoggid, evtt, rsni, dsni] = line.rstrip('\n').split('\t')
			if not leafblockid:
				leafblockid = lbid
				bevt = [repli, abid, lbid, evtt, rsni, dsni]
			if lbid == leafblockid:
				geneset.append(hoggid)
			else:
				if len(geneset) > 1:
					if not silent: print '\t', leafblockid, geneset
					# compute fun sim for the previous block and write it to output file
					dmetric_sim = groupSummarySim(lnannot, geneset, FA, G, aspects, repli, fbout, metrics=metrics, default=default, loc=-2, leafblockid=leafblockid)	#, taxid=int(taxid))
					# write block event info
					fbi.write('\t'.join(bevt)+'\n')
				geneset = [hoggid]
				leafblockid = lbid
				bevt = [repli, abid, lbid, evtt, rsni, dsni]
		fbg.close()
		fbi.close()
		fbout.close()
		# delete temporary file
		os.remove(nfbg)

def usage():
	s = "python score_gene_group_funsim.py [options]\n"
	s+= "-t int\n\ttaxon code or tax_id\n"
	s+= "-w range[,range[,...]]\n\tgene group window sizes to score, with 'range' as: int[:int]\n"
	s+= "-g path\n\tto GO term graph file\n"
	s+= "-r path\n\tto random gene pair similarity directory\n"
	s+= "-o path\n\tto random gene group similarity directory\n"
	s+= "-D str\n\tPostgreSQL database name\n"
	s+= "-H str\n\tPostgreSQL host server\n"
	s+= "-U str\n\tPostgreSQL user\n"
	s+= "-p str\n\tPostgreSQL user password\n"
	s+= "-b \n\tcompute block event gene group similarities\n"
	return s

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print usage()
		sys.exit(2)		
	options, args = getopt.getopt(sys.argv[1:], 't:g:r:o:D:H:U:p:w:b')
	if args:
		print "unused arguments:", args
	dopt = dict(options)
	if '-t' in dopt: taxid = dopt['-t']
	if '-g' in dopt: nfGOtermgraph = dopt['-g']
	if '-r' in dopt: randsimdir = dopt['-r']
	if '-o' in dopt: outdir = dopt['-o']
	if '-D' in dopt: dbname = dopt['-D']
	if '-U' in dopt: dbuser = dopt['-U']
	if '-H' in dopt: dbhost = dopt['-H']
	if '-p' in dopt: dbpwd = dopt['-p']
	if '-w' in dopt: 
		windowsizes = []
		winsz = dopt['-w'].split(',')
		for wins in winsz:
			if wins=='n': break
			try:
				windowsizes.append(int(wins))
			except ValueError:
				[a,b] = wins.split(':')
				windowsizes += range(int(a), int(b)+1)
	if '-b' in dopt: blockscores = True
				
	print "nfGOtermgraph:", nfGOtermgraph
	print "randsimdir:", randsimdir
	print "outdir:", outdir
	print "windowsizes:", windowsizes
	
	
	main(nfGOtermgraph, taxid, randsimdir, outdir, windowsizes, outcols, dbname, dbuser, dbhost, dbpwd, aspects=aspects, metrics=metrics, blockscores=blockscores, silent=True)	#reftree=reftree, 
