#!/usr/bin/python

import sys, os, tree2


dirspegeneannot = sys.argv[1]
nfreftree = sys.argv[2]
nfout = sys.argv[3]
constrastingclades = ['N15', 'N10','N5','N4']
if len(sys.argv)>4 : maxrelax = sys.argv[4]
else: maxrelax = 2

lnfspegeneannot = os.listdir(dirspegeneannot)
fout = open(nfout, 'w')
reftree = tree2.ReferenceTree(fic=nfreftree)
dtag = {'related':'nb_related_strains_possessing_the_subfamily', 'remote':'nb_remote_strains_possessing_the_subfamily'}
nullstr = '\N'

for nfspegeneannot in lnfspegeneannot:
	contrast = None
	othertag = 'remote'
	nf = nfspegeneannot.split('/')[-1]
	cladespe = nf.split('.')[0]
	abspres = nf.split('.')[1]
	clade = reftree[cladespe]
	for contrastspe in constrastingclades:
		contrast = reftree[contrastspe]
		if clade.is_child(contrast):
			othertag = 'related'
			break
	else:
		contrast = reftree
		
	print clade.label(), 'vs. other in', contrast.label(), 'mith max', maxrelax, 'occurence in', othertag
	if abspres=='specific_absence': nbother = contrast.nb_leaves() - clade.nb_leaves()
	else: nbother = 0
	
	fspegeneannot = open("%s/%s"%(dirspegeneannot, nfspegeneannot), 'r')
	subfams = set()
	header = fspegeneannot.readline().rstrip('\n').split('\t')
	for line in fspegeneannot:
		dannot = dict(zip(header, line.rstrip('\n').split('\t')))
		subfam = dannot['subfamily']
		if subfam in subfams: continue
		else: subfams.add(subfam)
		other = dannot[dtag[othertag]]
		try:
			relax = abs(nbother - int(other))
			if relax > maxrelax: continue
		except ValueError:
			relax = nullstr
		fout.write('\t'.join([subfam, cladespe, abspres, str(relax)])+'\n')
	fspegeneannot.close()
	print cladespe, abspres
	
fout.close()
