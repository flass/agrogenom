#!/usr/bin/python
"""to build genome.phylogenetic_profile table in agrogenomdb"""
import sys, os, tree2

nfphyloprofile = sys.argv[1]
nfout = sys.argv[2]

fphyloprofile = open(nfphyloprofile, 'r')
fout = open(nfout, 'w')

dboolstr = {True:'T', False:'F', '1':'T', '0':'F'}

header = fphyloprofile.readline().rstrip('\n').split('\t')[1:]
for line in fphyloprofile:
	lsp = line.rstrip('\n').split('\t')
	subfam = lsp[0]
	profile = lsp[1:]
	for i in range(len(profile)):
		try:
			fout.write('\t'.join([subfam, header[i], dboolstr[bool(int(profile[i]))], profile[i]])+'\n')
			#~ fout.write('\t'.join([subfam, header[i], profile[i]])+'\n')
		except KeyError, e:
			print [subfam, header[i], profile[i]]
			raise KeyError, e

fphyloprofile.close()
fout.close()
