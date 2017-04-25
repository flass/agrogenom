#!/usr/bin/python
# -*- coding: utf8 -*-

import sys
import os
import subprocess

if len(sys.argv) < 3:
	print "Usage: python pal2nal.py protalnlist path/nucfastadir path/outdir"
	sys.exit(2)

protalnlist = open(sys.argv[1], 'r')
nucfastadir = sys.argv[2]
outdir = sys.argv[3]
logfile = open("%s/pal2nal.err"%outdir, 'w')
nfnewfasta = "%s/new.fasta"%(outdir)
n = 0

for line in protalnlist:
	nfprotaln = line.rstrip('\n')
	fam = nfprotaln.split('/')[-1].rsplit('.', 1)[0]
	# checking protein alignment for 
	protaln = open(nfprotaln, 'r')
	lkeys = []
	for line in protaln:
		if line.startswith('>'):
			lkeys.append(line.split()[0]+'\n')
	protaln.close()
	# parse nucleic fasta file and builds dictionary of fasta headers to sequences
	nucfasta = open("%s/%s.fasta"%(nucfastadir, fam), 'r')
	d = {}
	for line in nucfasta:
		if line.startswith('>'):
			key = line.split()[0]+'\n'
		else:
			d[key] = d.setdefault(key, "") + line
	nucfasta.close()
	# create a temporary nucleic fasta file with sequences ordered as in protein alignment
	newfasta = open(nfnewfasta, 'w')
	for key in lkeys:
		newfasta.write(key+d[key])
	newfasta.close()
	# calls pal2nal
	fout = open("%s/%s.codon"%(outdir, fam), 'w')
	codon = subprocess.call(['/panhome/lassalle/pal2nal.v14/pal2nal.pl',nfprotaln, nfnewfasta, '-codontable', '11'], stdout=fout, stderr=logfile)
	fout.close()
	# updates counter
	n += 1
	sys.stdout.write("\r%d files treated"%n)
	sys.stdout.flush()
	
sys.stdout.write("\n")
logfile.close()
protalnlist.close()
os.remove(nfnewfasta)

