#!/usr/bin/python

import sys

fdatin = open(sys.argv[1], 'r')
fdiccodeid = open(sys.argv[2], 'r')
fdatout = open(sys.argv[3], 'w')
dcodeid = {}
for line in fdiccodeid:
	lsp = line.rstrip('\n').split('\t')
	dcodeid[lsp[1]] = lsp[0] 

fdiccodeid.close()

code = None
sourcetag = False
taxontag = False
for line in fdatin:
	if line.startswith('ID'):
		code = line.split()[1].split('_')[0]
	if line.startswith('FT   source'):
		sourcetag = True
	if line.startswith('FT                   /db_xref="taxon:'):
		taxontag = True
	if line.startswith('FT   gene') and sourcetag:
		if not taxontag:
			l = 'FT                   /db_xref="taxon:%s"\n'%(str(dcodeid[code]))
			fdatout.write(l)
			print l, code
		else:
			taxontag = False
		sourcetag = False
		
	fdatout.write(line)
	
fdatin.close()
fdatout.close()
