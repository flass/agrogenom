#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys

nftpms = sys.argv[1]
outdir = sys.argv[2]

famprefix = '"'

nftree = None

ftpms = open(nftpms, 'r')
for line in ftpms:
	if line.startswith(famprefix):
		nftree = line.strip('%s\n'%famprefix)
	else:
		fout = open("%s/%s"%(outdir, nftree), 'w')
		fout.write(line)
		fout.close()
