#!/usr/bin/python
import os, sys
if len(sys.argv) > 1:
	target = sys.argv[1]
else:
	target = os.getcwd()
lfile = os.listdir(target)
lfile.sort()
for f in lfile:
	print f.split('.')[0]
