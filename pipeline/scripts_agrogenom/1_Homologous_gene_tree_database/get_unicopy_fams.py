import re
import os
import sys, getopt

def usage():
	s =  "Usage: python get_unicopy_fams.py alndir outfile [--universal=list.codespe] [--exclude='excluded_species1,[excluded_species2,[...]]']\n"
	return s


	
opts, args = getopt.getopt(sys.argv[1:], 'h:', ['universal=', 'exclude='])
dopt = dict(opts)

if ('-h' in dopt) or ('--help' in dopt):
	print usage()
	sys.exit(0)


if len(argv)<2:
	print "not enough arguments, must specify 'alndir' and 'outfile'"
	print usage()
	sys.exit(2)

alndir = sys.argv[1]
fout = open(sys.argv[2], 'w')

lspe = []
nflspe = dopt.get('--universal')
if nflspe:
	fspe = open(sys.argv[3], 'r')
	for line in fspe:
		spe = line.rstrip('\n')
		lspe.append(spe)
	fspe.close()

excluded = dopt.get('--exclude')
if excluded:
	for spe in excluded.split(','):
		lspe.remove(spe)

lnfaln = os.listdir(alndir)

def product(l):
	if l:
		p=1
		for x in l:
			p *= x
		return p
		
for nfaln in lnfaln:
	faln = open("%s/%s"%(alndir,nfaln))
	aln = faln.read()
	faln.close()
	pres = [0]*len(lspe)
	for n in range(len(lspe)):
		ref = re.findall(lspe[n], aln)
		pres[n] = len(ref)
	if product(pres) == 1:
		print nfaln
		fout.write("%s/%s\n"%(alndir,nfaln))
			
fout.close()

				
			
			
			
			
