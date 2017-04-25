import re
import os
import sys
import acnuc, getpass

if len(sys.argv) < 5:
	print "Usage: python get_unicopy_fams_query.py alndir list.codespe outfile minnumspe [excluded_species1, [excluded_species2, [...]]]"
	sys.exit(2)

alndir = sys.argv[1]
fspe = open(sys.argv[2], 'r')
fout = open(sys.argv[3], 'w')
minnumspe = int(sys.argv[4])

lexcluded = []
for arg in sys.argv[5:]:
	lexcluded.append(arg)

lnfaln = os.listdir(alndir)
lspe = []
for line in fspe:
	spe = line.rstrip('\n')
	if spe and (not spe in lexcluded):
		lspe.append(spe)
fspe.close()

lfam = []

loadlist = raw_input('\Load a restricted list of gene families\n\
From a query on agrogenom ACNUC db? (q)\n\
From a file? (f)\n\
No (anything else)\n ')
if loadlist == 'q':
	agropwd = getpass.getpass("agrogenom ACNUC db password? ")
	acnuc.open_socket(server_ip='pbil.univ-lyon1.fr', port=5558)
	acnuc.opendb_pw(db_name='agrogenom', psswrd=agropwd)
	query = raw_input('logic query string (ACNUC language)? ')
	dlist = acnuc.proc_requete(query, 'fams')
	elt = acnuc.nexteltinlist(dlist['lrank'], 0)
	while elt['next']>0:
		lfam.append(elt['name'])
		elt = acnuc.nexteltinlist(dlist['lrank'], elt['next'])
elif loadlist == 'f':
	nflfams = raw_input('file path? ')
	flfams = open(nflfams, 'r')
	for line in flfams:
		fam = line.rstrip(' \t\n')
		lfam.append(fam)


def nonzero_product(l):
	if l:
		p=1
		for x in l:
			if x != 0:
				p *= x
		return p
		
def which(l):
	w = []
	for i in range(len(l)):
		if l[i]: w.append(i)
	return w
		
for nfaln in lnfaln:
	faln = open("%s/%s"%(alndir,nfaln))
	if lfam:
		fam = nfaln.split('.')[0]
		if not fam in lfam:
			continue
	aln = faln.read()
	faln.close()
	pres = [0]*len(lspe)
	for n in range(len(lspe)):
		ref = re.findall(lspe[n], aln)
		pres[n] = len(ref)
	if nonzero_product(pres)==1 and len(which(pres))>=minnumspe:
		print nfaln
		fout.write("%s/%s\n"%(alndir,nfaln))
			
fout.close()

				
			
			
			
			
