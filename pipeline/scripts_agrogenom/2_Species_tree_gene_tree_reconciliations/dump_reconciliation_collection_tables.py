import tree2, rec_to_db, os, sys, getopt

def usage():
	s =  "Usage: python dump_reconciliation_collection_tables.py "
	s += "Options (long | short equivalents):"
	s += "mandatory options \n"
	s += "{--reconciliation-collection= | -c} path/dir\n"
	s += "facultative options \n"
	s += "{--output-directory= | -o} path/dir, defaults to value of '--reconciliation-collection'\n"
	s += "\t[{--reconciliation-collection-id= | -i} 'ddmmyy[yy]']\n"
	s += "\t[{--reftree= | -r}path/file_newick]\n"
	s += "\t[{--output-dbtables | -t}]\n"
	s += "\t[{--output-xmls | -x}]\n"
	s += "\t[{--output-pickles | -p}]\n"
	s += "\t[{--verbose | -v}]\n"
	s += "\t[{--append | -a}]\n"
	return s

if len(sys.argv) < 2:
	print usage()
	sys.exit(2)

options, args = getopt.getopt(sys.argv[1:], 'c:i:r:o:atxpvh', \
["reconciliation-collection=", "reconciliation-collection-id=", "reftree=", "output-directory=", \
"output-dbtables", "output-xmls", "output-pickles", "verbose", "append"])
dopt = dict(options)

dirreccol = dopt.get('--reconciliation-collection', dopt['-c'])
dirout = dopt.get('--output-directory', dopt.get('-o', dirreccol))
reccolid = dopt.get('--reconciliation-collection-id', dopt.get('-i'))
if not reccolid:
	# assume the reconciliation_collection directory is in the form *_ddmmyy[yy]
	date = dirreccol.rstrip('/').rsplit('_', 1)[1]
	if isdigit(date):
		dd = date[0:2]
		mm = date[2:4]
		yy = date[4:]
		if len(yy)==2: YY = '20'+yy
		else: YY = yy
		reccolid = "%s-%s-%s"%(YY, mm, dd)
	else:
		reccolid = date
nfreftree = dopt.get('--reftree', dopt.get('-r'))
if nfreftree:
	reftree = tree2.ReferenceTree(fic=nfreftree)
	reftree.complete_node_ids()
outdbtab = ('--output-dbtables' in dopt) or ('-t' in dopt)
outxml = ('--output-xmls' in dopt) or ('-x' in dopt)
outpickles = ('--output-pickles' in dopt) or ('-p' in dopt)
integrate = (outdbtab or outxml or outpickles)
silent = (not (('--verbose' in dopt) or ('-v' in dopt)))
appendmode = ('--append' in dopt) or ('-a' in dopt)
if ('-h' in dopt):
	print usage()

dirgenetrees = "%s/reconciled_tree_pickles"%dirreccol
lnfgenetrees = os.listdir(dirgenetrees)
dbcon, dbcur = rec_to_db.dbconnect(dbclue='phylariane')
if integrate:
	tables = ['phylogeny.event', 'phylogeny.event_possible_species_node', "phylogeny.represented_event"]
	for table in tables:
		rec_to_db.createTempTable(dbcur, table, truncateOnCommit=True, silent=silent)
else:
	foutrepr = open("%s/represented_event_dump.tab"%dirout, 'w')
	
for nfgenetree in lnfgenetrees:
	fam = nfgenetree.rsplit('.', 1)[0]
	dnrnnodes, recgiid = rec_to_db.getNRNodeFamDict(fam, dbcur, outmode=3)
	genetree = tree2.load_pickle("%s/%s"%(dirgenetrees, nfgenetree))
	if integrate:
		rec_to_db.integrateReconciliations(fam, genetree, dirout, dbcur, reftree, reccolid=reccolid, complete=False, outxml=outxml, outpickles=outpickles, useTempTables=2, checkEvents=True, silent=silent)
	else:
		for node in genetree:
			nodeid = node.nodeid()
			nrnodeid = dnrnnodes[nodeid]
			eventid = node.eventid()
			if not eventid: raise ValueError, "no event_id for nr node %d (node %d in gene tree %s)"%(nrnodeid, nodeid, fam)
			foutrepr.write("%d\t%d\t%s\n"%(eventid, nrnodeid, reccolid))
	#~ sys.stdout.write('\r%s\t\t'%fam)
	

if outdbtab:	
	cwd = os.getcwd()
	nftmpdumpfile = "%s/tmp_dump.tab"%cwd
	for table in tables:
		temptable = table.split('.')[1]
		rec_to_db.dumpTableToFile(dirout, temptable, dbcur, delim='\t', append=appendmode, nftmpdumpfile=nftmpdumpfile, silent=silent)
else:			
	foutrepr.close()	
		
foutreccol = open("%s/reconciliation_collection_dump.tab"%dirout, 'w')
foutreccol.write("%s\t%s\n"%(reccolid, dirgenetrees))

dbcon.close()
