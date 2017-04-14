#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module serving as a specific API for AGROGENOMDB PostgreSQL database. 

It notably populates this database with descriptions of gene tree reconciliations 
based on Prunier software output for transfer detection.

"""

import sys
import os
import tree2
import numpy as np
import re
import copy
import cPickle
import psycopg2
import getpass, getopt
import time

check_bs = True

def dbg(fam):
	return 0
	#~ return (fam == "49RHIZOB_145_0")
	#~ return 1

### functions and table+columns definitions for filling AGROGENOMDB database

# dictionary of columns to fill in tables of AGROGENOMDB database
dtable_col = dict({ \
"phylogeny.gene_tree" : ['gene_tree_id', 'analysis_id', 'source_gene_tree_id', 'family_accession', 'family_type_id', 'file_path', 'type'], \
"phylogeny.reconciled_gene_tree" : ['rec_gi_id', 'filepath', 'tree_num', 'input_gene_tree_id', 'species_tree_id', 'analysis_id'], \
"phylogeny.reconciled_gene_tree_node" : ['nr_node_id', 'rec_gi_id', 'node_id', 'branch_support', 'branch_length'], \
"phylogeny.event" : ['event_id', 'rec_gi_id', 'event_type', 'anc_block_id', 'rec_gi_start_node', 'rec_gi_end_node'], \
"phylogeny.event_possible_species_node" : ['event_id', 'sp_node_id', 'characteristic', 'reference_node'], \
"phylogeny.represented_event" : ['event_id', 'nr_node_id', 'rec_col_id'], \
"phylogeny.reconciliation_collection" : ['rec_col_id', 'rec_tree_directory_path'], \
"analysis.prunier_inference" : ['prinf_id', 'unicopy_subtree_id', 'prunier_round', 'event_id', 'analysis_id', 'inverted_rec_don'], \
"analysis.prunier_inference_support" : ['prinf_id', 'branch_support'], \
"analysis.pruned_leaf" : ['prinf_id', 'gene_id'], \
"phylogeny.unicopy_subtree" : ['unicopy_subtree_id', 'subtree_name', 'subtree_file_path', 'prunier_outfile_path'], \
"phylogeny.subtree_implantation" : ['unicopy_subtree_id', 'prunier_congruent'], \
"phylogeny.node2subtree" : ['nr_node_id', 'unicopy_subtree_id'], \
"phylogeny.node2inference" : ['nr_node_id', 'prinf_id'], \
"phylogeny.gene2subtree" : ['gene_id', 'unicopy_subtree_id'], \
"blocks.ancestral_block" : ['anc_block_id', 'event_type'], \
"blocks.leaf_block" : ['leaf_block_id', 'anc_block_id'], \
"blocks.block_possible_species_node" : ['block_id', 'sp_node_id', 'characteristic', 'block_type', 'reference_node'], \
"blocks.gene2block" : ['leaf_block_id', 'gene_id'], \
"blocks.fam2block" : ['anc_block_id', 'family_accession'] \
})

devt_abrev = {'duplication':'D', 'transfer':'T', 'speciation':'S', 'gain':'G', 'loss':'L'}
dabrev_evt = {'D':'duplication', 'T':'transfer', 'S':'speciation', 'G':'gain', 'L':'loss'}

def getagrogenomdictdb(dbclue=None):
	if dbclue=='phylariane': 	
		dictdb = dict(dbhost='phylariane.univ-lyon1.fr', dbuser='agrogenomadmin', dbname='agrogenomdb', dbpwd=None)
	#~ elif dbclue in ['phylarianeInstall', '192.168.100.99']:
		#~ # agrogenomdb on phylarianneInstall
		#~ dictdb = dict(dbhost='192.168.100.99', dbuser='agrogenomadmin', dbname='agrogenomdb', dbpwd=None)
	elif dbclue in ['pbil-sgbd']:
		# agrogenomdb on pbil-sgbd
		dictdb = dict(dbhost='pbil-sgbd.univ-lyon1.fr', dbuser='lassalle', dbname='agrogenom', dbpwd=None)
	else:
		dictdb = None
	return dictdb
	
def dbconnect(**kw):
	"""return connection objects from connection parameters given separately or in a dictionary, or from a string indicating the DB to connect
	
	separated parameters override the ones in the dictionary that override the ones guessed from the string.
	"""
	dbhost = kw.get('dbhost')
	dbuser = kw.get('dbuser')
	dbname = kw.get('dbname')
	dbpwd = kw.get('dbpwd')
	dictdb = kw.get('dictdb')
	dbclue = kw.get('dbclue')
	if dictdb: ddb = dictdb
	else: ddb = getagrogenomdictdb(dbclue)
	dp = dict(zip(('dbhost', 'dbuser', 'dbname'),(dbhost, dbuser, dbname)))
	for p in dp:
		if not dp[p] and ddb:  
			dp[p] = ddb.get(p)
		if not dp[p]:
			raise KeyError, "Cannot guess database connection parameters...\n%s"%str(dp)
	if not dbpwd: dp['dbpwd']=getpass.getpass("password for %(dbuser)s on host %(dbhost)s for PostgreSQL database %(dbname)s? "%dp)
	else: dp['dbpwd'] = dbpwd
	dbcon = psycopg2.connect(dbname=dp['dbname'], host=dp['dbhost'], user=dp['dbuser'], password=dp['dbpwd'])
	dbcur = dbcon.cursor()
	return dbcon, dbcur

def createTempTable(dbcursor, table, truncateOnCommit=True, silent=True):
	"""Create session (transaction)-specific table"""
	tblsuffix = table.split('.')[1]
	s = "CREATE TEMPORARY TABLE %s ( LIKE %s )"%(tblsuffix, table)
	if truncateOnCommit: s += " ON COMMIT DELETE ROWS"
	s += ";"
	if not silent: print s
	dbcursor.execute( s )

### SQL db data import functions

def make_dvalues(lval, table, cols=dtable_col, justcheck=False):
	if isinstance(cols, dict):
		lcol = cols[table]
	elif isinstance(cols, list):
		lcol = cols
	if (len(lval) != len(lcol)):
		raise ValueError, "values list:\n%s\ndoes not match column list for table %s:\n%s"%(str(lval), table, str(lcol))
	if not justcheck:
		#~ d = {}
		#~ for i in range(len(lcol)):
			#~ d[lcol[i]] = lval[i]
		d = dict(zip(lcol, lval))
		return d
	
def parse_val_to_SQL(val, null='NULL', enumType=False):
	if val==None:
		sval = null
	elif isinstance(val, bool):
		if val==True:
			sval = 'TRUE'
		elif val==False:
			sval = 'FALSE'
	elif isinstance(val, str):
		if val in ['TRUE', 'FALSE', null, 'DEFAULT'] or enumType:
			sval = val
		else:
			sval = "'%s'"%val
	elif isinstance(val, long) or isinstance(val, float) or isinstance(val, int):
		sval = str(val)
	else:
		raise TypeError, "non handled type: %s"%str(type(val))
	return sval	
	
def SQLinsert(table, values={}, null='NULL', dtable_col=dtable_col, tmptable=False):
	"""make strings for PostgreSQL INSERT queries."""
	if values:
		if tmptable: s = "INSERT INTO %s "%table.split('.')[1]
		else: s = "INSERT INTO %s "%table
		
		if isinstance(values, str):
			# fills with a SELECT query
			if values.startswith('SELECT'):
				s += "%s;"%values
			else:
				raise ValueError, "'values' string is not a valid SELECT query: %"%value
		else:
			# fills with user values
			if isinstance(values, list):
				try:
					dval = make_dvalues(values, table, cols=dtable_col)
				except ValueError:
					lval = values
				else:
					return SQLinsert(table, dval, null=null, dtable_col=dtable_col, tmptable=tmptable)
			elif isinstance(values, dict):
				lcol = values.keys()
				lval = []
				for col in lcol:
					val = values[col]
					lval.append(val)
				s += "( %s ) "%(', '.join(lcol))
			else:
				raise TypeError, " %s is unvalid type for 'values', must be either a dictionary, list or a string (for a SELECT query)"%str(type(values))
			lsval = []
			for val in lval:
				sval = parse_val_to_SQL(val, null=null)
				lsval.append(sval)
#			print values, lsval
			s += "VALUES ( %s );"%(', '.join(lsval))
		return s	
		
def executeSQLinsert(dbcursor, table, values={}, null='NULL', dtable_col=dtable_col, tmptable=False, silent=True):
	"""dynamically execute insert queries on a PostgreSQL database"""
	s = SQLinsert(table, values=values, null=null, dtable_col=dtable_col, tmptable=tmptable)
	if not silent: print s
	dbcursor.execute(s)
	
def insertNewSerialVal(table, values, dbcursor, useTempTables=True, null='NULL', dtable_col=dtable_col, silent=True):
	"""Generic function to insert a new row in a table with auto-incrementing serial ids, with hard or temporary tables.
	
	'values' is a list of values following dtable_col order but ommiting the first collumn which corresponds to the auto-incremented serial.
	"""
	serialidseq = "%s_%s_seq"%(table, dtable_col[table][0])
	# fetch the serial id (and auto-increments it)
	dbcursor.execute( "SELECT nextval('%s');"%serialidseq)
	serialid = dbcursor.fetchone()[0]
	# creates an entry
	v = [serialid]+values
	executeSQLinsert(dbcursor, table, v, silent=silent, null=null, tmptable=useTempTables, dtable_col=dtable_col)
	return serialid
		
def writeSQLinsert(fdump, table, values={}, null='NULL'):
	"""function writing INSERT queries to a PostgreSQL dump file
	
	!!! makes huge files ; fillDataTableFile() must be prefered.
	"""
	s = SQLinsert(table, values=values, null=null)
	fdump.write(s+'\n')
		
def fillDataTableFile(ftable, table, values, force=False, delim='\t', null='\N'):
	"""function writing table to be read by PostgreSQL COPY queries"""
	if (not values) and (not force): 
		raise ValueError, "empty values are inserted"
	elif isinstance(values, list):
		# check user values match expected number of columns
		#~ make_dvalues(values, table, justcheck=True)
		#~ lsval = []
		#~ for val in values:
			#~ sval = parse_val_to_SQL(val, null=null)
			#~ lsval.append(sval)
		dval = make_dvalues(values, table)
	elif isinstance(values, dict):
		dval = values
		#~ # check user values contain expected columns and fix their order
		#~ lsval = []
		#~ for col in dtable_col[table]:
			#~ sval = parse_val_to_SQL(values[col], null=null, enumType=(col.endswith('type') or col=='characteristic'))
			#~ lsval.append(sval)
	else:
		raise TypeError, "expected a list or dict of values"
	lsval = []
	for col in dtable_col[table]:
		#~ sval = parse_val_to_SQL(values[col], null=null, enumType=(col.endswith('type') or col=='characteristic'))
		sval = parse_val_to_SQL(dval[col], null=null, enumType=(col.endswith('type') or col=='characteristic'))
		lsval.append(sval)
	s = "%s\n"%(delim.join(lsval))
	ftable.write(s)
		
	
def SQLupdate(table, where_clause, values={}, cols=dtable_col, null='NULL', tmptable=False):
	"""make strings for PostgreSQL UPDATE queries."""
	if values:
		if not tmptable: s = "UPDATE %s SET "%table
		else: s = "UPDATE %s SET "%table.split('.')[1]
		# fills with user values
		if isinstance(values, list):
			dval = make_dvalues(values, table, cols=cols)
			return SQLupdate(table, where_clause, values=dval)
		elif isinstance(values, dict):
			lsval = []
			for col in values:
				val = values[col]
				sval = parse_val_to_SQL(val, null=null)
				lsval.append("%s = %s"%(col, sval))	
			s += " %s %s;"%(', '.join(lsval), where_clause)
		else:
			raise TypeError, " %s is unvalid type for 'values', must be either a dictionary, list or a string (for a SELECT query)"%str(type(values))
		return s
	
def executeSQLupdate(dbcursor, table, where_clause, values={}, cols=dtable_col, null='NULL', tmptable=False, silent=True):
	"""dynamically execute update queries on a PostgreSQL database"""
	s = SQLupdate(table, where_clause, values=values, cols=cols, null=null, tmptable=tmptable)
	if not silent: print s
	dbcursor.execute(s)	
		
def writeSQLupdate(fdump, table, where_clause, values={}, cols=dtable_col, null='NULL'):
	"""function writing INSERT queries to a PostgreSQL dump file
	
	!!! makes huge files ; fillDataTableFile() must be prefered.
	"""
	s = SQLupdate(table, where_clause, values=values, cols=cols, null=null)
	fdump.write(s+'\n')
	
def SQLdelete(table, where_clause, tmptable=False):
	if not tmptable: s = "DELETE FROM %s "%table
	else: s = "DELETE FROM %s "%table.split('.')[1]
	s += "%s;"%(where_clause)
	return s
	
def executeSQLdelete(dbcursor, table, where_clause, tmptable=False, silent=True):
	s = SQLdelete(table, where_clause, tmptable=tmptable)
	if not silent: print s
	dbcursor.execute(s)	
	
### SQL db data export functions

def dumpTableToFile(dirsqlout, table, dbcursor, columns=None, delim='\t', null='\\N', append=True, nftmpdumpfile="tmp_dump.tab", silent=True):
	"""dump data content of a db table to a file using psql \copy command"""
	nfdumpfile = "%s/%s_dump.tab"%(dirsqlout, table)
	if append: nf = nftmpdumpfile
	else: nf = nfdumpfile
	f = open(nf, 'w')
	dbcursor.copy_to(f, table, sep=delim, null=null, columns=columns)
	f.close()
	if not silent:
		if columns: cols = " (%s)"%', '.join(list(columns))
		else: cols = ""
		print "\copy %s%s to '%s' with delimiter as '%s' null as '%s'"%(table, cols, nf, delim, null)
	if append:
		fdump = open(nfdumpfile, 'a')
		ftmpdump = open(nftmpdumpfile, 'r')
		tmpdump = ftmpdump.readlines()
		ftmpdump.close()
		fdump.writelines(tmpdump)
		if not silent: print "append to %s"%nfdumpfile
		fdump.close()
	
### AGROGENOM PostgreSQL db query functions
		
def colnames(cursor):
	"""get names of columns yielded by the last query executed with the cursor"""
	colnames = tuple()
	for col in cursor.description:
		colnames += (col[0],)
	return colnames

def characterize_event(eventid, cursor, useTempTables=False, silent=True):
	"""get the informations from an `event` entry that fully characterize it"""
	if not useTempTables:
		sevent = 'phylogeny.event'
		sepsn = 'phylogeny.event_possible_species_node'
	else:
		sevent = 'event'
		sepsn = 'event_possible_species_node'
	ssn = 'phylogeny.species_node'
	q =  "SELECT * "
	q += "FROM %s "%sevent
	q += "INNER JOIN %s USING (event_id) "%sepsn
	q += "INNER JOIN %s USING (sp_node_id) "%ssn
	q += "WHERE event_id=%d;\n"%eventid
	if not silent: print q
	cursor.execute( q )
	cucol = colnames(cursor)
	rows = cursor.fetchall()
	if not rows:
		if useTempTables > 1:
			# check also non-temporary table (in addition to temporary ones)
			return characterize_event(eventid, cursor, useTempTables=False, silent=silent)
		else:
			raise IndexError, "no entry for event_id %d in tables %s INNER JOIN %s"%(eventid, sevent, sepsn)
	refrecnodelab = None
	refdonnodelab = None
	lrec = []
	ldon = []
	d = None
	for row in rows:
		d = dict(zip(cucol, row))
		if d['characteristic'] in ['rec', 'location']:
			lrec.append(d['number'])
			if d['reference_node']:
				refrecnodelab = d['number']
		elif d['characteristic']=='don':
			ldon.append(d['number'])
			if d['reference_node']:
				refdonnodelab = d['number']
	startnode = d['rec_gi_start_node']
	endnode = d['rec_gi_end_node']
	evtype = d['event_type']
	recgiid = d['rec_gi_id']
	evchar  = (recgiid, evtype, startnode, endnode, refrecnodelab, refdonnodelab, set(lrec), set(ldon))
	return evchar
	
def probe_db_for_similar_events(evchar, cursor, useTempTables=False, silent=True):
	"""quickly get partial informations about `event` entries that resemble input event description"""
	if not useTempTables:
		sevent = 'phylogeny.event'
		sepsn = 'phylogeny.event_possible_species_node'
	else:
		sevent = 'event'
		sepsn = 'event_possible_species_node'
	valcols = zip(['rec_gi_id', 'event_type', 'rec_gi_start_node', 'rec_gi_end_node', 'recnode.number', 'donnode.number'], evchar[:6])
	whereclauses = []
	for col, char in valcols:
		sval = parse_val_to_SQL(char) 
		if sval=='NULL': whereclauses.append('%s IS %s'%(col, sval))
		else: whereclauses.append('%s = %s'%(col, sval))
	#~ dvalcols = dict(valcols)
	q =  "SELECT DISTINCT event_id "
	q += "FROM %s "%sevent
	q += "INNER JOIN ( "
	q += "	SELECT event_id, number FROM %s "%sepsn
	q += "	INNER JOIN phylogeny.species_node USING (sp_node_id) "
	q += "	WHERE characteristic IN ('rec', 'location') AND reference_node=TRUE "
	q += ") AS recnode USING (event_id) "
	q += "LEFT JOIN ( "
	q += "	SELECT event_id, number FROM %s "%sepsn
	q += "	INNER JOIN phylogeny.species_node USING (sp_node_id) "
	q += "	WHERE characteristic='don' AND reference_node=TRUE "
	q += ") AS donnode USING (event_id) "
	q += "WHERE %s;"%(' AND '.join(whereclauses))
	if not silent: print q
	cursor.execute( q )
	ltevents = cursor.fetchall()
	events = ()
	for tevent in ltevents:
		events += tevent
	if useTempTables > 1:
		# check also non-temporary table (in addition to temporary ones)
		events += probe_db_for_similar_events(evchar, cursor, useTempTables=False, silent=silent)
	return events
				
def getGeneInfo(apid, fields, dbcursor, returnDict=True):
	"""queries the gene informations from a database. Returns adictionary."""
#	fields = ['gene_id', 'replicon', 'genomic_begin', 'genomic_end', 'hogenom_gene_id', 'locus_tag', 'family_accession']
	query = "SELECT %s FROM genome.gene AS G "%(', '.join(fields))
	query += "WHERE G.hogenom_gene_id='%s'"%(apid)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	tgene = dbcursor.fetchone()
	if not returnDict:
		return tgene
	else:
		dgene = dict(zip(cols, tgene))
		return dgene
	
def getGeneTreeInfo(fam, fields, dbcursor, returnDict=True):
	"""queries the gene tree informations from a database. Returns a dictionary."""
	query = "SELECT %s FROM phylogeny.gene_tree INNER JOIN phylogeny.reconciled_gene_tree ON gene_tree_id=input_gene_tree_id "%(', '.join(fields))
	query += "WHERE family_accession='%s';"%(fam)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	tgenetree = dbcursor.fetchone()
	if not returnDict:
		return tgenetree
	else:
		dgenetree = dict(zip(cols, tgenetree))
		return dgenetree
	
def getFamilyInfo(fam, fields, dbcursor, returnDict=True):
	"""queries the gene tree informations from a database. Returns a dictionary."""
	query = "SELECT %s FROM genome.gene INNER JOIN phylogeny.gene_tree USING (family_accession) LEFT JOIN phylogeny.reconciled_gene_tree ON gene_tree_id=input_gene_tree_id "%(', '.join(fields))
	query += "WHERE family_accession='%s';"%(fam)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	if not returnDict:
		ttfamily = dbcursor.fetchall()
		return ttfamily
	else:
		tdfamily = ()
		for tfamily in dbcursor:
			tdfamily += (dict(zip(cols, tgenetree)),)
		return tdfamily
	
def getEventInfo(recgiid, startnodeid, fields, dbcursor, returnDict=True):
	"""queries the event informations from a database. Returns a tuple of dictionaries."""
	query = "SELECT %s FROM phylogeny.event "%(', '.join(fields))
	query += "WHERE rec_gi_id=%d AND rec_gi_start_node=%d;"%(recgiid, startnodeid)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	if not returnDict:
		ttevents = dbcursor.fetchall()
		return ttevents
	else:
		tdevent = ()
		for tevent in dbcursor:
			devent = dict(zip(cols, tevent))
			if 'event_type' in devent:	devent['eventtype'] = dabrev_evt[devent['event_type']]
			tdevent += (devent,)
		return tdevent
	
def getEventLocation(eventid, eventtype, dbcursor):
	"""perform query to get event location details in a format similar to the 'eventlocation' part of the return value of GeneTree.getdicevent()"""
	fields = ['number', 'sp_node_id', 'characteristic', 'reference_node']
	query = "SELECT %s "%', '.join(fields)
	query += "FROM phylogeny.species_node "
	query += "INNER JOIN phylogeny.event_possible_species_node USING (sp_node_id) "
	query += "WHERE event_id=%d ;"%eventid
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	ttlocs = dbcursor.fetchall()
	if eventtype=='transfer':
		lrec = []
		ldon = []
		reclab = None
		donlab = None
		# fetch transfer child id
		query = "SELECT rec_gi_end_node from phylogeny.event WHERE event_id=%d ;"%eventid
		dbcursor.execute( query )
		childid = dbcursor.fetchone()[0]
	else:
		lloc = []
		reflab = None
	for tloc in ttlocs:
		dloc = dict(zip(cols, tloc))
		c = dloc['characteristic']
		if c=='location':
			lloc.append(dloc['number'])
			if dloc['reference_node']==True:
				reflab = dloc['number']
		elif c=='rec':
			lrec.append(dloc['number'])
			if dloc['reference_node']==True:
				reclab = dloc['number']
		elif c=='don':
			ldon.append(dloc['number'])
			if dloc['reference_node']==True:
				donlab = dloc['number']
	if eventtype=='transfer':
		return (reclab, lrec, donlab, ldon, childid)
	else:
		return (reflab, lloc)
	
		
def getEventLocsInfo(eventid, fields, dbcursor, returnDict=True):
	"""queries the event informations from a database. Returns a tuple of dictionaries."""
	query = "SELECT %s FROM phylogeny.event_locs "%(', '.join(fields))
	query += "WHERE event_id=%d;"%(eventid)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	if not returnDict:
		ttevent = dbcursor.fetchall()
		return ttevent
	else:
		tdevent = ()
		for tevent in dbcursor:
			tdevent += (dict(zip(cols, tevent)),)
		return tdevent
		
def getTreeEventsInfo(fam, startnodeid, fields, dbcursor, returnDict=True):
	"""queries the event informations from a database. Returns a tuple of dictionaries."""
	query = "SELECT %s FROM phylogeny.tree_events "%(', '.join(fields))
	query += "WHERE family_accession='%s' AND rec_gi_start_node=%d;"%(fam, startnodeid)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	if not returnDict:
		ttevent = dbcursor.fetchall()
		return ttevent
	else:
		tdevent = ()
		for tevent in dbcursor:
			tdevent += (dict(zip(cols, tevent)),)
		return tdevent
		
def getTreeEventRefLocsInfo(fam, startnodeid, fields, dbcursor, returnDict=True):
	"""queries the event informations from a database. Returns a tuple of dictionaries."""
	query = "SELECT %s FROM phylogeny.tree_event_reflocs "%(', '.join(fields))
	query += "WHERE family_accession='%s' AND rec_gi_start_node=%d;"%(fam, startnodeid)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	if not returnDict:
		ttevent = dbcursor.fetchall()
		return ttevent
	else:
		tdevent = ()
		for tevent in dbcursor:
			tdevent += (dict(zip(cols, tevent)),)
		return tdevent
		
def getBlockEventRefLocsInfo(blockid, dbcursor, fields=['*'], returnDict=True):
	"""queries the event informations from a database. Returns a tuple of dictionaries."""
	query = "SELECT %s FROM block.block_event_reflocs "%(', '.join(fields))
	query += "WHERE block_id=%d;"%(blockid)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	if not returnDict:
		ttevent = dbcursor.fetchall()
		return ttevent
	else:
		tdevent = ()
		for tevent in dbcursor:
			tdevent += (dict(zip(cols, tevent)),)
		return tdevent
		
def getTreeEventBlockRefLocsInfo(eventid, dbcursor, fields=['*'], returnDict=True):
	"""queries the event informations from a database. Returns a tuple of dictionaries."""
	query = "SELECT %s FROM phylogeny.event "%(', '.join(fields))
	query += "INNER JOIN block.block_event_reflocs ON event.anc_block_id=block_event_reflocs.block_id "
	query += "WHERE event_id=%d;"%(eventid)
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	if not returnDict:
		ttevent = dbcursor.fetchall()
		return ttevent
	else:
		tdevent = ()
		for tevent in dbcursor:
			tdevent += (dict(zip(cols, tevent)),)
		return tdevent
		
def getNRNode(fam, nodeid, dbcursor):
	query = "SELECT nr_node_id "
	query += "FROM phylogeny.gene_tree "
	query += "INNER JOIN phylogeny.reconciled_gene_tree ON gene_tree_id=input_gene_tree_id "
	query += "INNER JOIN phylogeny.reconciled_gene_tree_node USING (rec_gi_id) "
	query += "WHERE family_accession='%s' AND node_id=%d; "%(fam, nodeid)
	dbcursor.execute( query )
	nrnodeid = dbcursor.fetchone()[0]
	return nrnodeid
		
def getNRNodeFamDict(fam, dbcursor, outmode=1):
	query = "SELECT rec_gi_id, node_id, nr_node_id "
	query += "FROM phylogeny.gene_tree "
	query += "INNER JOIN phylogeny.reconciled_gene_tree ON gene_tree_id=input_gene_tree_id "
	query += "INNER JOIN phylogeny.reconciled_gene_tree_node USING (rec_gi_id) "
	query += "WHERE family_accession='%s'; "%fam
	dbcursor.execute( query )
	ttnodes = dbcursor.fetchall()
	dnrnodes = {}
	recgiid = None
	for tnode in ttnodes:
		recgiid, nodeid, nrnodeid = tnode
		if outmode==1: dnrnodes[(recgiid, nodeid)] = nrnodeid
		elif outmode==2	: dnrnodes[nrnodeid] = (recgiid, nodeid)
		elif outmode==3	: dnrnodes[nodeid] = nrnodeid
	if outmode in [1,2]: return dnrnodes
	elif outmode==3: return dnrnodes, recgiid
		
def getNodeObservation(nrnodeid, dbcursor):
	query = "SELECT tested FROM phylogeny.node_observation "
	query += "WHERE nr_node_id=%d;"%(nrnodeid)
	dbcursor.execute( query )
	obs = dbcursor.fetchone()[0]
	return obs
	
def getNodeObservingReplicates(nrnodeid, dbcursor):
	query = "SELECT distinct(subtree_name) FROM phylogeny.node_observing_replicates "
	query += "WHERE nr_node_id=%d;"%(nrnodeid)
	dbcursor.execute( query )
	ttsup = dbcursor.fetchall()
	tsups = ()
	for tsup in ttsup:
		tsups += tsup
	return tsups
		
def getTransferObservation(eventid, dbcursor):
	query = "SELECT infered FROM phylogeny.transfer_observation "
	query += "WHERE event_id=%d;"%(eventid)
	dbcursor.execute( query )
	tobs = dbcursor.fetchone()
	if tobs: obs = tobs[0]
	else: obs = None
	return obs
	
def getTransferSupportingReplicates(eventid, dbcursor):
	query = "SELECT distinct(subtree_name) FROM phylogeny.transfer_supporting_replicates "
	query += "WHERE event_id=%d;"%(eventid)
	dbcursor.execute( query )
	ttsup = dbcursor.fetchall()
	tsups = ()
	for tsup in ttsup:
		tsups += tsup
	return tsups
		
def getEventLinkedFams(eventid, dbcursor):
	"""the tuple of families associated to this event in an ancestral block event"""
	query = "SELECT fam_accession FROM phylogeny.event "
	query += "INNER JOIN blocks.fam2block USING (anc_block_id) "
	query += "WHERE event_id=%d;"%(eventid)
	dbcursor.execute( query )
	ttfams = dbcursor.fetchall()
	tfams = ()
	for tfam in ttfams:
		tfams += tfam[0]
	return tfams
		
def getEventLinkedGeneCounts(eventid, dbcursor):
	"""the tuple of gene counts in leaf block events associated to this event"""
	query = "SELECT count(gene_id) FROM phylogeny.event "
	query += "INNER JOIN blocks.leaf_block USING (anc_block_id) "
	query += "INNER JOIN blocks.gene2block USING (leaf_block_id) "
	query += "WHERE event_id=%d "%(eventid)
	query += "GROUP BY leaf_block_id;"
	dbcursor.execute( query )
	ttcounts = dbcursor.fetchall()
	tcounts = ()
	for tcount in ttcounts:
		tcounts += tcount
	if not tcounts: tcounts = (1,)
	return tcounts
	
def get_code(taxid, dbcursor):
	dbcursor.execute("SELECT number FROM phylogeny.species_node WHERE tax_id=%d;"%taxid)
	code = dbcursor.fetchone()[0]
	return code
	
def get_taxid(code, dbcursor):
	dbcursor.execute("SELECT tax_id FROM phylogeny.species_node WHERE number='%s';"%code)
	taxid = dbcursor.fetchone()[0]
	return taxid
	
def get_organism_name(taxid, dbcursor):
	dbcursor.execute("SELECT name_txt FROM taxonomy.names WHERE tax_id=%d AND name_class='scientific name';"%taxid)
	nametax = dbcursor.fetchone()[0]
	return nametax
	
def get_drepli_genes(taxid, dbcursor, idcol='hogenom_gene_id', famWithTreeOnly=False):
	d = {}
	q = "SELECT chromosome, %s FROM genome.gene "%(idcol)
	if famWithTreeOnly: q += "INNER JOIN phylogeny.gene_tree USING (family_accession) "
	q += "WHERE tax_id=%d ORDER BY genomic_beg;"%(taxid)
	dbcursor.execute( q )
	ltcg = dbcursor.fetchall()
	for tcg in ltcg:
		d.setdefault(tcg[0], []).append(tcg[1])
	return d
	
def get_drepli_topology(taxid, dbcursor):
	dbcursor.execute("SELECT chromosome, topology FROM genome.chromosome where tax_id=%d;"%taxid)
	ltrt = dbcursor.fetchall()
	d = dict(ltrt)
	return d
	
def getRepliconLoc(dbcursor, hogenomid=None, subfam=None, fam=None):
	if hogenomid: whereclause = "hogenom_gene_id='%s'"%hogenomid
	elif subfam: whereclause = "subfam_id='%s'"%subfam
	elif fam: whereclause = "family_accession='%s'"%fam
	else: raise ValueError, "must specify at least one argument among 'hogenomid', 'subfam' or 'fam'"
	q =  "SELECT hogenom_gene_id, tax_id, family_accession, subfam_id, replicon "
	q += "FROM genome.gene_location "
	q += "WHERE %s ;"%whereclause
	dbcursor.execute(q)
	ttloc = dbcursor.fetchall()
	return ttloc
	
def get_all_gene_annotation(dbcursor, cols=None, join_clause=[], where_clause=None):
	if not cols: cols = ["name_txt", "family_accession", "hogenom_gene_id", "name", "locus_tag", "new_locus_tag", "chromosome", "genomic_beg", "genomic_end", "description"]
	q =  "SELECT %s "%(', '.join(cols))
	q += "FROM genome.gene "
	q += "LEFT JOIN genome.old2new_labels AS onl ON gene.locus_tag=onl.old_locus_tag "
	q += "INNER JOIN taxonomy.names USING (tax_id) "
	if join_clause:
		for table, col in join_clause:
			q += "INNER JOIN %s USING (%s) "%(table, col)
	q += "WHERE name_class='scientific name' "
	if where_clause: q += "AND " + where_clause.replace("WHERE ", "").replace("where ", "") 
	q += ";"
	print q+'\n'
	dbcursor.execute(q)
	#~ cols = colnames(dbcursor)
	ihg = cols.index("hogenom_gene_id")
	ttgenes = dbcursor.fetchall()
	dgene_annot = {}
	for annot in ttgenes:
		dgene_annot[annot[ihg]] = annot
	return dgene_annot

def write_gene_annot(fout, hogenomid, dgene_annot, fields, supvalues=[]):
	values = [str(val) for val in dgene_annot[hogenomid]]
	line = '\t'.join(values+supvalues)+'\n'
	fout.write(line)
	
def subfamspecificityclause(specificity=None, logicalop='AND'):
	"""generate join clause(s) for restricting a query of sufamily sets given clade specificity and/or phylogenetic pattern profile constrains
	this is to be specified by passing a [list of] tuples (clade, specificity) where 
	  clade is a species/ancestor code
	  specificity is  specifying
		either the presence or absence state coded by booleans: True or False, respectively
		or their exact count in the species/ancestor, coded by integers
		or more specific patterns, coded by strings: 'gain', 'loss', 'specific_presence', 'specific_absence'
	conditions defined by each tuple are linearly combined (whith the logical operator logicalop, which can be 'AND' or 'OR').
	can also include user-defined clauses passed as strings.
	"""
	dboolstr = {True:'t', False:'f'}
	q = ""
	#~ lw = []
	if isinstance(specificity, tuple): lspec = [specificity]
	elif isinstance(specificity, list): lspec = specificity
	elif isinstance(specificity, str): return specificity
	else: raise ValuError, "unproper specificity descriptor, needs a [list of] tuples (clade, specificity)"
	nspec = 0
	for tspec in lspec:
		if isinstance(tspec, str):
			q += tspec
		else:
			nspec += 1
			tag = 'spec%d'%nspec
			clade, spec = tspec
			if isinstance(clade, str): cladeclause = "number='%s'"%clade
			elif isinstance(clade, tuple): cladeclause = "number in %s"%str(clade)
			else: raise ValuError, "unproper clade descriptor, needs a single or tuple of clade code strings"
			
			if isinstance(spec, bool): 
				q += "inner join (select subfam_id from genome.phylogenetic_profile where %s AND present='%s') as %s using (subfam_id) "%(cladeclause, dboolstr[spec], tag)
			elif isinstance(spec, int):
				q += "inner join (select subfam_id from genome.phylogenetic_profile where %s AND count=%d) as %s using (subfam_id) "%(cladeclause, spec, tag)
			else:
				q += "inner join (select subfam_id from genome.specific_gene where %s AND specificity='%s') as %s using (subfam_id) "%(cladeclause, spec, tag)
	return q
	
def getSubfamFromPhyloPattern(dbcursor, specificity=None, getGenes=False, logicalop='AND', tempTable=None):
	"""
	fetch a set of subfamilies that can be defined by specific presence/absence profiles:
	this is to be specified by passing a [list of] tuples (clade, specificity) (see subfamspecificityclause())
	result can be stored in a temporary table of name tempTable.
	"""
	q = ""	
	if tempTable: q += "create temporary table %s as ( "%tempTable
	if not getGenes:
		q += "select distinct subfam_id from genome.gene2subfam "
	else:
		q += "select locus_tag, new_locus_tag, subfam_id "
		q += "from genome.gene2subfam "
		q += "inner join genome.gene using (hogenom_gene_id) "
		q += "left join genome.old2new_labels as onl on gene.locus_tag=onl.old_locus_tag "
	if specificity: q += subfamspecificityclause(specificity=specificity, logicalop=logicalop) 
	if tempTable: q += ") "
	q += ";"
	print q+'\n'
	dbcursor.execute(q)
	if tempTable:
		qt = "select * from %s ;"%tempTable
		dbcursor.execute(qt)
		print qt+'\n'
	ttsubfams = dbcursor.fetchall()
	if not getGenes:
		tsubfams = ()
		for tsubfam in ttsubfams:
			tsubfams += tsubfam
		return tsubfams
	else:
		return ttsubfams
	
def getSubfamOccurence(dbcursor, subfam, occurence='pres', returnDict=False):
	"""fetch the presence/absence profile of a subfamily.
	if abspres=True, return a list of tuples with (speciescode, abspresbool)
	else (default), return the list of species codes where the subfamily is present.
	"""
	if occurence=='pres':
		q = "SELECT number FROM genome.phylogenetic_profile WHERE subfam_id='%s' AND present IS TRUE ;"%subfam
	elif occurence=='abspres':
		q = "SELECT number, present FROM genome.phylogenetic_profile WHERE subfam_id='%s' ;"%subfam
	elif occurence=='count':
		q = "SELECT number, count FROM genome.phylogenetic_profile WHERE subfam_id='%s' ;"%subfam
	else:
		raise ValueError, 'unproper occurence definition: %s'%occurence
	dbcursor.execute(q)
	ttpres = dbcursor.fetchall()
	if occurence=='pres':
		lpres = []
		for tpres in ttpres:
			lpres.append(tpres[0])
		return lpres
	else:
		if returnDict:
			docc = {}
			for tpres in ttpres:
				docc[tpres[0]] = tpres[1]
			return docc
		else:
			return ttpres
	
def getSubfamSymbols(dbcursor): #, fromTempTable=None):	
	"""
	extract the symbol of a gene subfamily from the majoritary name among the subfamily genes, favoring the shortest symbol. 
	the symbol are uniquely attributed per family, if necessary by generating numerical suffixes.
	"""
	lastNone = 0
	def genUnicSymbol(sym, tempsymset, totsymset, lastNone=lastNone):
		if sym not in totsymset:
			return (sym, lastNone)
		for tsym in tempsymset:
			if tsym not in totsymset:
				return (tsym, lastNone)
		if sym=='None': n = lastNone
		else: n = 0
		usym = "%s%d"%(sym, n)
		while usym in totsymset:
			n += 1
			usym = "%s%d"%(sym, n)
		if sym=='None': lastNone = n
		return (usym, lastNone)
		
	dfamsym = {}
	q = "select subfam_id, symbol, count(hogenom_gene_id)"
	q += " from ( select hogenom_gene_id, upper(name) as symbol from genome.gene ) as genesym"
	q += " inner join genome.gene2subfam using (hogenom_gene_id)"
	#~ if fromTempTable: q += " inner join %s using (subfam_id)"%fromTempTable
	q += " group by subfam_id, symbol"
	q += " ;"
	print q+'\n'
	dbcursor.execute(q)
	tsym = dbcursor.fetchone()
	currsym = tsym
	totsymset = set([])
	tempsymset = set([])
	while tsym:
		if tsym[0]==currsym[0]:
			# same subfamily
			if tsym[1] and (tsym[1] not in totsymset):
				if (not currsym[1]) or (tsym[2] > currsym[2]) or ((tsym[2]==currsym[2]) and len(tsym[1])<len(currsym[1])):
					currsym = tsym
				else:
					tempsymset.add(tsym[1])
		else:
			# another subfamily
			# store previous subfamily symbol
			usym, lastNone = genUnicSymbol(currsym[1], tempsymset, totsymset, lastNone=lastNone)
			dfamsym[currsym[0]] = usym
			totsymset.add(currsym[1])
			# start exploring symbols within the next subfamily
			currsym = tsym
			tempsymset = set([])
			
		tsym = dbcursor.fetchone()
	
	return dfamsym
		
		
def getSubfamDescriptions(dbcursor): #, fromTempTable=None):	
	"""
	extract the functional description of a gene subfamily from the majoritary description among the subfamily genes, favoring the longest description.
	"""

	dfamdesc = {}
	q = "select subfam_id, descr, count(hogenom_gene_id)"
  	q += " from ( select hogenom_gene_id, lower(description) as descr from genome.gene ) as genedesc"
  	q += " inner join genome.gene2subfam using (hogenom_gene_id)"
	#~ if fromTempTable: q += " inner join %s using (subfam_id)"%fromTempTable
	q += " group by subfam_id, descr"
	q += " ;"
	print q+'\n'
	dbcursor.execute(q)
	tdesc = dbcursor.fetchone()
	currdesc = tdesc
	while tdesc:
		if tdesc[0]==currdesc[0]:
			# same subfamily
			if tdesc[1]:
				if (not currdesc[1]) or (tdesc[2] > currdesc[2]) or ((tdesc[2]==currdesc[2]) and len(tdesc[1])>len(currdesc[1])):
					currdesc = tdesc
		else:
			# another subfamily
			# store previous subfamily description
			dfamdesc[currdesc[0]] = currdesc[1]
			# start exploring descriptions within the next subfamily
			currdesc = tdesc
			
		tdesc = dbcursor.fetchone()
	
	return dfamdesc

def getSubfamGOAnnotationFile(dbcursor, dirout, prefix="all", specificity=None, logicalop='AND', fromTempTable=None, dfamdesc=None, dfamsym=None, taxid=0, annotator='FL'):
	"""
	extract data of GO annotation of subfamilies in a GAF text format.
	results from a previous query selecting a subset of subfamilies can be loaded from a temporary table of name fromTempTable. 
	
	!!! the (symbol, description) combination does not necessarily match an actual gene annotation from the subfamily gene set.
	"""
	def makeGAF2lines(famgos, dfamsym, dfamdesc):
		if famgos[0][0]==None:
			return ['!gaf-version: 2.0\n']
		else:
			subfam = famgos[0][0]
			desc = dfamdesc[subfam]
			sym = dfamsym[subfam]
			annots = []
			for tgo in famgos:
				annot = ['Agrogenom', subfam, sym, '', tgo[1], 'NO_REF', tgo[2], 'NO_EVID_REF', tgo[3], str(desc), '',	'protein', 'taxon:%d'%taxid, t,	annotator, '', '']
				annots.append('\t'.join(annot)+'\n')
			return annots
			
	# get date of the day
	t = time.strftime("%Y%m%d", time.localtime())
	# get subfamily descriptions and symbols
	if not dfamdesc: dfamdesc = getSubfamDescriptions(dbcursor) #, fromTempTable=fromTempTable)
	if not dfamsym: dfamsym = getSubfamSymbols(dbcursor) #, fromTempTable=fromTempTable)
	# get subfamily GO annotations
	q =  "select subfam_go.subfam_id, go.accession, go.evidence, go.aspect" 
	q += " from genome.subfam_go"
	q += " inner join genome.go using (go_id)"
	if fromTempTable: q += " inner join %s using (subfam_id)"%fromTempTable
	if specificity: q += subfamspecificityclause(specificity=specificity, logicalop=logicalop) 
	q += " ;"
	lines = []
	print q+'\n'
	dbcursor.execute(q)
	currfamgos = [(None,)]
	tgo = dbcursor.fetchone()
	while tgo: 
		if tgo[0]==currfamgos[0][0]:
			# same subfam annotation 
			currfamgos.append(tgo)
		else:
			# another subfam annotation 
			# complete and save the previous subfam annotation
			lines += makeGAF2lines(currfamgos, dfamsym, dfamdesc)
			# load the new subfam annotation
			currfamgos = [tgo]
		tgo = dbcursor.fetchone()
				
	ffaout = open("%s/%s_subfamFA.gaf"%(dirout, prefix), 'w')
	ffaout.writelines(lines)
	ffaout.close()	
			
def completeLeafLabels(subtree, fam, nsubtree, subtreeid, dirsubtree2leafdict, dbcursor, useTempTables=False, fillDB=True):
	"""updates the subtree with the dictionary of species in the unicopy subtree to the genetree leaves"""
	dspe_leaves = {}
	dleaf_spe = {}
	fsubtree2leafdict = open("%s/%s.dict"%(dirsubtree2leafdict, fam), 'r')
	for line in fsubtree2leafdict:
		if line.startswith(nsubtree+'\t'):
			lleaves = line.rstrip('\n').split('\t')[1].split(',')
			for leaf in lleaves:
				spe = leaf.split('_')[0]
				dspe_leaves[spe] = leaf
				dleaf_spe[leaf] = spe
				# edits unicopy subtree leaf labels to protein names
				subtree[spe].edit_label(leaf)
				# get gene id from db
				if fillDB: 
					dbcursor.execute(" SELECT genome.gene.gene_id FROM genome.gene \
										WHERE genome.gene.hogenom_gene_id=%s", (leaf,))
					geneid = dbcursor.fetchone()[0]
					# store in database
					executeSQLinsert(dbcursor, 'phylogeny.gene2subtree', [geneid, subtreeid], tmptable=useTempTables)
			break # the for line loop
	else:
		raise IndexError, "no entry was found for subtree %s in dictionary file %s/%s.dict"%(nsubtree, dirsubtree2leafdict, fam)	
	fsubtree2leafdict.close()
	return dspe_leaves, dleaf_spe
	
def get_scientific_name_dict(reftree, dbcur, childtags=True):
	dcode_names = {} 	
	for node in reftree.get_sorted_children(order=0):
		taxid = get_taxid(node.label(), dbcur)
		sciname = get_organism_name(taxid, dbcur)
		if sciname in dcode_names.values():
			childscinames = get_child_scientific_names(node, sciname, dbcur, childscinames=[], childtags=childtags)
			childscinames.sort()
			sciname += ' subgroup [%s]'%(', '.join(childscinames))
		dcode_names[node.label()] = sciname
	return dcode_names

def get_child_scientific_names(node, sciname, dbcur, childscinames=[], childtags=True):
	
	def get_lone_strain_species(species, strain, fatsciname):
		if not (species in fatsciname) and not (' sp. ' in species):
			# strains that are the only representative of a species
			for sep in ['genomosp.', 'genomovar']:
				if sep in species:
					tag = species[species.index(sep):]
					break
			else:
				tag = species
		else:
			tag = strain
		return tag
	
	if childscinames: csn = childscinames
	else: csn = []
	fatsciname = get_organism_name(get_taxid(node.label(), dbcur), dbcur)	
	# explore the scientific names below
	for child in node.get_children():
		chisciname = get_organism_name(get_taxid(child.label(), dbcur), dbcur)		
		if chisciname != sciname:
			#~ print '\t', child.label(), chisciname, '*'
			if not childtags:
				if not chisciname in csn: csn.append(chisciname)
			else:
				if ' str. ' in chisciname:
					species, strain = chisciname.split(' str. ')
					#~ tag = get_lone_strain_species(species, strain, fatsciname)
					tag = strain
				else:
					for sep in ['genomosp.', 'genomovar', 'CFBP', 'NCPPB', 'ICPPB', 'Kerr', 'CFN', 'CIAT', 'ATCC']:
						if sep in chisciname:
							tag = chisciname[chisciname.index(sep):]
							break
					else:
						speciesstrain = chisciname.rsplit(' ', 1)
						if len(speciesstrain)>1:
							species, strain = speciesstrain
							#~ tag = get_lone_strain_species(species, strain, fatsciname)
							tag = strain
						else:
							tag = chisciname
				
				if not tag in csn: csn.append(tag.strip(' '))
		else:
			#~ print '\t', child.label(), chisciname,
			csn = get_child_scientific_names(child, sciname, dbcur, childscinames=csn, childtags=childtags)
	return csn
	
def getAllTransferMatrix(dbcursor, reftree, eventtable="phylogeny.event", blocks=False, tofile=None, treeorder=4):
	"""generates donor-to-reeiver transfer matrix"""
	if blocks: fields = ['anc_block_id']
	else: fields = ['event_id']
	fields += ['number', 'characteristic']
	query = "SELECT %s "%', '.join(fields)
	query += "FROM %s "%eventtable
	# transfers observed in block or as single events
	if blocks:
		query += "INNER JOIN blocks.ancestral_block USING (anc_block_id) "
		query += "INNER JOIN blocks.block_possible_species_node ON anc_block_id=block_id "
	else:
		query += "INNER JOIN phylogeny.event_possible_species_node USING (event_id) "
	query += "INNER JOIN phylogeny.species_node USING (sp_node_id) "
	query += "WHERE %s.event_type='T' "%(eventtable.split(' as ')[-1].split('.')[-1])
	# order by clauses
	if blocks: query += "ORDER BY anc_block_id "
	else: query += "ORDER BY event_id "
	query += ";"
	print query
	dbcursor.execute( query )
	cols = colnames(dbcursor)
	dmatindex = {}
	i = 0
	orderedchildrenlabels = [node.label() for node in reftree.get_sorted_children(order=treeorder)]
	for lab in orderedchildrenlabels:
		dmatindex[lab] = i
		i += 1
		
	transfermat = [[0.0]*i for j in range(i)]
	nevent = 0
	curreventid = ""
	scurrrec = set([])
	scurrdon = set([])
	tevent = dbcursor.fetchone()
	while tevent:
		eventid, code, char = tevent
		if eventid!=curreventid:
			if curreventid != "":
				# record the previous event
				# observed frequency of event (1) is distributed among the possible rec/don couples
				try:
					f = 1.0 / (len(scurrrec) * len(scurrdon))
				except ZeroDivisionError, e:
					print curreventid, "scurrrec", scurrrec, "scurrdon", scurrdon
					#~ raise ZeroDivisionError, e
				else:
					for rec in scurrrec:
						reci = dmatindex[rec]
						for don in scurrdon:
							doni = dmatindex[don]
							transfermat[reci][doni] += f
			curreventid = eventid
			scurrrec = set([])
			scurrdon = set([])
			nevent += 1
		if char=='rec':
			scurrrec.add(code)
		elif char=='don':
			scurrdon.add(code)
		tevent = dbcursor.fetchone()
		
		
	print nevent, "transfer events"
	print "blocks =", blocks
	print "eventtable =", eventtable
	if tofile:
		fout = open(tofile, 'w')
		fout.write('\t'.join(['donor']+orderedchildrenlabels)+'\n')
		for i in range(len(orderedchildrenlabels)):
			fout.write('\t'.join([orderedchildrenlabels[i]]+[str(int(f)) for f in transfermat[i]])+'\n')
		fout.close()
	return transfermat
		
### Prunier output parsing and interpretation functions
	
def readPrunierOut(fprout, nroot=None, oneroot=False):
	"""reads in Prunier ouput file"""
	if isinstance(fprout, file): prout = fprout.readlines()
	elif isinstance(fprout, list): prout = fprout
	else: raise TypeError, "wrong object type '%s' provided for 'fprout', must be an (newly) opened file or a list (of lines)\n%s"%(type(fprout), repr(fprout))
	n = -1
	y = 0
	ntrans = 0
	minbs = None
	ny = None
	kroot = -1
	dntrans = {}
	dny = {}
	for line in prout:
		n += 1
		if line.startswith("\tsupport value threshold for conflict: "):
			minbs = float(line.rstrip('\n ').split(': ')[1])
		elif oneroot and line.startswith("Corresponding Root"):
			nroot = int(line.rstrip('\n').split(' : ')[1])
		elif line.startswith("The two trees are the same"):
			break
		### searches for relevant root entry	
		if line.split(" ")[0] == "=====":
			kroot += 1
			if prout[n+1].startswith("***"): y = 2
			else: y = 1
			ntrans = int(prout[n+y].rstrip("\n").split(" ")[-1])
			ny = (n, y)
			if nroot:
				if (line == "===== Root %d EVENTS =====\n"%(nroot)):
					break
			else:
				dntrans[kroot] = ntrans
				dny[kroot] = ny
	if nroot: return prout, ntrans, minbs, ny
	else: return prout, dntrans, minbs, dny, kroot
	
def parsePrunierOut(prout, ntrans, ny):
	"""identifies transfered leaves in gene tree"""
	n, y = ny
	ltrans = []
	linfsup = []
	ltransSpe = []
	for i in range(ntrans):
		recspe = (prout[n+y+1+i*2].rstrip("\n").split(" : ")[-1]).split()
		ltrans.append(recspe)
		ltransSpe += recspe
		infsupport = (prout[n+y+2+i*2].rstrip("\n").split(" : ")[-1]).split()
		linfsup.append(infsupport)
	return ltrans, linfsup, ltransSpe
	
def invertPrunierRecDon(recleaves, subtree, dleaf_spe, ltrans, t):
	"""Swap donor and receptor clades because Prunier's choice of pruning does not fit the rooting of the gene tree.
	
	'recleaves' role change from receptor clade to donor clade ; 
	the new donor clade is the remaining part of 'subtree' (the other side of the branch where Prunier has cut)
	identification of this new receptor clade need appropriate rooting 
	in order to have 'recleaves' monophyletic and under the root (so that the other node under the root is the new donor)
	"""

	print "Prunier's choice of pruning does not fit the rooting of the gene tree: must swap donor and receptor clades"
	# i.e. 'recleaves' role change from receptor clade to donor clade ; the new donor clade is the remaining part of 'subtree' (the other side of the branch where Prunier has cut)
	# identification of this new receptor clade need appropriate rooting in order to have 'recleaves' monophyletic and under the root (so that the other node under the root is the new donor)
	altdonleaves = recleaves
	for initspe in altdonleaves:
		initnode = subtree[initspe]
		upnode = initnode.go_father()	
		# starting from one pruned leaf, find the upper node that contain only pruned leaves (last value of 'initnode')
		while set(upnode.get_leaf_labels()) <= set(altdonleaves):
			initnode = upnode
			upnode = initnode.go_father()
		if upnode != subtree:
			# brother of 'initnode' will serve as outgroup
			outgroupnode = initnode.go_brother()
			altrootst = copy.deepcopy(subtree)
			altrootst.newOutgroup(altrootst[outgroupnode.label()], branch_lengths=False)
			altrootdoncla = altrootst.map_to_node(altdonleaves)
			if altrootst.is_monophyletic(altdonleaves) and (altrootdoncla in altrootst.get_children()):
				break # the for initspe loop
	else:
		raise IndexError, "could not find the branch where to root"
	# find the new receptor clade
	altrootreccla = altrootdoncla.go_brother()
	recleaves = altrootreccla.get_leaf_labels()
	recspe = []
	for rl in recleaves:
		recspe.append(dleaf_spe[rl])
	# check there is no overlap between these species and those removed during subsequent pruning steps	
	lcommonspe = []
	tlastcommon = None
	for it in range(t+1, len(ltrans)):
		lspe = ltrans[it]
		commonspe = list(set(lspe) & set(recspe))
		if commonspe:
			lcommonspe += commonspe
			tlastcommon = it
	return recleaves, recspe, lcommonspe, tlastcommon
	
### Reconciliation mapping functions
	
def mapPrunierToGeneTree(fam, nfgt, genetree, lnfsubtrees, dirsubtree2leafdict, dirprout, dbcur, reftree, oneroot, checkEvents=True, fillDB=True, useTempTables=True, onlyFullGeneTrees=False, silent=True):
	"""parses a list of Prunier output file and map the inferred events to a gene tree (main function of the module)
	
	mapping data are outputed to a GeneTree object and to a database via dynamic queries.
	"""
	
	
	tgtid = None
	trecgiid = None
	### creates the database entry for the gene tree
	# check for existence of an entry for this fam
	if not fillDB:
		dbcur.execute( "SELECT gene_tree_id FROM phylogeny.gene_tree \
							WHERE family_accession=%s", (fam,))
		tgtid = dbcur.fetchone()
	if tgtid:
		gtid = tgtid[0]
	else:
		gtid = insertNewSerialVal('phylogeny.gene_tree', [None, None, fam, 1, nfgt, None], dbcur, useTempTables=useTempTables, silent=silent)
	### creates the database entry for the reconciled gene tree
	if not fillDB:
		dbcur.execute( "SELECT rec_gi_id FROM phylogeny.reconciled_gene_tree \
							WHERE input_gene_tree_id=%s", (gtid,))
		trecgiid = dbcur.fetchone()
	if trecgiid:
		recgiid = trecgiid[0]
	else:
		recgiid = insertNewSerialVal('phylogeny.reconciled_gene_tree', [None, None, gtid, 1, 1], dbcur, useTempTables=useTempTables, silent=silent)
		### creates database entries for all nodes of the reconciled gene tree
		for node in genetree:
			insertNewSerialVal('phylogeny.reconciled_gene_tree_node', [recgiid, node.nodeid(), node.bs(), node.lg()], dbcur, useTempTables=useTempTables, silent=silent)
		
	for nfst in lnfsubtrees:
		### loads unicopy subtree that was used as input for Prunier
		subtree = tree2.GeneTree(fic=nfst)
		nsubtree = nfst.rstrip('\n').split('/')[-1].rsplit('.', 1)[0]
		### searches the Prunier output file
		nfprout = "%s/%s.prout"%(dirprout, nsubtree)
		if os.access(nfprout, os.F_OK):
			fprout = open(nfprout, "r")
			registerednfprout = nfprout
		else:
			print "No unicopy subtree found for family %s"%nsubtree
			registerednfprout = None
			fprout = None
		if fillDB:
			subtreeid = insertNewSerialVal('phylogeny.unicopy_subtree', [nsubtree, nfst, registerednfprout], dbcur, useTempTables=useTempTables, silent=silent)
		else: 
			dbcur.execute( "SELECT unicopy_subtree_id FROM phylogeny.unicopy_subtree WHERE subtree_name=%s", (nsubtree,))
			subtreeid = dbcur.fetchone()[0]
		
		print "# # # # # # # # # # # # # # # # "
		print "%d\t%s"%(subtreeid, nsubtree)

		outtrans = []
		allspe = set(subtree.get_leaf_labels())			
		
		# loads dictionary of species in the unicopy subtree to the genetree leaves
		dspe_leaves, dleaf_spe = completeLeafLabels(subtree, fam, nsubtree, subtreeid, dirsubtree2leafdict, dbcur, fillDB=fillDB, useTempTables=useTempTables)
		
		# unicopy subtree was unrooted (trifurcated): must root it as in gene tree
		schildren = subtree.get_children()
		lsc = len(schildren)
		maxd = 0
		outgroup = None
		for i in range(lsc):
			# find the pair of subroot nodes that coalesce in genetree the further from root = the closest from each other in gene tree ; the outgroup is the remaining node
			scoal = genetree.coalesce([schildren[i], schildren[(i+1)%lsc]])
			if scoal.depth() > maxd:
				maxd = scoal.depth()
				outgroup = schildren[(i+2)%lsc]
		subtree.defineNode([outgroup])
		subtree.complete_internal_labels()
		
		if fillDB: 
			## map all nodes from subtree to one or several nodes of the full gene tree
			## thus get the list of node/branches of the full gene tree observed in the unicopy subtree
			dnodecollapsed = genetree.map_collapsed_nodes(subtree)
			for node in dnodecollapsed:
				for observednodlab in dnodecollapsed[node]:
					observednode = genetree[observednodlab]
					observednodeid = observednode.nodeid()
					# increment nrnodeid if never met before
					if useTempTables: rgtn = 'reconciled_gene_tree_node'
					else: rgtn = 'phylogeny.reconciled_gene_tree_node'
					q = "SELECT nr_node_id FROM %s "%rgtn
					q += "WHERE rec_gi_id=%s "%recgiid
					q += "AND node_id=%s;"%observednodeid
					dbcur.execute( q )
					tnrnodeid = dbcur.fetchone()
					if tnrnodeid:
						nrnodeid = tnrnodeid[0]
					else:
						raise IndexError, "nr_node_id should already exist in database for rec_gi_id=%d AND node_id=%d"%(recgiid, observednodeid)
					executeSQLinsert(dbcur, 'phylogeny.node2subtree', [nrnodeid, subtreeid], tmptable=useTempTables, silent=silent)
		
		if fprout:
			#### reads in Prunier output file
			prout, ntrans, minbs, ny = readPrunierOut(fprout, nroot, oneroot)
			fprout.close()
						
			#### parses Prunier relevant output				
			if ntrans > 0:
				
				### identifies transfered leaves in gene tree
				ltrans, linfsup, ltransSpe = parsePrunierOut(prout, ntrans, ny)
				lpostponedinverted = []
				
				if not onlyFullGeneTrees:
					# check if original implantation of the subtree in the full gene tree is congruent with that implied by reconciliation by Prunier
					cpst = copy.deepcopy(subtree)
					try:
						cpst.rootGeneTreeAsRefTree(reftree, ltransSpe, silent=1)
					except IndexError:
						# topology between reference tree and subtree still differ after having pruned the transfers, because of non-supported topological conflict
						congrimplant = None
					else:
						congrimplant = True
						if not cpst.hasSameRoot(subtree):
							congrimplant = False
							print "subtree implantation not congruent with Prunier"
						del cpst
					if fillDB: executeSQLinsert(dbcur, 'phylogeny.subtree_implantation', [subtreeid, congrimplant], tmptable=useTempTables, silent=silent)
						
				### inference of donor clade in gene tree, as referenced in reference tree
				t = 0
				print " - - - - - - - -"
				while t < len(ltrans):
					print "prunier_round:", t
					recspe = ltrans[t]
					# mapping of Prunier transfer receptor and donor on gene tree
					recleaves = []
					for spe in recspe:
						leaflab = dspe_leaves[spe]
						recleaves.append(leaflab)
					#~ print "recleaves:", recleaves
					
					# check if recleaves is a monophyletic clade
					if not subtree.is_monophyletic(recleaves):
						print "Prunier's choice of pruning does not fit the rooting of the gene tree: must swap donor and receptor clades"	
						recleaves, recspe, lcommonspe, tlastcommon = invertPrunierRecDon(recleaves, subtree, dleaf_spe, ltrans, t)
						
						if lcommonspe:
							# reduce the set of species to that not pruned elsewhere
							remainingrs = list(set(recspe) - set(lcommonspe))
							# must postpone this pruning step after the last pruning involving common species
							ltrans[tlastcommon+1:tlastcommon+1] = [remainingrs]
							linfsup[tlastcommon+1:tlastcommon+1] = [linfsup[t]]
							lpostponedinverted.append(tlastcommon+1)
							# goes to the next inference
							t += 1
							print "new receptor clade %s has other infered transfer events; must postpone this inference to round %d"%(str(recspe), tlastcommon+1)
							print " - - - - - - - -"
							continue # the while t loop	
						else:			
							inverted_rec_don = True
					elif t in lpostponedinverted:
						inverted_rec_don = True						
					else:
						inverted_rec_don = False				
						
					# mapping of Prunier transfer receptor and donor on input subtree
					subreccla = subtree.map_to_node(recleaves)
					subdoncla = subreccla.go_brother()		 
					donleaves = subdoncla.get_leaf_labels()		
					
					if not congrimplant	and subreccla.go_father().is_root():
						# transfer inferred at the basis of the badly rooted subtree is false ; there is a speciation at this node
						# but it means that there is a remaining conflict not treated by Prunier: a transfer should be inferred on the path between the node and the implicit root of Prunier
						# transfer inference will neither be stored in database nor annotated on gene tree object (and its phyloXML representation)
						eventid = None
						genetrans = None
						
					else:			
					
						# could use here the tree2.GeneTree.addTransferEventFromSubtree() function
					
						# mapping of Prunier transfer receptor and donor on full gene tree			
						genereccla = genetree.map_to_node(recleaves)
						genedoncla = genetree.map_to_node(donleaves)			
						donspe = genedoncla.dictLeafLabelsToSpecies().values()
			#			print ['mapped clade in gene tree:']+genereccla.get_leaf_labels()
								
						## mapping of prunier transfer receptor and donor nodes on reference tree 
						# simple/reference mapping (most recent common ancestors of leaf sets) 
						# and cautious mapping (takes into account species representation in gene tree as well as topological uncertainties given branch support)
						refreccla, lrec = genereccla.possibleReceptors(reftree, recspe, subtree)
						refdoncla, ldon = genedoncla.possibleDonors(reftree, donspe, recspe, subtree, minbs=minbs)
												
						## filling transfer info in normed tree
						# annotates the node above the donor and receptor clade
						# finds transfer node in full gene tree which is the MRCA in gene tree of receptor and donor nodes (as found by Prunier in unicopy subtree)
						genetrans = genetree.coalesce([genereccla, genedoncla])
						transid = genetrans.nodeid()
						# finds transfer recipient branch which is the one under transfer node leading to receptor node
						for gtc in genetrans.get_children():
							if (genereccla == gtc) or genereccla.is_child(gtc):
								generecline = gtc
								childid = generecline.nodeid()
								break
						else:
							raise IndexError, "%s (gtc) is not a descendant of %s (genetrans)"%(gtc.label(), genetrans.label())
						# !!!!! 'fatnorm' qui est le noeud de transfer annote doit etre le coalescent de 
						# node annotation
						genetrans.set_transfer(reclab=refreccla.label(), lrec=lrec, donlab=refdoncla.label(), ldon=ldon, childid=childid, reftree=reftree)
						#~ dfam_nodeinfereces[fam] = dfam_nodeinfereces.setdefault(fam, []) + [dict(node=genetrans, recid=refreccla.nodeid(), lrec=lrec, donid=refdoncla.nodeid(), ldon=ldon, childid=childid)]
						print "transfer node:", genetrans.label()
						print "recleaves:", recleaves
						print "donleaves:", donleaves
						print "lrec:", lrec, "ldon:", ldon
						print "inverted_rec_don:", inverted_rec_don
						print ""
						
						## storing information relative to Prunier inference into database
						evchar = (recgiid, 'T', transid, childid, refreccla.label(), refdoncla.label(), set(lrec), set(ldon))
						if checkEvents:
							# check that no identical event is already stored in database
							if not silent: print "evchar:", evchar
							similarevents = probe_db_for_similar_events(evchar, dbcur, useTempTables=useTempTables, silent=silent)
							if not silent: print "similarevents:", similarevents
						else:
							similarevents = ()
						for simeventid in similarevents:
							sevchar = characterize_event(simeventid, dbcur, useTempTables=useTempTables, silent=silent)
							if not silent: print "\tsevchar:", sevchar
							if evchar == sevchar:
								# current event is excatly the same than one previously recorded in table `event` ;
								# will use this event_id for following insertions ; skip creation of `event` and `event_possible_species_node` entries
								eventid = simeventid 
								if not silent: print "\t\teventid = simeventid:", eventid
								break
						else:
							# current event is different from previous records ; storing event in the database
							eventid = insertNewSerialVal('phylogeny.event', [recgiid, 'T', None, transid, childid], dbcur, useTempTables=useTempTables)
						
							# storing possible rec/don nodes in the database
							for don in ldon:
								epsn = reftree[don]
								epsnid = epsn.nodeid()
								if epsn == refdoncla: refnodebool = 'TRUE'
								else: refnodebool = 'FALSE'
								executeSQLinsert(dbcur, 'phylogeny.event_possible_species_node', [eventid, epsnid, 'don', refnodebool], tmptable=useTempTables, silent=silent)
							for rec in lrec:
								epsn = reftree[rec]
								epsnid = epsn.nodeid()
								if epsn == refreccla: refnodebool = 'TRUE'
								else: refnodebool = 'FALSE'
								executeSQLinsert(dbcur, 'phylogeny.event_possible_species_node', [eventid, epsnid, 'rec', refnodebool], tmptable=useTempTables, silent=silent)
					
					## storing Prunier inference in the database
					prinfid = insertNewSerialVal('analysis.prunier_inference', [subtreeid, t, eventid, 1, inverted_rec_don], dbcur, useTempTables=useTempTables)
					# storing leaves pruned at this round in the database
					for leaflab in recleaves:
						dbcur.execute( "SELECT gene_id FROM genome.gene \
											WHERE hogenom_gene_id=%s", (leaflab,))
						leafgeneid = dbcur.fetchone()[0]
						executeSQLinsert(dbcur, 'analysis.pruned_leaf', [prinfid, leafgeneid], tmptable=useTempTables, silent=silent)
					# storing bs supporting this inference in the database
					for support in linfsup[t]:
						executeSQLinsert(dbcur, 'analysis.prunier_inference_support', [prinfid, support], tmptable=useTempTables, silent=silent)		
					
					if genetrans:
						# find the set of branches in gene tree that may match the transfer branch
						for node in genetrans.path_to(genereccla):
							postransid = node.nodeid()
							if useTempTables: rgtn = 'reconciled_gene_tree_node'
							else: rgtn = 'phylogeny.reconciled_gene_tree_node'
							q = "SELECT nr_node_id FROM %s "%rgtn
							q += "WHERE rec_gi_id=%s "%recgiid
							q += "AND node_id=%s;"%postransid
							dbcur.execute( q )
							nrnodeid = dbcur.fetchone()[0]
							# link them to this inference of Prunier in the database
							executeSQLinsert(dbcur, 'phylogeny.node2inference', [nrnodeid, prinfid], tmptable=useTempTables, silent=silent)

					## pruning of transfered species for further donor inferences
					subtree.pop( subreccla )
					t += 1		
					print " - - - - - - - -"
					
	return recgiid
	
		
def integrateReconciliations(fam, genetree, outfiles, dbcur, reftree, reccolid=None, recgiid=None, inferTPMStransfers=False, duplications=True, speciations=True, checkEvents=True, complete=True, outxml=True, outpickles=True, **kw):
	"""complete the reconciliation of a gene tre with transfer detected by mapping duplication and speciations at remaining nodes
	
	mapping data are outputed to a GeneTree object and to a database via dynamic queries.
	"""
	silent = kw.get('silent', True)
	useTempTables = kw.get('useTempTables', 2)
	minbs = kw.get('minbs', 0.9)
	if minbs: minbs = float(minbs)
	
	### complete reconciliation
	if complete: genetree.completeReconciliation(minbs, reftree, inferTPMStransfers=inferTPMStransfers, duplications=duplications, speciations=speciations, silent=silent)
	
	### output
	
	## prepare directories
	normtreedir = "%s/normed_phyloxml_trees"%outfiles
	if outxml:
		if not os.access(normtreedir, os.F_OK):
			os.mkdir(normtreedir)
	nfnormtree = "%s/%s.xml"%(normtreedir, fam)
	if outpickles:
		dirtreepickles = "%s/reconciled_tree_pickles"%outfiles
		if not os.access(dirtreepickles, os.F_OK):
			os.mkdir(dirtreepickles)
		nfpickle = "%s/%s.pickle"%(dirtreepickles, fam)

	## database output
	# get the non-redundant node ids
	dnrnnodes, rgi = getNRNodeFamDict(fam, dbcur, outmode=3)
	if not recgiid: recgiid = rgi
	#~ # update the database reconciled_gene_tree entry for the file path
	#~ executeSQLupdate(dbcur, 'phylogeny.reconciled_gene_tree', "WHERE rec_gi_id=%d"%recgiid, [nfnormtree], ['filepath'], tmptable=useTempTables, silent=silent)
	# create event entries for duplications
	for node in genetree:
		nodeid = node.nodeid()
		if node.eventtype() in ['speciation', 'duplication','gain']:
			refloclab, llocs = node.eventloc()
			refdonlab = childid = None
			ldon = []
		elif node.eventtype()=='transfer':
			refloclab, llocs, refdonlab, ldon, childid = node.eventloc()
		else:
			raise ValueError, 'unknown event type'
		etype = devt_abrev[node.eventtype()]
			
		if not node.eventid():
			if checkEvents:
				# check that no identical event is already stored in database
				evchar = (recgiid, etype, nodeid, childid, refloclab, refdonlab, set(llocs), set(ldon))
				if not silent: print "evchar", evchar
				similarevents = probe_db_for_similar_events(evchar, dbcur, useTempTables=useTempTables) #, silent=silent)
				if not silent: print "similarevents", similarevents
			else:
				similarevents = ()
			for simeventid in similarevents:
				sevchar = characterize_event(simeventid, dbcur, useTempTables=useTempTables)
				if evchar == sevchar:
					# current event is excatly the same than one previously recorded in table `event` ;
					# will use this event_id for following insertions ; skip creation of `event` and `event_possible_species_node` entries
					eventid = simeventid 
					break
			else:	
				eventid = insertNewSerialVal('phylogeny.event', [recgiid, etype, None, nodeid, childid], dbcur, useTempTables=useTempTables) #, silent=silent)
				# storing possible species tree location nodes in the database
				for i, l in enumerate((llocs, ldon)):
					for epsnlab in l:
						if i==0:
							if ldon: char = 'rec'
							else: char = 'location'
							refnodebool = (epsnlab == refloclab)
						else:
							char = 'don'
							refnodebool = (epsnlab == refdonlab)
						epsn = reftree[epsnlab]
						epsnid = epsn.nodeid()
						executeSQLinsert(dbcur, 'phylogeny.event_possible_species_node', [eventid, epsnid, char, refnodebool], tmptable=useTempTables) #, silent=silent)
			node.set_eventid(eventid)
		else:
			#~ # current event was already recorded in table `event` ; skip creation of `event` and `event_possible_species_node` entries
			eventid = node.eventid()
		if reccolid:
			# record reference event
			nrnodeid = dnrnnodes[nodeid]
			executeSQLinsert(dbcur, 'phylogeny.represented_event', [eventid, nrnodeid, reccolid], tmptable=useTempTables) #, silent=silent)
							
	## file output
	if outxml:
		# write output normed phyloXML tree file
		genetree.write_phyloXML(nfnormtree, treename=fam, normparams='default')
	if outpickles:
		# store (assign or update 'fam' entry) annotated genetrees into the dictionary of the new persistent storage shelf
		tree2.dump_pickle(genetree, nfpickle)	

	
	for node in genetree:
		if not node.eventtype(): print "!!! no event inferred at", node.label(), ':\n', node

def main(lnfgt, jobid, tasktable, dirsubtrees, dirsubtree2leafdict, dirprout, nfreftrees, nroot, oneroot, outfiles, reccolid, starttime, tl=0, agropwd=None, force_db_overwrite=False, test_session=False, **kw):
	"""make the overall reconciliation for a list of gene tree(s) and store the result into a database and as tree objects"""
	
	silent = kw.get('silent', True)
	useTempTables = kw.get('useTempTables', 2)
	
	def select_gt(lnfgt, jobid, tasktable, dbconnec, dbcursor, test_session=False):
		"""Search for next gene tree task to do. First looks at potential task list, if not querries the task survey table in database"""
		nfgt = None
		if lnfgt!=None:
			if lnfgt: nfgt = lnfgt.pop(0)
		else:
			dbcursor.execute("LOCK %s ;"%tasktable)
			dbcursor.execute("SELECT gene_tree_file_path FROM %s WHERE current_status=0 LIMIT 1;"%tasktable)
			tnfgt = dbcursor.fetchone()
			if tnfgt: nfgt = tnfgt[0]
		if nfgt:
			if test_session: curflag=3
			else: curflag=1
			dbcursor.execute("UPDATE %s SET current_status=%d, date=now(), job_id=%d where gene_tree_file_path='%s'"%(tasktable, curflag, jobid, nfgt))
			dbconnec.commit()
		return lnfgt, nfgt
	
	ftrees = open(nfreftrees, "r")
	reftrees = ftrees.readlines()
	ftrees.close()
	reftree = tree2.ReferenceTree(newick=reftrees[nroot]) #, branch_lengths=False)
	reftree.complete_internal_labels(prefix = '.N')
	reftree.complete_node_ids()
	
	def check_time(starttime, tl, prevelapsedtime):
		"""perform remaining computation time check"""
		currtime = time.time()
		elapsedtime = currtime - starttime
		print "elapsed time", elapsedtime
		if tl and (elapsedtime - prevelapsedtime) > (tl - elapsedtime):
			# time to parse one family is higher than remaining time, must stop
			stop = True
		else:
			stop = False
		prevelapsedtime = elapsedtime
		return stop, prevelapsedtime
		
	## set AGROGENOM database connection
	agrocon, agrocur = dbconnect(dbhost=dbhost, dbuser=dbuser, dbname=dbname, dbpwd=dbpwd)
	
	if useTempTables:
		# must create temporary tables like ones in the current db
		for table in dtable_col:
			schema = table.split('.')[0]
			if schema in ['phylogeny', 'analysis']:
				createTempTable(insertcur, table, truncateOnCommit=True, silent=silent)

	# initiate time count
	prevelapsedtime = 0
	
	lnfgt, nfgt = select_gt(lnfgt, jobid, tasktable, agrocon, insertcur, test_session=test_session)
	while nfgt:
		
		fam = nfgt.split('/')[-1].split('.')[0]
		
		# time check
		stop, prevelapsedtime = check_time(starttime, tl, prevelapsedtime)
		if stop: break
		
		lnfst = []
		for nfst in os.listdir(dirsubtrees):
			if nfst.split('.')[0]==fam:
				lnfst.append("%s/%s"%(dirsubtrees, nfst))
				
		### loads gene tree 
		genetree = tree2.GeneTree(fic=nfgt)
		### grooms gene tree 
		genetree.complete_internal_labels()
		genetree.complete_node_ids()
		
		print "\n%s:\nMapping transfer event inferences from unicopy subtrees to full gene tree"%fam	
		recgiid = mapPrunierToGeneTree(fam, nfgt, genetree, lnfst, dirsubtree2leafdict, dirprout, agrocur, reftree, oneroot, silent=silent)	
		
		print "\nCompletion of reconciliation of gene tree knowing transfer events"
		integrateReconciliations(fam, genetree, outfiles, agrocur, reftree, reccolid=reccolid, recgiid=recgiid, **kw)
		
		if useTempTables:
			print "Export of SQL database to dump files"
			dirsqlout = "%s/sql_dump.%d"%(outfiles, jobid)
			if not os.access(dirsqlout, os.F_OK):
				os.mkdir(dirsqlout)
			cwd = os.getcwd()
			nftmpdumpfile = "%s/tmp_dump.%d.tab"%(cwd, jobid)
			for table in dtable_col:
				schema, temptable = table.split('.')
				if schema in ['phylogeny', 'analysis']:
					dumpTableToFile(dirsqlout, temptable, insertcur, delim='\t', append=True, nftmpdumpfile=nftmpdumpfile, silent=silent)
					
		if test_session and not useTempTables:
			agrocon.rollback()
		else:
			insertcur.execute("UPDATE %s SET current_status=2, date=now() where gene_tree_file_path='%s'"%(tasktable, nfgt))
			agrocon.commit()
			
		lnfgt, nfgt = select_gt(lnfgt, jobid, tasktable, agrocon, insertcur, test_session=test_session)
	else:
		# all reconciliations have been done, record the collection
		executeSQLinsert(dbcur, 'phylogeny.reconciliation_collection', [reccolid, outfiles], tmptable=False, silent=silent)
		
	if test_session:
		# resets flags of tasks done during test
		insertcur.execute("UPDATE %s SET current_status=0, date=now(), job_id=0 where current_status=3"%tasktable)
		agrocon.commit()
		
	agrocon.close()
		

def usage():
		s =  "Usage: python rec_to_db.py [options]\n"
		s += "Script for populating AGROGENOM PostgreSQL database with results of reconciliation of gene tree histories\n\n"
		s += "Mandatory Options:\n"
		s += " Specify task schedule with:\n"
		s += " --task.survey.table=db_table_name\n\t\tname of the table in the database completion of gene tree tasks is documented;\nNew tasks are queried to this table unless '--list.full.gene.trees' is specified\n"
		s += " --dir.subtrees=path\n\t\tpath to directory containing unicopy subtrees\n"
		s += " --dir.subtree2leaf.dict=path\n\t\tpath to directory containing dictionary files of species <-> leaf labels in each subtree\n"
		s += " --dir.prunier.out=path\n\t\tpath to directory containing Prunier output files\n"
		s += " --reference.tree=path\n\t\tpath to reference\n"
		s += " --dir.output=path\n\t\tpath to output directory where phyloXML trees, pickles of gene tree objects and SQL data dump will be writen\n"
		s += " --jobid=int\n\t\tunique indentifier for parrallel jobs\n"
		s += " --queue=int\n\t\tnumber of allowed calculation hours\n"
		s += "Facultative Options:\n"
		s += " --reconciliation-collection-id=date\n\t\treference in the database to the coherent set of event recorded in output trees\n"
		s += " --list.full.gene.trees=path\n\t\tpath to file listing the paths of all full gene trees to be treated (for one sequential job)\n"
		s += " -n int\t\tthe number of the line in reference tree file corresponding to root desired for analysis\n\t\t(for multiroot Prunier transfer inference only);if root was fixed, provide any non-integer\n\n"
		s += " -b float\tminimum bootstrap threshold for event consideration\n\n"
		s += " -w database_password\n\t\tif not provided, the password is asked interactively\n\n"
		s += " -i\t\tpopulate database with INSERT commands directly on true tables;\n\t\tdefault behaviour uses INSERT and COPY commands on temporary tables to export dumps\n\n"
		s += " -t, --test\tset the script in test mode: override database checks for redundency in family entries\n\t\tand no permanent change to the database can be done (when '-i' option is set)\n"
		s += " -v, --verbose\tverbose mode\n"
		s += " -h, --help\tprint this help message and exit\n"
		return s 
	
if __name__=="__main__":
	
	try:
		options, args = getopt.getopt(sys.argv[1:], 'n:b:w:itvh', ["list.full.gene.trees=", "dir.subtrees=", "dir.subtree2leaf.dict=", "dir.prunier.out=", "reference.tree=", "dir.output=", "task.survey.table=", "reconciliation-collection-id=", "queue=", "jobid=", "test", "verbose", "help"])
	except getopt.GetoptError, err:
		print err
		print usage()
		sys.exit(2)	
	dopt = dict(options)
	#~ print dopt
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(1)
	
	
	# mandatory options	
	tasktable = dopt["--task.survey.table"]
	dirsubtrees = dopt["--dir.subtrees"]
	dirsubtree2leafdict = dopt["--dir.subtree2leaf.dict"]
	dirprout = dopt["--dir.prunier.out"]
	nfreftrees = dopt["--reference.tree"]
	outfiles = dopt["--dir.output"]
	jobid = int(dopt["--jobid"])
	
	# facultative options
	reccolid = dopt.get("--reconciliation-collection-id", time.strftime("%Y-%m-%d"))
	
	if "--list.full.gene.trees" in dopt:
		flnfgt = open(dopt["--list.full.gene.trees"], 'r')
		lnfgt = []
		for line in flnfgt:
			lnfgt.append(line.rstrip('\n'))
		flnfgt.close()
	else:
		lnfgt = None
	
	queue = dopt.get("--queue")
	starttime = time.time()
	print "Start job %d at %s"%(jobid, time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
	if queue:
		tl = float(queue) * 3600
		print "time limit = %fs (%sh)"%(tl, queue)
	else:
		tl = 0


	if '-n' in dopt:
		try:
			nroot = int(dopt['-n'])
			oneroot = False
		except ValueError :
			nroot = 0
			oneroot = True
	else:
		nroot = 0
		oneroot = True
		
	minbs = float(dopt.get('-b', 0.9))
	agropwd = dopt.get('-w')
	useTempTables = (not ('-i' in dopt))
	test_session = (('-t' in dopt) or ('--test' in dopt))
	silent = (not (('-v' in dopt) or ('--verbose' in dopt)))
	
	main(lnfgt, jobid, tasktable, dirsubtrees, dirsubtree2leafdict, dirprout, nfreftrees, nroot, oneroot, outfiles, reccolid, starttime, tl, minbs=minbs, agropwd=agropwd, useTempTables=useTempTables, test_session=test_session, silent=silent)
