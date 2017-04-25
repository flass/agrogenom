import tree2
import sys
import copy
if len(sys.argv) < 3:
	print 'Usage: make_reftree.py treefile dic_spe_code dic_code_taxid prefix'
	sys.exit(2)

nfintree = sys.argv[1]
nfdspecode = sys.argv[2]
nfdcodetaxid = sys.argv[3]
outprefix = sys.argv[4]

reftree = tree2.ReferenceTree(fic=nfintree)
# no internal labels and neither branch length nor branch support (compatible with Prunier 2.0 usage)
reftree.write_newick("%s.reftree.prunier"%outprefix, mode='write', branch_lengths=False, ignoreBS=True)

# naming internal nodes
reftree.complete_internal_labels()
reftree.complete_node_ids()

# standard reference tree
reftree.write_newick("%s.reftree.phb"%outprefix, mode='write', ignoreBS=True)

# write distance matrix of all nodes
reftree.write_matrix("%s.distance_nodes.mat"%outprefix, mode='write')
# write distance matrix of all leaves only
reftree.write_matrix("%s.distance_leaves.mat"%outprefix, mode='write', leavesOnly=True)
# write distance matrix of all lineages (from mid-branches)
reftree.write_matrix("%s.distance_lineages.mat"%outprefix, mode='write', distMethod='lineage_distance')
# write node ages tale
fnodeage = open("%s.age_nodes.tab"%outprefix, 'w')
flineage = open("%s.age_lineages.tab"%outprefix, 'w')
for node in reftree:
	lab = node.label()
	nid = node.nodeid()
	fnodeage.write("%s\t%d\t%f\n"%(lab, nid, node.mean_leaf_distance()))
	flineage.write("%s\t%d\t%f\n"%(lab, nid, node.lineage_age()))
fnodeage.close()
flineage.close()

# latin binomial instead of UniProt code
fdicspecode = open(nfdspecode,'r')
dcode_spe = {}
for line in fdicspecode:
	lsp = line.rstrip('\n').split('\t')
	dcode_spe[lsp[1]] = lsp[0]
fdicspecode.close()
latintree = copy.deepcopy(reftree)
for leaf in latintree.get_leaf_labels():
	latintree[leaf].edit_label(dcode_spe[leaf])
latintree.write_newick("%s.reftree.names"%outprefix,mode='write', ignoreBS=True)
latintree.write_newick("%s.reftree.names.bs"%outprefix,mode='write')

# quoted labels and ROOT mention (TPMS standard for reference tree)
tpmstree = copy.deepcopy(reftree)
for node in tpmstree:
	lab = node.label()
	node.edit_label('"%s"'%lab)
tpmstree.edit_label('ROOT')
tpmstree.write_newick("%s.reftree.tpms"%outprefix, mode='write', ignoreBS=True)#, branch_lengths=False)

# NCBI taxid
fdiccodetaxid = open(nfdcodetaxid,'r')
dcode_taxid = {}
for line in fdiccodetaxid:
	lsp = line.rstrip('\n').split('\t')
	dcode_taxid[lsp[1]] = int(lsp[0])
fdiccodetaxid.close()

# SQL dump for database 
analysisid = 1
i, dleftright = reftree.get_leftright_index()
lq = []

q = "INSERT INTO phylogeny.species_tree (species_tree_id, name, file_path, analysis_id, source_species_tree_id) VALUES (DEFAULT, '%s', '%s', %d, NULL);"%(outprefix, nfintree, analysisid)
lq.append(q)

q = "INSERT INTO phylogeny.species_node (sp_node_id, number, parent_sp_node_id, species_tree_id, tax_id, support, uniprot_id, left_num, right_num, branch_length) VALUES "
lq.append(q)

lv = []
for node in reftree:
	left, right = dleftright[node.nodeid()]
	if node.go_father(): parentid = node.go_father().nodeid()
	else: parentid = None
	node.set_taxid(dcode_taxid.get(node.label()))
	n = node
	while not node.taxid() and n.go_father():
		n = n.go_father()
		node.set_taxid(dcode_taxid.get(n.label()))
	if len(node.label())==5 and not node.label().startswith('ATU'):	uniprotid = node.label()
	else: uniprotid = None
	v = "(%d, '%s', %s, %d, %s, %s, '%s', %d, %d, %s)"%(node.nodeid(), node.label(), str(parentid), analysisid, str(node.taxid()), str(node.bs()), uniprotid, left, right, str(node.lg()))
	nv = v.replace("'None'", "NULL")
	nnv = nv.replace("None", "NULL")
	lv.append(nnv)
q = ('\t'+',\n\t'.join(lv)+';')
lq.append(q)

nfxmlout = "%s.reftree.xml"%outprefix
reftree.write_phyloXML(nfxmlout, treename=outprefix.rsplit('/', 1)[1], normparams=[])
fsqlout = open("%s.species_tree_node_tables.sql"%outprefix, 'w')
fsqlout.write('\n'.join(lq)+'\n')
fsqlout.close()
