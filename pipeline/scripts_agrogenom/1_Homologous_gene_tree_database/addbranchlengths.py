import tree2
import sys
if len(sys.argv) < 4:
	print "Usage: python addbranchlengths.py supports_tree lengths_tree output_tree"
	sys.exit(2)
	
nfsuptree =  sys.argv[1]
nflentree = sys.argv[2]
nfout = sys.argv[3]

suptree = tree2.Node(fic=nfsuptree, branch_lengths=False)
lentree = tree2.Node(fic=nflentree)
for node in lentree:
	lleaves = node.get_leaf_labels()
	bsn = suptree.map_to_node(lleaves)
	bsn.set_lg(node.lg())
suptree.set_lg(0)
suptree.write_newick(nfout, mode='write')
