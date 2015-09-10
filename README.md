# agrogenom
Python modules and scripts for the construction of Agrogenom database

'rec_to_db.py' is an API for the loading and construction of the phylogenomic database Agrogenom from a collection of reconciled gene trees.
[depend on tree2 Python package https://github.com/flass/tree2]

'blockevents.py' (library) and 'getblockevents.py' (execution script) are used to build block evolutionary events given a collection of reconciled gene trees and a phylogenomic database documenting gene neighbourhoods.
[depend on tree2 Python package https://github.com/flass/tree2]

'score_genegroup_funsim.py' is used to compute Functional Similarity metric (derived from Gene Ontology annoatation graph) within groups of genes given a list of gene GO annoatations.
[depend on AIGO Python package https://pypi.python.org/pypi/AIGO/0.1.0]
