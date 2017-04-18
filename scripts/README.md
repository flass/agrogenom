Here are maintained the main Python modules and scripts for the construction of [Agrogenom database](http://phylariane.univ-lyon1.fr/db/agrogenom/3/). The full pipeline for construction of the database is described in the [pipeline section](https://github.com/flass/agrogenom/pipeline/).

['rec_to_db.py'](https://github.com/flass/agrogenom/rec_to_db.py) is an API for the loading and construction of the phylogenomic database Agrogenom from a collection of reconciled gene trees.
[depend on tree2 Python package https://github.com/flass/tree2]

['blockevents.py'](https://github.com/flass/agrogenom/blockevents.py) (library) and ['getblockevents.py'](https://github.com/flass/agrogenom/getblockevents.py) (execution script) are used to build block evolutionary events given a collection of reconciled gene trees and a phylogenomic database documenting gene neighbourhoods.
[depend on tree2 Python package https://github.com/flass/tree2]

['score_genegroup_funsim.py'](https://github.com/flass/agrogenom/score_genegroup_funsim.py) is used to compute Functional Homogeneity metric (derived from Gene Ontology annoatation graph) within groups of genes given a list of gene GO annoatations.
[depend on AIGO Python package https://pypi.python.org/pypi/AIGO/0.1.0]
