# agrogenom
This github repository explains the procedure and makes available the required software for the construction of a phylogenomic database from a (bacterial) pangenome dataset, and its use for the inference of genome histories, including the estimation of events of gene transfer, duplication and loss and the reconstruction of ancestral genome gene contents. The specificity of this procedure is to trace the co-evolution of genes across genomes through events of co-transfer and co-duplications, modelled as 'block events'. 

The application of this procedure to a dataset of 47 Rhizobiaceae genomes is presented in the manuscript [Lassalle et al. 2016], and results are interactively browsable on the website of the [Agrogenom database](http://phylariane.univ-lyon1.fr/db/agrogenom/3/).

The [pipeline/] section decribes extensively the procedure for generating the database and performing the ancestral genome reconstruction, whereas the [scripts/] section contains the main programs (in which are implemented original algorithms) that were adapted towards further use and will be maintained.

Authors: Florent Lassalle, Rémi Planel.
contact: [f.lassalle@imperial.ac.uk](mailto:f.lassalle@imperial.ac.uk)

[Lassalle et al. 2016]: http://biorxiv.org/content/early/2016/10/20/034843
[pipeline/]: https://github.com/flass/agrogenom/blob/master/pipeline/
[scripts/]: https://github.com/flass/agrogenom/blob/master/scripts/
