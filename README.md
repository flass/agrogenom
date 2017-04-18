## Introduction

Here is described the programatic details of the procedure for generating and analysing the Agrogenom phylogenomic database, as described in [Lassalle et al. 2016].  
This software pipeline is semi-automatized, i.e. there are scripts and routines that have been developped and made (fairly) flexible for generic usage, but it cannot be deployed in one click or command, because of the many idiosyncracies of everyone's project.  
Typically, the way input data - from genome sequence and annotation up to homologous gene family trees - is generated may depend on contingent factors such as the technology that generated the sequences, the version of annotation formats, etc., but also on someone's taste (e.g. ML vs. bayesian phylogenetic trees).  
Also, the way parallel calculation is implemented is very specific to anyone's computational environment - e.g. if using a computer cluster, what kind of job scheduler system is used? - and should be adapted accordingly. This should be fairly straightforward, given the parallelization relies on the fondamentally parallel structure of the data, typically the many gene families are considered independent for several long computational steps like sequence alignement, gene tree inference and gene tree/species tree reconciliation - in fact only the last step of block event reconstruction consider the gene families jointly.  
Finally, one should consider the following code as a loose manual to perform inferences and analyses as described in [Lassalle et al. 2016]; to replicate the work from this study, please refer to the more detailed script [pipeline_agrogenom.csh] (which is also  more a blueprint for manual execution of the pipeline with steps to adapt than an all-in-one executable).  

The pipeline is divided in three parts:
  
1. Homologous gene tree database  
*this describes the generation of the sequence alignement and phylogenetic tree database*
2. Species tree/gene tree reconciliations
3. Ananalysis of genome histories



## I. Homologous gene tree database

The procedure is the same than for [HOGENOM databases], which is described in the [Penel et al. 2009] paper and with updated details of the database releases [here](http://doua.prabi.fr/databases/hogenom/home.php?contents=methods).  
For that reason, the procedure will not be described on this page. For the details of the procedure that was used in [Lassalle et al. 2016], please refer to the [pipeline_agrogenom.csh] script.  

## II. Species tree/gene tree reconciliations

In this section, we will infer evolutionary events that affect gene families: gene duplication, horizontal gene transfer (HGT or just transfer) and gene loss, i.e. DTL events.  
We will first reconstruct the evolutionary scenarios for each gene family independently, and then consider them jointly by infering co-events that involve blocks of several contiguous genes.  
This reconstruction is performed via the reconciliation of gene trees with the species trees, and lead to the inference of the gene contents of genome of ancestral species (putative ancestors of the extant species placed at the nodes of the species phylogeny).  

For this, we will rely on a step-wise procedure for the gradual reconciliation of gene trees based on a series of criteria based on (supported) topological conflict with the species tree and species distribution in clades.  
It first defines paralogous lineages in gene trees and extracts the corresponding subtrees from the ull gene tree. Within each of these lineages, statistically supported topological conflict is then resolved by inferring HGT events, using [Prunier], a parsimony-based algorithm. Histories of all lineages are then mapped back to the full gene family tree and integrated as one gene family scenario; this involves the additional inferrence of HGT events to resolve remaining topological conflict that would otherwise be explained by costly combinations of duplication and losses. At this point though, several alternative (possibly conflicting) reconciliation scenarios generated by independent inferences on overlapping lineage subtrees ("replicates") are all recorded as potential solutions.  
Then, block events - unique DTL events that are thought to have involved blocks of several genes - are inferred across gene families, by searching similar events (that took place in the same ancestors in the species tree) in trees of genes that are contiguous in extant genomes. The size of these block events is then used asa criterion to evaluate which of the multiple solutions from different replicates are optimal, assuming that the larger the blocks of similar events, the more likely the component events are to be rightly estimated, i.e. the inferred single events are "confirmed" by their many similar neighbours. With the optimal transfer events chosen, the remaining topological conflict is resolved by again inferring transfer or duplication events.  
Whith optimal transfer and duplication evenA last step of completion of reconciliations is then performed within orthologous lineages (deriving from any duplication or "additive" transfer event), to detect transfer from the distribution of species in the lineage, where heterogeneous ditribution of the gene in the species tree would otherwise have to be explained by complicated pattern sof ancient duplication and covergent losses. Block events are ultimately recomputed, to capture co-events in this final set of DTL events.  

### 1. Root gene trees

The reconciliation of gene and species tree requires both trees to be rooted and highly depends on where the root are placed. Usual approaches for rootig based on balancing the tree topology (e.g. using midpoint rooting) can be inappropriate in case of certain lineages evolving at different.  
In our context, we will favour rooting that already consider the gene tree together with the species tree, in order to minimize the impact of the rooting on the subsequently  inferred scenario, i.e. to find root that are parsimonious in terms of implied DTL events.  





### References:
Lassalle F et al.(2016) "Ancestral genome reconstruction reveals the history of ecological diversification in Agrobacterium.", bioRxiv, p. 034843. doi: 10.1101/034843.  
Penel S et al. (2009) "Databases of homologous gene families for comparative genomics" BMC Bioinformatics, 10 (Suppl 6):S3  

[Lassalle et al. 2016]: http://biorxiv.org/content/early/2016/01/13/034843
[pipeline_agrogenom.csh]: https://github.com/flass/agrogenom/blob/master/pipeline/pipeline_agrogenom.csh
[HOGENOM databases]: http://doua.prabi.fr/databases/hogenom/home.php
[Penel et al. 2009]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-S6-S3
[Prunier]: http://pbil.univ-lyon1.fr/software/prunier/
