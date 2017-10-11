## Introduction

Here is described the programatic details of the procedure for generating and analysing the Agrogenom phylogenomic database, as described in [Lassalle et al. 2016].  

This software pipeline is semi-automatized, i.e. there are scripts and routines that have been developped and made (fairly) flexible for generic usage, but it cannot be deployed in one click or command, because of the many idiosyncracies of everyone's project.  

Typically, the way input data - from genome sequence and annotation up to homologous gene family trees - is generated may depend on contingent factors such as the technology that generated the sequences, the version of annotation formats, etc., but also on someone's taste (e.g. ML vs. bayesian phylogenetic trees).  

Also, the way parallel calculation is implemented is very specific to anyone's computational environment - e.g. if using a computer cluster, what kind of job scheduler system is used? - and should be adapted accordingly. This should be fairly straightforward, given the parallelization relies on the fondamentally parallel structure of the data, typically the many gene families are considered independent for several long computational steps like sequence alignement, gene tree inference and gene tree/species tree reconciliation - in fact only the last step of block event reconstruction consider the gene families jointly.  

Thus, one should consider the following code (gatherred in [pipeline_simple.sh]) as a template manual to perform inferences and analyses as described in [Lassalle et al. 2016]. To exactly replicate the work from this study please refer to the more detailed script [pipeline_agrogenom.csh], which is also  annotated with details on intermediate results of that particular study.  

The pipeline is divided in three parts:
  
1. Homologous gene tree database  
*the generation of the genome database and sequence alignements and phylogenetic trees for each gene family*
2. Species tree/gene tree reconciliations  
*the modification and annotation of gene trees to obtain scenarios of gene family evolution over the species tree, involving gene duplication, transfer and loss (co-)events*
3. Functional analysis of genome histories  
*the study of distribution of functional annotation of gene sets of interest, namely clade-specific gene sets*


## I. Homologous gene tree database

The procedure is the same than for [HOGENOM databases], which is described in the [Penel et al. 2009] paper and with updated details of the database releases [here](http://doua.prabi.fr/databases/hogenom/home.php?contents=methods).  
For that reason, the procedure will not be described on this page. For the details of the procedure that was used in [Lassalle et al. 2016], please refer to the [pipeline_agrogenom.csh] script.  
The source dataset and intermediary result files can be found at the corresponding [Figshare project] (items).

## II. Species tree/gene tree reconciliations

In this section, we will infer evolutionary events that affect gene families: gene duplication, horizontal gene transfer (HGT or just transfer) and gene loss, i.e. DTL events.  
We will first reconstruct the evolutionary scenarios for each gene family independently, and then consider them jointly by infering co-events that involve blocks of several contiguous genes.  
This reconstruction is performed via the reconciliation of gene trees with the species trees, and lead to the inference of the gene contents of genome of ancestral species (putative ancestors of the extant species placed at the nodes of the species phylogeny).  

For this, we will rely on a step-wise procedure for the gradual reconciliation of gene trees based on a series of criteria based on (supported) topological conflict with the species tree and species distribution in clades:

Let's assume we start with a rooted gene tree. (if not, one can use [TPMS] to root it with species tree-aware criteria, cf. [subsection II.1](https://github.com/flass/agrogenom/tree/master/pipeline#1-root-gene-trees))

![fig0]

First, putative paralogous lineages are identified based on species multiplicity in clades. Within those clades, gene tree topology updates (SPR moves that do not disrupt well-supported branches) were attempted, and retained when they decreased the incidence of duplications. For remaining paralogous lineage, the corresponding subtrees are extracted from the full gene tree. (cf. [subsection II.2](https://github.com/flass/agrogenom/tree/master/pipeline#2-find-putative-duplications-and-extract-unicopy-subtrees))

![fig1]

Within each of these lineages, statistically supported topological conflict is then resolved by inferring HGT events, using [Prunier], a parsimony-based algorithm [Abby et al. 2010]. (cf. [subsection I.3](https://github.com/flass/agrogenom/tree/master/pipeline#3-detect-hgt-based-on-statistically-suported-topological-conflict))

![fig2] 

Histories of all lineages are then mapped back to the full gene family tree and integrated as one gene family scenario. (cf. [subsection II.4](https://github.com/flass/agrogenom/tree/master/pipeline#4-integration-of-reconciliation-scenarios))

![fig3]  

At this point though, several alternative (possibly conflicting) reconciliation scenarios generated by independent inferences on overlapping lineage subtrees ("replicates") are all recorded as potential solutions.  

Then, block events - unique DTL events that are thought to have involved blocks of several genes - are inferred across gene families, by searching similar events (that took place in the same ancestors in the species tree) in trees of genes that are contiguous in extant genomes.  (cf. [subsection II.5](https://github.com/flass/agrogenom/tree/master/pipeline#5-build-genomic-blocks-of-evolutionary-events))

![fig4] 

The size of these block events is then used as a criterion to evaluate which of the multiple solutions from different replicates are optimal, assuming that the larger the blocks of similar events, the more likely the component events are to be rightly estimated, i.e. the inferred single events are "confirmed" by their many similar neighbours. This optimal choice of strongly supported transfer events will then be used as a backbone for reconciliations, to which more events will be added in order to obtain scenarios that are more parsimonious in duplication and losses.

Due to the previous reliance on strong statistical support for topological conflict to characterize HGT, many potential transfer events may have come unrecognize yet, due to their not inducing gene/species tree topological conflict, or it not occuring on well-supported branches of the gene tree. Failiure to recognize those events may induce the inferrence of many spurious ancient duplication and loss events.  
We thus use the algorithm for detection of taxonomic incongruencies described in [Bigot et al. 2013] (re-implemented in [rec_to_db.py] module) to identify subtrees with abnormal combination of species; transfer events are then inferred within these subtree when they avoid a too high number of duplication and loss events to be induced otherwise this conflictis first used for the inferrence of HGT events. Finally, duplication and losses are inferred where necessary to obtain consistent reconciliations, i.e. have no unexplained topological incongruency left between gene and species trees. 

![fig5] 

Once optimal transfer and duplication events are fixed, a last step of completion of reconciliations is then performed within orthologous lineages (deriving from any duplication or "additive" transfer event), to detect transfer from the distribution of species in the lineage, where heterogeneous occurrence pattern of the gene in the species tree would otherwise have to be explained by complicated patterns of ancient duplication and covergent losses. 

![fig6] 

Block events are ultimately recomputed, to capture co-events in this final set of DTL events.  

The pangenome-wide history of genomes is now reconstructed! Now it's time to study and to interpret the patterns of genome evolution...

![figspetreegenome]

-----------

Below is the simplified version of code supporting this procedure:

### 0. Set-up relational database

This pipeline requires setting up a PostgreSQL database. For efficiency reasons (use of many dynamic queries at certain stages), it is advised to use it locally, i.e. with `localhost` as host server, but using a host on a local network will do as well.

```bash
sudo apt-get install postgresql postgresql-client postgresql-contrib postgresql-doc python-psycopg2 
# then create a database (see https://www.postgresql.org/docs/9.6/static/tutorial-createdb.html for troubleshooting)
sudo createdb yourdbname
# and load in it the agrogenom schema
psql -h yourservername -U yourusername -d yourdbname < scripts_agrogenom/agrogenomdb_schema.sql
```

### 1. Root gene trees

The reconciliation of gene and species tree requires both trees to be rooted and highly depends on where the root are placed. Usual approaches for rooting based on balancing the tree topology (e.g. using midpoint rooting) can be inappropriate in case of certain lineages evolving at different rate.  
In our context, we will favour rooting that already consider the gene tree together with the species tree, in order to minimize the impact of the rooting on the subsequently inferred scenario, i.e. to find a root that is parsimonious in terms of implied DTL events. This is achieved by minimizing a score that combines to criteria: the first one aims at maximizing the size of subtrees with unicopy sequences while the second one tries to maximize the accuracy of taxonomic affectations for nodes [Bigot et al. 2013].
(for this step, it can be handy to refer to [TMPS manual][TPMS])

```bash
mkdir ./tpms_db/
# prepare a TPMS tree collection database from a folder containing all the pangenome's gene trees
tpms_mkdb -extract-seqnames -families-trees=./phyml_trees -output=./tpms_db/seqNames2spe
# fill in the sequence list template
python scripts_agrogenom/fill_tpms_seqlist.py ./tpms_db/seqNames2spe ./tpms_db/seqNames2species
# have to manually create a modified species tree file bearing the label 'ROOT:0' at the root (both are in Newick format)
sed -e 's/;/ROOT:0;/' ./reftree > ./reftree.tpms
# create database
tpms_mkdb -families-trees=phyml_trees -seq-to-species=./tpms_db/seqNames2species -sp-tree=./reftree.tpms -output=./tpms_db/genetrees.unrooted
# root gene trees
tpms_computations -collection=./tpms_db/genetrees.unrooted -output-dir=./tpms_db/ -threads=8 << EOF
RC
S
EOF
mv ./tpms_db/Collection ./tpms_db/genetrees.rooted
# split the result into unique tree files
mkdir ./rooted_trees/
python scripts_agrogenom/split_tpms_db.py ./tpms_db/genetrees.rooted ./rooted_trees
```
Because TPMS roots poorly trees in absence of duplication and presence of many transfers, one can use alternative method: find tranfers in unicopy gene trees using [Prunier], testing every root, to then re-root them consistently to the parsimonious transfer scenario.

```bash
# identify unicopy gene families from nucleic alignments 
python scripts_agrogenom/get_unicopy_fams.py ./nuc_alns ./unicopy_nuc_aln_list

# prunier detection of transfers
mkdir unicopy_prunier_out
# use a reference species tree devoid of branch lengths and internal labels or node supports
python -c "import tree2 ; t = tree2.Node(file='reftree') ; t.write_newick('reftree.prunier', branch_lengths=False, ignoreBS=True)"
# [TO ADAPT FOR COMPUTATION IN PARALLEL]
for famaln in `cat unicopy_nuc_aln_list` ; do
fam=${famaln%%.*}
famgt=phyml_trees/$fam.nwk
echo $famgt >> unicopy_gene_tree_list
PrunierFWD_linux64 input.tree.file=reftree.prunier aln.file=$famaln genetree.file=$famgt sequence.type=dna fwd.depth=2 aln.type=FASTA boot.thresh.conflict=.90 max.bp=1.00 multi_root.bool=false > unicopy_prunier_out/$fam.prout
done
```
Note `multi.root.bool=false` refers to trying the reconciliation with an unique rooting of the *species tree*, which is assumed to be known; Prunier always tries all the roots of the unicopy gene tree, and returns the reconciliation given the most parsmonious rooting, hence our use of this reconciliation program here for rooting purposes.

```bash
# rooting of trees as by Prunier
mkdir unicopy_rooted_trees
python scripts_agrogenom/root_as_prunier.py ./unicopy_gene_tree_list ./reftree ./unicopy_prunier_out ./unicopy_rooted_trees
# replace TMPS-rooted unicopy gene tres by the Prunier-rooted ones
cp -f ./unicopy_rooted_trees/* ./rooted_trees/
```
Find example results for Agrogenom study on [Figshare](https://doi.org/10.6084/m9.figshare.4907066).

### 2. Find putative duplications and extract unicopy subtrees

Here the script finds potentially duplication nodes in gene trees based on the unicity criterion described in [Bigot et al. 2013] (reimplemented in [tree2 module][tree2] and used in [find_ancestral_duplications.py] script).

When the support for the node is below a threshold (in this exemple 0.9) and the number of losses implied by this duplication (and the leaf set below the node) is higher than a threshold (here 0.2, in number of implied loss per branches in the duplicated lineage), an attempt to modify the local tree topology is made. 
Targeted subtree pruning and regrafting (SPR) moves are proposed (without disrupting nodes with strong support - in this exemple above or equal to 0.9) in order to merge duplication events or to bring the events closer to the tips, so that duplication events can be inferred more recently in the species tree, which imply less subsequent loss events.
If a modification of the topology has been made, the branch length and supports of the new subtree are then re-evaluated with PhyML.

```bash
mkdir ./duplications
# [TO ADAPT FOR COMPUTATION IN PARALLEL]
for famgt in `ls phyml_trees/*` ; do
python scripts_agrogenom/find_ancestral_duplications.py u $famgt ./nuc_alns ./duplications reftree 0.9 0.2
done
```
As an output of running this script, one obtains a set of unicopy subtrees extracted from th gene tree. 
These subtrees are generated by pruning branches from the gene tree so that at most one sequence per species is present. 
The pruning is guided by the annotation of duplication nodes, and tries to sample different (possibly overlapping) combinations of leaves in presence of in-paralogs, i.e. when some species are duplicated but not all (cf. [fig2]).

Find example results for Agrogenom study and detail of output files on [Figshare](https://doi.org/10.6084/m9.figshare.4907381).

### 3. Detect HGT based on statistically suported topological conflict

The program [Prunier] is run on these unicopy subtree, to detect transfer events based on statistically suported topological conflict with the reference species tree.
Prunier outputs a list of gene tree pruning steps (analogous to the SPR moves in topology search, without regrafting) that generates a forest of subtrees that are in statistical agreement with the reference species tree ([Abby et al. 2010]), i.e. that there is no topological conflict that remains in the gene tree that involves node with support higher than a threshold (here 0.9, on a scale of 0 to 1.0).

```bash
# generate unicopy sub-family alignments from previous family alignments
mkdir ./duplications/subalns_fasta
python scripts_agrogenom/extract_subaln.py ./duplications/subtree2leaves ./nuc_alns_fasta fasta ./duplications/subalns_fasta

# perform horizontal transfer detection by Prunier on unicopy subtrees (06/07/12)
mkdir -p prunier_out/
# [TO ADAPT FOR COMPUTATION IN PARALLEL]
for subfamaln in `ls duplications/subalns_fasta/*` ; do
subfam=${subfamaln%.*}
subfamgt=./duplications/subtrees/$subfam.nwk
PrunierFWD_linux64 input.tree.file=reftree.prunier aln.file=$subfamaln genetree.file=$subfamgt sequence.type=dna fwd.depth=2 aln.type=FASTA boot.thresh.conflict=.90 max.bp=1.00 multi_root.bool=false > prunier_out/$subfam.prout
done
```
Find example results for Agrogenom study on [Figshare](https://doi.org/10.6084/m9.figshare.4910177).

### 4. Integration of reconciliation scenarios

Prunier outputs a forest of pruned subtree, which needs to be interpreted in terms of HGT events matched to a gene tree node, and with donor and recipient located in the species tree. 
This translation of Prunier output is done for every unicopy subtree; those many independent reconciliations are then integrated using [rec_to_db.py] script, by creating unique records of events and merging those that are redundant, thus preparing the integration of reconciliations into a global scenario for the whole gene family. 

```bash
# it is possible for parallel jobs to run the script: a stack of tasks will be distributed to parallel workers via dynamic query of the database.
# generate task list and corresponding db table
scripts_agrogenom/lsfullpath.py ./duplications/modified_trees > ./duplications/modified_tree_list
psql -h yourservername -U yourusername -d yourdbname << EOF
CREATE TABLE IF NOT EXISTS reconciliation_to_db (gene_tree_file_path varchar(500) PRIMARY KEY, current_status int DEFAULT 0, date timestamp with time zone DEFAULT now(), job_id DEFAULT 0);
\copy reconciliation_to_db (gene_tree_file_path) from './duplications/modified_tree_list'
UPDATE reconciliation_to_db SET current_status=0, date=now(), job_id=0;
EOF
# [TO ADAPT FOR COMPUTATION IN PARALLEL]
mkdir ./integrate_rec
for i in `seq 1 20` ; do
# every job will query a task table in the SQL db to fetch a new task; concurrent access of jobs to task table is avoided with explicit table lock
python rec_to_db.py --task.survey.table=public.reconciliation_to_db --dir.subtrees=./duplications/subtrees --dir.subtree2leaf.dict=./duplications/subtree2leaves --dir.prunier.out=./prunier_out --reference.tree=./reftree --jobid=$i --dir.output=./integrate_rec -b 0.9 &
done
```

This information will be then uploaded to a relational (SQL) database that gathers data relating to genome annotation, gene content and order and gene family scenarios of evolution (see [database/](https://github.com/flass/agrogenom/tree/master/pipeline/database) section and [pipeline_agrogenom.csh] for setting up the SQL db). 
It is important to note that at this point the plurality of solutions that may have been found in independent replicates is preserved.

```bash
# import of data (with interactive validation of commits)
# copy dumps created by multiple jobs
python scripts_agrogenom/import_sql_dumps.py -c -i -d --job.dumps.dir=./integrate_rec
```
Find example results for Agrogenom study and detail of output files on [Figshare](https://doi.org/10.6084/m9.figshare.4910210).

### 5. Build genomic blocks of evolutionary events

The script [getblockevents.py] finds block events and create the corresponding objects in the database.

It first explores the contemporary genomes (replicon by replicon), recognizing events on their lineage that are similar between neighbouring genes and aggregate them into "leaf block events". This search is performed with a greedy algorithm covering the replicon; "gap" genes without signal for a congruent event are first allowed, and then tested for compatibility and eventually rejected if negative (see a schematic representation of the [leaf block building algorithm][figleafblockalgo]). 


```bash
# similarly to above, parallel jobs can be working on independent tasks, distributed from a task table in the SQL db
psql -h yourservername -U yourusername -d yourdbname << EOF
CREATE TABLE public.buildingblocks AS (SELECT chromosome, count(*) as nb_genes, 0 as job_id, 0 as current_status, now() as date FROM genome.gene GROUP BY chromosome);
EOF

# leaf blocks reconstruction
mkdir ./blockevents 
# [TO ADAPT FOR COMPUTATION IN PARALLEL]
for i in `seq 1 20` ; do
# every job will query a task table in the SQL db to fetch a new task; concurrent access of jobs to task table is avoided with explicit table lock
scripts_agrogenom/getblockevents.py ./integrate_rec/reconciled_tree_pickles reftree ./blockevents $i -d yourdbname -c public.buildingblocks l
done
```

In a second step, leaf block events from different genomes referring to the same evolutionary (duplication, transfer) events are aggregated into an "ancestral block event". Because the compatibility of member events is not transitive across connected leaf blocks, some leaf blocks must be split and recomposed to reach a consistent assignation of leaf blocks to ancestral blocks (see a schematic representation of the [ancestral block building algorithm][figancblockalgo]). 

```bash
# this is done in one sequential job
i=21
scripts_agrogenom/getblockevents.py ./integrate_rec/reconciled_tree_pickles reftree ./blockevents $i -d yourdbname -c public.buildingblocks a
```

Ancestral blocks are thus the units of actual evotutionary events that gather homologous single gene events set in the past (i.e. located at internal gene tree nodes) - or rather, gathering independent observations in different gene families of the same evolutionary event.  
Leaf blocks are the descendent of the ancestral blocks, i.e. their realization in contemporary genomes.

Find example results for Agrogenom study and detail of output files on [Figshare](https://doi.org/10.6084/m9.figshare.4910576).

### 6. Refinement and completion of reconciliations

When multiple transfer events have been proposed by independent Prunier inferences on overlapping unicopy subtrees ( = replicates)
for a same non-redundant gene tree node, it is possible to discriminate the best reconciliation among all inferences using the regional signal, i.e. block event information.  
The events that are found in the largest blocks of genes are retained as optimal. When this criterion is not discriminant (most events involve single genes), the (equivalent) transfer event(s) with the highest frequency amongst the independent replicates is chosen as optimal.

This choice of optimal events may create the cohabitation in the same gene tree of events inferred in different replicates, which does not guaratee the consistency of the global scenario.  
This potential conflict is resolved by an a posteriori top-down parcours of the gene tree where events are selected from the event pool so to ensure that the species tree assignation of the gene tree node is coherent with its parent's assignation.

Additional HGT are then inferred, based on the recognition of taxonomic incongruence (using the scoring algorithm from [Bigot et al. 2013] reimplemented in [tree2 module][tree2]). This algorithm computes the most recent common ancestor of species represented in a gene clade and assigns this ancestral taxon to the node; a large increase in taxonomic depth (i.e. closeness to species tree root) from a node to its father node is seen as an incongruence, than can be interpreted as evidence for HGT event(s) occurred in the subtree. "Alien" parts of the subtree generating the excess taxonomic depth are identified and transfer events are inferred at their stem.

```bash
mkdir refined_trees
scripts_agrogenom/refine_transfer_annotation.py -d yourdbname duplications/modified_trees reftre refined_trees

# import refined scenarios into database
psql -h yourservername -U yourusername -d yourdbname << EOF
\copy phylogeny.event from './refined_trees/event_dump.tab'
\copy phylogeny.event_possible_species_node from './refined_trees/event_possible_species_node_dump.tab'
INSERT INTO phylogeny.reconciliation_collection VALUES ('reco_col_1', '/full/path/to/refined_trees');
\copy phylogeny.represented_event FROM './refined_trees/represented_event_dump.tab'
EOF
```

Find example results for Agrogenom study and detail of output files on [Figshare](https://doi.org/10.6084/m9.figshare.4910699).

At this point, it is assumed that we identified all events that lead to the creation of a group of  orthologous genes, i.e. any duplication or additive HGT that created respectively paralogous or xenologous lineages.

### 7. Ancestral content reconstruction

Gene families are thus separated into orthologous subfamilies in order to reconstruct the presence/absence state of each ortholog in ancestral genomes.
For this, we use the Wagner parsimony algorithm implemented in [COUNT] program ([Csuros 2008]), which is called by [ancestral_content.py] script.
Doing so, a final round of reconciliation is made, inferring HGT events to explain heterogeneous representation of species within othologuous gene subfamilies.

```bash
mkdir -p ./ancestral_content/tmp
scripts_agrogenom/lsfullpath.py ./nuc_alns > ./nucfam_fasta_list
python scripts_agrogenom/ancestral_content.py -d yourdbname --infer='Count.AsymmetricWagner' -o '-gain 1.5' --tmp-dir=./ancestral_content/tmp --output-subfamily-trees=true --output-pickled-family-histories=true ./nucfam_fasta_list ./refined_trees/reconciled_tree_pickles ./reftree ./ancestral_content
```

The script [ancestral_content.py] generates new trees (with additional transfers), which need to be gatherred with previously calculated trees into a single collection (used notably for generating the [Agrogenom] web interface).

```bash
# merge phyloXML file collections
mkdir ./reconciliation_collection
cp -r ./modified_trees/* ./reconciliation_collection/
cp -np ./refined_trees/normed_phyloxml_trees/*.xml ./reconciliation_collection/normed_phyloxml_trees/
cp -np ./refined_trees/reconciled_tree_pickles/*.pickle ./reconciliation_collection/reconciled_tree_pickles/

# export data as a table dump
python scripts_agrogenom/dump_reconciliation_collection_tables.py -c ./reconciliation_collection/ -r reftree -i 'reco_col_1' -t -a

# final storage the reconciliation collection into database
psql -h yourservername -U yourusername -d yourdbname << EOF
\copy phylogeny.event from './ancestral_content/Count.AsymmetricWagner.gain_15/reconciliation_collection/event_dump.tab'
\copy phylogeny.event_possible_species_node from './ancestral_content/Count.AsymmetricWagner.gain_15/reconciliation_collection/event_possible_species_node_dump.tab'
\copy phylogeny.reconciliation_collection FROM './ancestral_content/Count.AsymmetricWagner.gain_15/reconciliation_collection/reconciliation_collection_dump.tab'
\copy phylogeny.represented_event FROM './ancestral_content/Count.AsymmetricWagner.gain_15/reconciliation_collection/represented_event_dump.tab'
EOF
```
Find example results for Agrogenom study and detail of output files on [Figshare](https://doi.org/10.6084/m9.figshare.4910702).

### 8. Finalizing reconciliation and synthesis

Now to finalize the phylogenomic database, we want to repeat the block event reconstruction with the final set of DTL events.
Thus we need to clean the database from previous block event records:

```bash
# connect to database
psql -h yourservername -U yourusername -d yourdbname
```
```sql
CREATE TABLE public.buildingblocks AS (SELECT chromosome, count(*) as nb_genes, 0 as job_id, 0 as current_status, now() as date FROM genome.gene GROUP BY chromosome);
UPDATE public.buildingblocks SET job_id=0, current_status=0, date=now();
-- empties the pre-existing block data
TRUNCATE blocks.block_possible_species_node;
TRUNCATE blocks.leaf_block CASCADE ;
TRUNCATE blocks.fam2block ;
-- CAUTION when making TRUNCATE ... CASCADE on blocks.ancestral_block  !!!
ALTER TABLE phylogeny.event DROP CONSTRAINT event_anc_block_id_fkey ;
UPDATE phylogeny.event SET anc_block_id = NULL ;
TRUNCATE blocks.ancestral_block CASCADE ;
ALTER TABLE phylogeny.event ADD FOREIGN KEY (anc_block_id) REFERENCES blocks.ancestral_block (anc_block_id) ON DELETE RESTRICT ;
```
and generate new block events:

```bash
# sequential run (could be parallel, but not as many events to consider now)
# leaf blocks reconstruction
python scripts_agrogenom/getblockevents.py ./reconciliation_collection/reconciled_tree_pickles reftree ./blockevents 99 -c public.buildingblocks 
# ancestral blocks reconstruction !!! high memory use
python scripts_agrogenom/getblockevents.py ./reconciliation_collection/reconciled_tree_pickles reftree ./blockevents 100 -a -v -c public.buildingblocks 
# note job ids 99 and 100 are arbitrary and the only important thingis that they are different
# import into db
psql -h yourservername -U yourusername -d yourdbname < ./blockevents/getblockevent_out.a.100/blocks_db_dump/allblocktables_sqldump.sql
```
You can now generate synthesis files, including annotated Newick tree files, for rapid vizualization in a tree viewer like [seaview](http://doua.prabi.fr/software/seaview) (available as a package in most Debian distributions) of a particular gene family history, or of the genome-wide synthesis of geneome gene content evolution.
Output also includes tables for presence/absenceprofiles of all gene families (and most interestingly, ortholous subfamilies) at all species tree's node, ancestors and extant. There are also files listing clade-specific genes - orthologous subfamilies specifically gained by an ancestor and conserved by all its descendants.

```bash
# generate synthesis files: genome history and clade-specific gene tables
python scripts_agrogenom/ancestral_content.py -d phylariane -p $dbpwd --input-family-history-pickles=$anccont/family_history_pickles --input-subfamily-trees=$anccont/ortho_subtrees --output-all --output-subfamily-trees=false --output-pickled-family-histories=false ./nucfam_fasta_list ./refined_trees/reconciled_tree_pickles reftree ./synthesis
```

Find example results for Agrogenom study and detail of output files on [Figshare](https://doi.org/10.6084/m9.figshare.4924439).  

You can also export the genome-wide synthesis of genome gene contents as a nice SVG output graphics!

```bash
python scripts_agrogenom/draw_genome_contents+events.py ./synthesis/genome_synthesis.pickle
```

Find example of rendering for Agrogenom study on [Figshare](https://doi.org/10.6084/m9.figshare.4924433).  

------------------

That's all, done!

Here is a graphical summary of that (long!) pipeline for ancestral genome reconstruction (see [here][figallpdf] for higher resolution):

![figall]

## III. Functional analysis of genome histories  

An exemple of application of this pipeline is presented in the [Lassalle et al. 2016] study, where we introduce the statistical testing of the preferential co-transfer and conservation of genes with more related biochemical functions.
The bioinformatics methods involved there depend quite strongly on the specific datset, notably the sources of functional annotation that were used, so no generic pipeline is detailed here; for the specific use of the program [score_genegroup_funsim] for the [Lassalle et al. 2016] study,  please see command details in the corresponding section of the script [pipeline_agrogenom.csh].



### References:
 [Lassalle F et al. (2016)][Lassalle et al. 2016] "Ancestral genome reconstruction reveals the history of ecological diversification in Agrobacterium.", bioRxiv, p. 034843. doi: 10.1101/034843.  
 [Penel S et al. (2009)][Penel et al. 2009] "Databases of homologous gene families for comparative genomics" BMC Bioinformatics, 10(S6):S3.  
 [Bigot T et al. (2013)][Bigot et al. 2013] "TPMS: a set of utilities for querying collections of gene trees", BMC Bioinformatics 14:109. doi: 10.1186/1471-2105-14-109  
 [Abby SS et al. (2010)][Abby et al. 2010] "Detecting lateral gene transfers by statistical reconciliation of phylogenetic forests". 11:324-324. doi: 10.1186/1471-2105-11-324.  
 [Csűrös M (2008)][Csuros 2008] "Ancestral Reconstruction by Asymmetric Wagner Parsimony over Continuous Characters and Squared Parsimony over Distributions", in Nelson, C. E. and Vialette, S. (eds) Comparative Genomics. Springer Berlin Heidelberg (Lecture Notes in Computer Science, 5267), pp. 72–86.  

[Lassalle et al. 2016]: http://biorxiv.org/content/early/2016/10/20/034843
[HOGENOM databases]: http://doua.prabi.fr/databases/hogenom/home.php
[Penel et al. 2009]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-S6-S3
[Prunier]: http://pbil.univ-lyon1.fr/software/prunier/
[TPMS]: https://github.com/tbigot/tpms
[Bigot et al. 2013]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-109
[Abby et al. 2010]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-324
[Agrogenom]: http://phylariane.univ-lyon1.fr/db/agrogenom/3/
[Csuros 2008]: http://link.springer.com/chapter/10.1007/978-3-540-87989-3_6
[COUNT]: http://www.iro.umontreal.ca/~csuros/gene_content/count.html

[pipeline_agrogenom.csh]: https://github.com/flass/agrogenom/blob/master/pipeline/pipeline_agrogenom.csh
[pipeline_simple.sh]: https://github.com/flass/agrogenom/blob/master/pipeline/pipeline_simple.sh
[rec_to_db.py]: https://github.com/flass/agrogenom/blob/master/scripts/rec_to_db.py
[find_ancestral_duplications.py]: https://github.com/flass/agrogenom/blob/master/scripts/find_ancestral_duplications.py
[getblockevents.py]: https://github.com/flass/agrogenom/blob/master/scripts/getblockevents.py
[ancestral_content.py]: https://github.com/flass/agrogenom/blob/master/scripts/ancestral_content.py
[score_genegroup_funsim]: https://github.com/flass/agrogenom/blob/master/scripts/score_genegroup_funsim.py
[tree2]: https://github.com/flass/tree2
[Figshare project]: https://figshare.com/projects/Ancestral_genome_reconstruction_reveals_the_history_of_ecological_diversification_in_Agrobacterium/20894

[fig0]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline-0.png
[fig1]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline-1.png
[fig2]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline-2.png
[fig3]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline-3.png
[fig4]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline-4.png
[fig5]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline-5.png
[fig6]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline-6.png
[figspetreegenome]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/species_tree_map_events_gene_contents.png
[figall]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline_loop.png
[figallpdf]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/reconciliation_pipeline_loop.pdf
[figleafblockalgo]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/leaf_block_algorithm.pdf
[figancblockalgo]: https://github.com/flass/agrogenom/blob/master/pipeline/figures/anc_block_algorithm.pdf

