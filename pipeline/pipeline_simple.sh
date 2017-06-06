#!bash

# variables to define:
scripts_agrogenom=/path/to/scripts_agrogenom

cd /path/to/yourreconciliationfolder
# this folder must contain the following files/folders:
# - nuc_alns     : the nucleotidic alignments (preferentially aligned by codons) in FASTA format for all homologous gene families
# - phyml_trees/ : containing gene trees in Newick format from all gene homologous families with >= 3 gene sequences;
#                  it can be computed by any ML/bayesian tree inference software (e.g. PhyML) provided they carry branch supports
# - reftree      : the reference species tree in Newick format

### 0. Set-up relational database

# it is necesasry to create a database following the agrogenom schema
# for that, first set up a SQL database server
# set up a PostgreSQL database (advise using locally, i.e. with localhost as host server)
sudo apt-get install postgresql postgresql-client postgresql-contrib postgresql-doc python-psycopg2 
# then create a database (see https://www.postgresql.org/docs/9.6/static/tutorial-createdb.html for troubleshooting)
sudo createdb yourdbname
# and load in it the agrogenom schema 
# generate table yourself based on the model, or use the perl module in pipeline/perl_utils/generate_genome_schema_tables.tar.gz
psql -h yourservername -U yourusername -d yourdbname < $scripts_agrogenom/agrogenomdb_schema.sql


### 1. Root gene trees
## general case of multicopy gene families (i.e. with paralogous lineages)
mkdir ./tpms_db/
# prepare a TPMS tree collection database from a folder containing all the pangenome's gene trees
tpms_mkdb -extract-seqnames -families-trees=./phyml_trees -output=./tpms_db/seqNames2spe
# fill in the sequence list template
python $scripts_agrogenom/fill_tpms_seqlist.py ./tpms_db/seqNames2spe ./tpms_db/seqNames2species
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
python $scripts_agrogenom/split_tpms_db.py ./tpms_db/genetrees.rooted ./rooted_trees
## case of unicopy gene families (i.e. likely orthologs or xenologs)
# identify unicopy gene families from nucleic alignments 
python $scripts_agrogenom/get_unicopy_fams.py ./nuc_alns ./unicopy_nuc_aln_list
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
# rooting of trees as by Prunier
mkdir unicopy_rooted_trees
python $scripts_agrogenom/root_as_prunier.py ./unicopy_gene_tree_list ./reftree ./unicopy_prunier_out ./unicopy_rooted_trees
# replace TMPS-rooted unicopy gene tres by the Prunier-rooted ones
cp -f ./unicopy_rooted_trees/* ./rooted_trees/

### 2. Find putative duplications and extract unicopy subtrees

mkdir ./duplications
# [TO ADAPT FOR COMPUTATION IN PARALLEL]
for famgt in `ls phyml_trees/*` ; do
python $scripts_agrogenom/find_ancestral_duplications.py u $famgt ./nuc_alns ./duplications reftree 0.9 0.2
done

### 3. Detect HGT based on statistically suported topological conflict

# generate unicopy sub-family alignments from previous family alignments
mkdir ./duplications/subalns_fasta
python $scripts_agrogenom/extract_subaln.py ./duplications/subtree2leaves ./nuc_alns_fasta fasta ./duplications/subalns_fasta
# perform horizontal transfer detection by Prunier on unicopy subtrees (06/07/12)
mkdir -p prunier_out/
# [TO ADAPT FOR COMPUTATION IN PARALLEL]
for subfamaln in `ls duplications/subalns_fasta/*` ; do
subfam=${subfamaln%.*}
subfamgt=./duplications/subtrees/$subfam.nwk
PrunierFWD_linux64 input.tree.file=reftree.prunier aln.file=$subfamaln genetree.file=$subfamgt sequence.type=dna fwd.depth=2 aln.type=FASTA boot.thresh.conflict=.90 max.bp=1.00 multi_root.bool=false > prunier_out/$subfam.prout
done

### 4. Integration of reconciliation scenarios

# it is possible for parallel jobs to run the script: a stack of tasks will be distributed to parallel workers via dynamic query of the database.
# generate task list and corresponding db table
$scripts_agrogenom/lsfullpath.py ./duplications/modified_trees > ./duplications/modified_tree_list
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
# import of data (with interactive validation of commits)
# copy dumps created by multiple jobs
python $scripts_agrogenom/import_sql_dumps.py -c -i -d --job.dumps.dir=./integrate_rec

### 5. Build genomic blocks of evolutionary events

# similarly to above, parallel jobs can be working on independent tasks, distributed from a task table in the SQL db
psql -h yourservername -U yourusername -d yourdbname << EOF
CREATE TABLE public.buildingblocks AS (SELECT chromosome, count(*) as nb_genes, 0 as job_id, 0 as current_status, now() as date FROM genome.gene GROUP BY chromosome);
EOF
# leaf blocks reconstruction
mkdir ./blockevents 
# [TO ADAPT FOR COMPUTATION IN PARALLEL]
for i in `seq 1 20` ; do
# every job will query a task table in the SQL db to fetch a new task; concurrent access of jobs to task table is avoided with explicit table lock
$scripts_agrogenom/getblockevents.py ./integrate_rec/reconciled_tree_pickles reftree ./blockevents $i -d yourdbname -c public.buildingblocks l
done
# this is done in one sequential job
i=21
$scripts_agrogenom/getblockevents.py ./integrate_rec/reconciled_tree_pickles reftree ./blockevents $i -d yourdbname -c public.buildingblocks a

### 6. Refinement and completion of reconciliations

mkdir refined_trees
scripts_agrogenom/refine_transfer_annotation.py -d yourdbname duplications/modified_trees reftre refined_trees

# import refined scenarios into database
psql -h yourservername -U yourusername -d yourdbname << EOF
\copy phylogeny.event from './refined_trees/event_dump.tab'
\copy phylogeny.event_possible_species_node from './refined_trees/event_possible_species_node_dump.tab'
INSERT INTO phylogeny.reconciliation_collection VALUES ('reco_col_1', '/full/path/to/refined_trees');
\copy phylogeny.represented_event FROM './refined_trees/represented_event_dump.tab'
EOF

### 7. Ancestral content reconstruction

# merge phyloXML file collections
mkdir ./reconciliation_collection
cp -r ./modified_trees/* ./reconciliation_collection/
cp -np ./refined_trees/normed_phyloxml_trees/*.xml ./reconciliation_collection/normed_phyloxml_trees/
cp -np ./refined_trees/reconciled_tree_pickles/*.pickle ./reconciliation_collection/reconciled_tree_pickles/

# export data as a table dump
python $scripts_agrogenom/dump_reconciliation_collection_tables.py -c ./reconciliation_collection/ -r reftree -i 'reco_col_1' -t -a

# final storage the reconciliation collection into database
psql -h yourservername -U yourusername -d yourdbname << EOF
\copy phylogeny.event from './ancestral_content/Count.AsymmetricWagner.gain_15/reconciliation_collection/event_dump.tab'
\copy phylogeny.event_possible_species_node from './ancestral_content/Count.AsymmetricWagner.gain_15/reconciliation_collection/event_possible_species_node_dump.tab'
\copy phylogeny.reconciliation_collection FROM './ancestral_content/Count.AsymmetricWagner.gain_15/reconciliation_collection/reconciliation_collection_dump.tab'
\copy phylogeny.represented_event FROM './ancestral_content/Count.AsymmetricWagner.gain_15/reconciliation_collection/represented_event_dump.tab'
EOF

### 8. Finalizing reconciliation and synthesis

psql -h yourservername -U yourusername -d yourdbname << EOF
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
EOF

# sequential run (could be parallel, but not as many events to consider now)
# leaf blocks reconstruction
python $scripts_agrogenom/getblockevents.py ./reconciliation_collection/reconciled_tree_pickles reftree ./blockevents 99 -c public.buildingblocks 
# ancestral blocks reconstruction !!! high memory use
python $scripts_agrogenom/getblockevents.py ./reconciliation_collection/reconciled_tree_pickles reftree ./blockevents 100 -a -v -c public.buildingblocks 
# note job ids 99 and 100 are arbitrary and the only important thingis that they are different
# import into db
ppython $scripts_agrogenom/draw_genome_contents+events.py ./synthesis/genome_synthesis.pickle
sql -h yourservername -U yourusername -d yourdbname < ./blockevents/getblockevent_out.a.100/blocks_db_dump/allblocktables_sqldump.sql

# generate synthesis files: genome history and clade-specific gene tables
python $scripts_agrogenom/ancestral_content.py -d phylariane -p $dbpwd --input-family-history-pickles=$anccont/family_history_pickles --input-subfamily-trees=$anccont/ortho_subtrees --output-all --output-subfamily-trees=false --output-pickled-family-histories=false ./nucfam_fasta_list ./refined_trees/reconciled_tree_pickles reftree ./synthesis


