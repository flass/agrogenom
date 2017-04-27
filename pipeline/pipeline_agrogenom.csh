#!/bin/csh

###################################
# I # Homologous gene tree database
###################################

# This part of the pipeline is derived from the HOGENOM database pipeline (Penel et al. 2009) 
# It can be made completely differently, as long as it delivers:
# - a bank of coding sequences (CDSs) grouped into homologous gene families (not just orthologous!)
# - nucleic acid alignments of those gene families
# - phylogenies derived from those alignements (gene trees), with branch supports.
# - a species phylogeny (species tree)
# (all trees need to be single trees, e.g. ML point estimates or consensuses of bayesian analyses)

 
################################
# I.1 # create sequence database
################################

## define environment variables
setenv genomes /path/to/data/storage/genomes
setenv agrogenom /path/to/data/storage/agrogenom
setenv acnucraw $agrodata/acnuc/raw
setenv acnucfam $agrodata/acnuc/hfxfams
setenv families $agrodata/families
setenv penel /panhome/penel/penel/svn/bin/Debian
setenv TPMS /panhome/lassalle/tpms/build
setenv dbname 49RHIZOB
setenv subdbname 47RHIZOB

# create directories
mkdir $genomes
mkdir $agrodata
mkdir $agrodata/acnuc
mkdir $agrodata/simon

mkdir $acnucraw
mkdir $acnucfam
mkdir $families

# SQL database credentials
setenv sqldbserver yourservername
setenv sqldbuser yourusername
setenv sqldbname yourdbname

## [TO DO YOURSELF / TO ADAPT TO YOUR PROJECT] 
# dowloaded 49 flat files of annotated genomic sequences in EMBL format (14/04/2012)
cd $acnucraw/flat_files/
# made availble on figshare.com
wget https://ndownloader.figshare.com/files/8252768

# set 5-letter code identifier (by matching taxid from UniProt speclist file or by creating a new one for new sequences)
cp pipeline/database/code5/$dbname.code $agrodata/code5/

## use of ACNUC database for (sub-)sequence query and retrieval ; references:
# Gouy M et al. 1985. ACNUC–a portable retrieval system for nucleic acid sequence databases: logical and physical designs and usage. Comput. Appl. Biosci. 1:167–172.
# Gouy M and Delmotte S. 2008. Remote access to ACNUC nucleotide and protein sequence databases at PBIL. Biochimie 90(4):555-562.
# can be dowloaded from http://doua.prabi.fr/databases/acnuc in `Local installation' section
# [TO DO YOURSELF] install ACNUC database system

## build an ACNUC database from genomic flat files (03/05/2012)
cd $agrodata
scripts/build_agrogenom_acnuc.csh $acnucraw >& $acnucraw/build.log 
# extracts protein sequences from annotated CDS
mkdir $agrodata/blastdb
scripts/extr_prot.csh $acnucraw $agrodata/blastdb/$dbname.prot.fasta
# 281223 extracted sequences.

# spot CDS having label without .EXT extension 
# (comes from accession bearing only one CDS ; they most often span the entire DNA molecule accession to which they belong and were predicted on basis of blastx hit in databases : most should not be true or complete CDS )
mkdir $agrodata/check_missed_annots
grep ">" $agrodata/blastdb/$dbname.prot.fasta | grep -v "\." | cut -d">" -f2 | cut -d" " -f1 > check_missed_annots/acc_without_subseqs
# 37 "false CDS"

## generate all vs. all protein BLAST matrix
# build blast db
cd $agrodata/blastdb
formatdb -i $dbname.prot.fasta
# [TO DO YOURSELF] here one need to implement his own routine for completing in parallel the (massive: 281223^2 ~ 8e+11) set of blast jobs.
# parameters are filtering on e-value <= 1e-5; tabular output, no header (old -m 8)
# agregate all the results into one file: $agrodata/families/$dbname.blastout
cd $agrodata
mkdir $families
foreach bo (`ls BLAST/blastout/res.1/`)
cat BLAST/blastout/res.1/$bo >> families/$dbname.blastout
end

## build protein sequence families with Silix/Hifix (04/05/2012) ; references:
# Miele V, Penel S, Duret L. 2011. Ultra-fast sequence clustering from similarity networks with SiLiX. BMC Bioinformatics. 12:116. doi: 10.1186/1471-2105-12-116.
# Miele V et al. 2012. High-quality sequence clustering guided by network topology and multiple alignment likelihood. Bioinformatics. 28(8):1078-1085. doi: 10.1093/bioinformatics/bts098.
cd $agrodata
cp blastdb/$dbname.prot.fasta $families/$dbname.fasta
$scripts/hifix.csh $families $dbname.fasta $dbname.blastout

##### 16/04/12
# number of protein families found by Hifix: 42171

# edit genome flat file to introduce gene family feature
# fetch all the molecule accesion names
cd $agrodata
setenv gcgacnuc $acnucraw/flat_files
setenv acnuc $acnucraw/index
query << EOF
term
n
select
sp=@
save
list1
$dbname.acc
stop
EOF
# produce a new flat file
mkdir $acnucfam/flat_files
mkdir $acnucfam/index
utils/add_fam_annot_genomes $agrodata/acnuc/$dbname.acc $families/$dbname.HFX.fnodes $acnucfam/flat_files/$dbname.dat
# 2086 extracted entries.
# 281223 CDS described.
# 281186 family written.
# -> match the 37 "false CDS"

# add NCBI taxon id if not present in 'FT /source' field of flat files
cp $acnucfam/flat_files/$dbname.dat acnucfam/flat_files/$dbname.dat_091112
$scritps/complete_taxid_in_embl.py $acnucfam/flat_files/$dbname.dat_091112 $agrodata/code5/dic_taxid_code $acnucfam/flat_files/$dbname.dat

# build a new ACNUC database
cd $agrodata
# require 'custom_qualifier_policy' file to specify that the gene family field must be indexed
# to be copied from 'pipeline/acnuc_utils/custom_qualifier_policy' into '$agrodata/acnuc/hfxfams/index/' directory
cat $agrodata/acnuc/hfxfams/index/custom_qualifier_policy
$scripts/build_agrogenom_acnuc.csh $acnucfam >& $acnucfam/build.log 


########################################
# I.2 # calculate gene family alignments
########################################

## retrieve sequences

# create one protein fasta file per gene family
mkdir $agrodata/protfam_fasta
$scripts/retrieve_fam_fasta.csh $acnucfam families/$dbname.flist protfam_fasta aa

# create one nucleotide fasta file per gene family
mkdir $agrodata/nucfam_fasta
$scripts/retrieve_fam_fasta.csh $acnucfam families/$dbname.flist nucfam_fasta nuc

# check for complete coverage of the families
cd $agrodata/check_missed_annots
ls protfam_fasta/ | cut -d'.' -f1 | sort -u > protfam_fasta_list.sortu
cut -f1 families/$dbname.HFX.fnodes | sort -u > $dbname.HFX.fnodes.sortu
diff $dbname.HFX.fnodes.sortu protfam_fasta_list.sortu | grep '<' | cut -d' ' -f2 > missed_annot_fams
foreach fam ( `cat missed_annot_fams` )
grep $fam families/$dbname.HFX.fnodes | cut -f2 >> missed_annot_cds
end
# compare the list of unannotated CDS to the list of those with accession names (i.e. missing the .EXT extension in their label) found in 'acc_without_subseqs'
sort -u acc_without_subseqs > acc_without_subseqs.sortu
sort -u missed_annot_cds > missed_annot_cds.sortu
diff acc_without_subseqs.sortu missed_annot_cds.sortu | grep '>'
# if nothing (only '<' marks in diff output) : OK
# prepares checks post alignments and get an idea of family size distribution
R --no-save << eof
famsizes = read.table("$families/49RHIZOB.HFX.fsizes", stringsAsFactors=FALSE)
distrsizes = sapply( list(famsizes$V2==1, famsizes$V2>=3, famsizes$V2>=4, famsizes$V2>=500, famsizes$V2>=1000, famsizes$V2>=2000), function(x){ length(which(x)) })
print(distrsizes)
toalign = famsizes$V1[famsizes$V2>=3]
length(toalign)
write(toalign, "fams_to_align")
eof
sort fams_to_align > fams_to_align.sort
##### 16/04/12
# number of gene families annotated in ACNUC db : 					42261
# number of ORFan gene families : 							27564
# number of family of size >= 3 (for which an alignment/tree will be computed) : 	10813	(listed in $agrodata/check_missed_annots/fams_to_align)
# number of family of size >= 4 (for which an alignment/tree will be reconciled) :	 9069
# number of family of size >= 500 : 	 						   19
# number of family of size >= 1000 : 	 						    4
# number of family of size >= 2000 : 	 						    2

## alignment of proteic families with MUSCLE
cd $agrodata
mkdir protfam_aln_muscle
mkdir logs/protfam_aln_muscle_logs
# [TO ADAPT] computation in parallel - here using a PBS/Torque queue submission system
$scripts/lsfullpath.py $PWD/protfam_fasta > protfam_fasta_list
# involve call to task manager utils/getnextbig: a stack of tasks is distributed to parallel workers;
# parallel jobs are cast with predefined walltime and as many as possible tasks are executed by each job.
$scripts/muscle.qsub $PWD/protfam_fasta_list $PWD/protfam_aln_muscle 1 500 $PWD/logs/protfam_aln_muscle_logs q1day

# muscle.csh does not make alignment for families < 3 seqs : 10813 alignment expected
ls protfam_aln_muscle | cut -d"." -f1 | sort > check_missed_annots/fams_aligned.sort
diff fams_to_align.sort check_missed_annots/fams_aligned.sort
# if 0 diff: OK.

## reverse translate the proteic alignments into codon alignments
# (for further dN/dS analysis and to avoid mis-alignment due to nucleic signal saturation) NB: create alignment file in Clustal format.
cd $agrodata
mkdir codonfam_aln_clustal
$scripts/lsfullpath.py $PWD/protfam_aln_muscle > protfam_aln_muscle_list
python $scripts/pal2nal.py $PWD/protfam_aln_muscle_list $PWD/nucfam_fasta $PWD/codonfam_aln_clustal
mv codonfam_aln_clustal/pal2nal.err logs/

###################################
# I.3 # calculate gene family trees
###################################

# !!!!!! must check first the species to include : Liberibacter? -> make MLSA phylogeny first
# -> Liberibacter sp. strains (LIBAP and LIBSC) will not be included in further analyses
# because their high divergence (one order of magnitude larger) will lead to systematic errors

cd $agrodata
# exclude Liberibacter sp. from codon alignments
python $scripts/subsample_clustalaln.py $PWD/codonfam_aln_clustal $PWD/$subdbname\_codon_aln_clustal LIBAP LIBSC
# verify their absence
grep 'LIBAP\|LIBSC' $subdbname\_codon_aln_clustal/*
$scripts/lsfullpath.py $PWD/$subdbname\_codon_aln_clustal > $subdbname\_codon_aln_clustal_list
# 10787 files remain in '$subdbname\_codonfam_aln/'
# convert back to fasta
mkdir $subdbname\_codon_aln_fasta
# use of bioscripts.convert Python package, available at https://pypi.python.org/pypi/bioscripts.convert/0.4
convalign -i clustal -e fasta fasta $subdbname\_codon_aln_clustal/*

## compute trees of gene families on nucleic acid alignements
# 1) with PhyML if further use of Prunier.2.0 (uses GTR+G8+I model)
cd $agrodata
mkdir phyml_trees
mkdir nucgb_alns
mkdir logs/phyml_tree_logs
# [TO ADAPT] computation in parallel - here using a PBS/Torque queeu submission system
# parallel jobs are cast with predefined walltime and as many as possible tasks are executed by each job:
# a stack of tasks is distributed to parallel workers; involve call to task manager utils/getnextbig.
$scripts/lsfullpath.py $PWD/$subdbname\_codon_aln_fasta > $subdbname\_codon_aln_fasta_list
$scripts/phyml.qsub $agrodata $subdbname\_codon_aln_fasta_list $PWD/nucgb_alns $PWD/phyml_trees 1 500 $PWD/logs/phyml_tree_logs q1day $worker fasta


# 2) with RAxML if further use of Prunier.2.1 (uses GTR+G4+I model)
cd $agrodata
mkdir raxml_trees
mkdir logs/raxml_tree_logs
# [TO ADAPT] computation in parallel - here using a PBS/Torque queeu submission system
# parallel jobs are cast with predefined walltime and as many as possible tasks are executed by each job:
# a stack of tasks is distributed to parallel workers; involve call to task manager utils/getnextbig.
scripts/lsfullpath.py $PWD/$subdbname\_codonfam_aln > codonfam_aln_list
scripts/raxml.qsub $PWD/$subdbname\_codonfam_aln_list $PWD/raxml_trees 1 1000 $PWD/logs/raxml_tree_logs q1day


###################################
# I.4 # calculate species phylogeny
###################################

# selects faimilies that are in exactly one copy in each genome of the dataset (that are mostly devoid of duplications/transfers/losses)
cd $agrodata
mkdir reftree
mkdir codongb_alns_fasta

$scripts/gblocks.csh $PWD/$subdbname\_codon_aln_fasta_list $PWD/codongb_alns_fasta codon
python $scripts/get_unicopy_fams.py $PWD/codongb_alns_fasta universal_unicopy_alns_list --universal=$agrodata/code5/$dbname.code --exclude='LIBAP,LIBSC'
# 455 families

## 1) core genome concatenated alignement, tree with RAxML - computation of bootstraps in parallel
mkdir species_phylogeny_raxml
mkdir species_phylogeny_raxml/universal
cp universal_unicopy_alns_list species_phylogeny_raxml/universal/
python $scripts/concat.py species_phylogeny_raxml/universal/universal_unicopy_alns_list species_phylogeny_raxml/universal/universal_unicopy_alns.concat
# ! concatenated aln is in phylip format
cd $agrodata
mkdir logs/raxml_species_phylogeny
mkdir logs/raxml_species_phylogeny/universal
# do RAxML with CAT approximation (50 categories of rate heterogeneoty, no invariant sites)
echo "-m GTRCAT -c 50 -f a" > raxmloptfile
# launch several replicates in parallel with different random seed for rapid bootstrap analysis
scripts/raxml_simple.qsub $PWD/species_phylogeny_raxml/universal/universal_unicopy_alns.concat $PWD/species_phylogeny_raxml/universal/ $PWD/logs/raxml_species_phylogeny/universal 20 10 q1day phylip $PWD/raxmloptfile
# synthesis of parrallel results
cd species_phylogeny_raxml/universal/
mkdir synthesis
cat RAxML_bootstrap.universal_unicopy_alns.* > synthesis/RAxML_bootstrap.universal_unicopy_alns.concat
cat RAxML_bestTree.universal_unicopy_alns.* > synthesis/RAxML_bestTree.universal_unicopy_alns.concat
# check all best ML trees are the same... OK
raxmlHPC -t RAxML_bestTree.universal_unicopy_alns.1 -z synthesis/RAxML_bootstrap.universal_unicopy_alns.concat -f b -m GTRCATI -n synthesis/RAxML_bestTree.universal_unicopy_alns.bs
python $scripts/addbranchlengths.py synthesis/RAxML_bestTree.universal_unicopy_alns.bs RAxML_bestTree.universal_unicopy_alns.1 synthesis/RAxML_bestTree.universal_unicopy_alns.len.bs


# 2) core ribosomal protein concatenated alignement, tree with RAxML - computation of bootstraps in parallel
# get ribosomal protein families that are unicopy and found in at least 45 species (over 47)
mkdir species_phylogeny_raxml/ribosomal45/
python $scripts/get_universal_unicopy_fams_query.py $PWD/codongb_alns_fasta $agrodata/code5/$dbname.code species_phylogeny_raxml/ribosomal45/ribosomal45_unicopy_alns_list 45 LIBAP LIBSC << EOF
q
XXXXX
pk(k=@ribosomal protein@) and kd(nk=gene_family)
EOF
python $scripts/concat.py species_phylogeny_raxml/ribosomal45/ribosomal45_unicopy_alns_list species_phylogeny_raxml/ribosomal45/ribosomal45_unicopy_alns.concat
# ! concatenated aln is in phylip format
mkdir logs/raxml_species_phylogeny/ribosomal45
# do RAxML with CAT approximation (50 categories of reta heterogeneoty, no invariant sites)
echo "-m GTRCAT -c 50 -f a" > raxmloptfile
# launch several replicates in parallel with different random seed for rapid bootstrap analysis
scripts/raxml_simple.qsub $PWD/species_phylogeny_raxml/ribosomal45/ribosomal45_unicopy_alns.concat $PWD/species_phylogeny_raxml/ribosomal45/ $PWD/logs/raxml_species_phylogeny/ribosomal45 20 50 q1day phylip $PWD/raxmloptfile
# synthesis of parrallel results
cd species_phylogeny_raxml/ribosomal45/
mkdir synthesis
cat RAxML_bootstrap.ribosomal45_unicopy_alns.* > synthesis/RAxML_bootstrap.ribosomal45_unicopy_alns.concat
cat RAxML_bestTree.ribosomal45_unicopy_alns.* > synthesis/RAxML_bestTree.ribosomal45_unicopy_alns.concat
# check all best ML trees are the same... OK
raxmlHPC -t RAxML_bestTree.ribosomal45_unicopy_alns.1 -z synthesis/RAxML_bootstrap.ribosomal45_unicopy_alns.concat -f b -m GTRCATI -n synthesis/RAxML_bestTree.ribosomal45_unicopy_alns.bs
python $scripts/addbranchlengths.py synthesis/RAxML_bestTree.ribosomal45_unicopy_alns.bs RAxML_bestTree.ribosomal45_unicopy_alns.1 synthesis/RAxML_bestTree.ribosomal45_unicopy_alns.len.bs


# 3) jacknife trees from sampling of core genome
mkdir jk_fams
mkdir jk_alns
mkdir jk_trees
mkdir logs/jk_tree_logs
python $scripts/mk_jk_fam.py $PWD/universal_unicopy_alns_list 500 25 $PWD/jk_fams
scripts/concat.sh jk_fams/ jk_alns/
# ! concatenated alns are in phylip format and already treated with GBlocks
scripts/lsfullpath.py $PWD/jk_alns > jk_aln_list
scripts/phyml_simple.qsub $PWD jk_aln_list NOGB $PWD/jk_trees $PWD/logs/jk_tree_logs 1 q1day phylip
# make a MRE consensus tree
mkdir species_phylogeny_jacknife
cat jk_trees/* > species_phylogeny_jacknife/jk_trees.concat
# root with Parvibaculum = 8th label to appear in first tree
phylip consense << EOF
species_phylogeny_jacknife/jk_trees.concat
O
8
Y
EOF
mv species_phylogeny_jacknife/outtree species_phylogeny_jacknife/jk_trees.mreconsensus
python $scripts/addbranchlengths.py species_phylogeny_jacknife/jk_trees.mreconsensus species_phylogeny_jacknife/jk_trees.concat species_phylogeny_jacknife/jk_trees.mreconsensus.len
raxmlHPC -t species_phylogeny_jacknife/jk_trees.mreconsensus.len -z species_phylogeny_jacknife/jk_trees.concat -f b -m GTRCATI -n species_phylogeny_jacknife/jk_trees.mreconsensus.len.bs

mkdir reftree
# open the several species trees with SEAVIEW and save them as rooted with Parvibaculum/PARL1 in directory 'reftree/'
python $scripts/make_reftree.py $agrodata/reftree/jk_trees.mreconsensus.len.bs $agrodata/code5/dic_spe_code $agrodata/code5/dic_taxid_code $agrodata/reftree/$subdbname

# choose the consensus of jackknife-sampled core gene concatenate (method 3) as reference species tree from this point
setenv reftree $agrodata/reftree/jk_trees.mreconsensus.len.bs

# all done for speciecies phylogeny at 18/05/12

#############################################
# II # Species tree/gene tree reconciliations
#############################################

# this part of the procedure is specific to the making of Agrogenom database


########################################################
# II.0 # Initiate storage of information in SQL database
########################################################


## Use of a PostGreSQL database ; schema inspired from Phylariane scheme, adding block information and differentiating DTLSU events in separate tables
# [TO ADAPT] change server and database/user names
# create database (server side) ; require system admin rights and log in as 'postgres' usercluster
ssh pbil-sgbd
psql -h localhost -U postgres -c 'CREATE DATABASE agrogenom ; CREATE USER lassalle' # additionally, need to alloacte the appropriate admin rights to the user
logout
# loads the 'ltree' module
psql -h localhost -U postgres -d $yourdbname -f /usr/share/postgresql/8.3/contrib/ltree.sql
# builds the database from modified Phylariane scheme
cd $agrodata
psql -h $yourservername -U $yourusername -d $yourdbname -f database/agrogenom_schema.sql
bunzip2 database/agrogenom_gene_table.csv.bz2
psql -h $yourservername -U $yourusername -d $yourdbname
\copy genome.gene from '/pandata/lassalle/agrogenom/database/agrogenom_gene_table.csv' with delimiter as '|'
UPDATE genome.gene SET tax_id = s.tax_id FROM (SELECT tax_id, number FROM phylogeny.species_node ) as s WHERE split_part(gene.chromosome, '_', 1)=s.number ;

# loads the `species_tree` and `species_node` table from reference tree
python $scripts/make_reftree.py $agrodata/reftree/jk_trees.mreconsensus.len.bs $agrodata/code5/dic_spe_code $agrodata/code5/dic_taxid_code $agrodata/reftree/$subdbname
psql -h $yourservername -U $yourusername -d $yourdbname -f $agrodata/reftree/$subdbname.species_tree_node_tables.sql

## map back transfer events on full annotated gene trees and refine duplication annotations on a forest of subtrees from gene trees pruned at transfer nodes
scripts/lsfullpath.py $currentrec/prunier_out > $currentrec/prunier_out_list
# 38165 reconciled unicopy subtrees

#############################
# II.1 # (re-)root gene trees
#############################

## use TPMS to root the gene trees using the combo criterion, so that, knowing the species phylogeny, both the species multiplicity and the taxonomic depth of subtrees are minimized.
# A root minimizing these criteria favours reconciliation scenarios more parsimonious in ancient gain (duplication and transfer) events (which implies more losses).
# (19/05/12)
cd $agrodata
mkdir tpms_db
# prepare database
$TPMS/tpms_mkdb -extract-seqnames -families-trees=$agrodata/phyml_trees -output=$agrodata/tpms_db/seqNames2spe
python $scripts/fill_tpms_seqlist.py tpms_db/seqNames2spe tpms_db/seqNames2species
# have to create a modified species tree file bearing the label 'ROOT:0' at the root ; save as $agrodata/reftree/$subdbname.reftree.tpms
# create database
$TPMS/tpms_mkdb -families-trees=$agrodata/phyml_trees -seq-to-species=$agrodata/tpms_db/seqNames2species -sp-tree=$agrodata/reftree/$subdbname.reftree.tpms -output=$agrodata/tpms_db/$subdbname.unrooted
# root gene trees
$TPMS/tpms_computations -collection=$agrodata/tpms_db/$subdbname.unrooted -output-dir=$agrodata/tpms_db/ -threads=8 << EOF
RC
S
EOF
mv tpms_db/Collection tpms_db/$subdbname.rooted
# split the result into unique tree files
mkdir rooted_trees
python $scripts/split_tpms_db.py $agrodata/tpms_db/$subdbname.rooted $agrodata/rooted_trees

# TPMS roots poorly trees in absence of duplication and presence of transfers

# Use alternatie method: find tranfers in unicopy gene trees to then re-root them consistently to the parcionious transfer scenario
# identify unicopy gene families ; use alignments that were already filtered by GBlocks to positions with phylogenetic signal
python $scripts/get_unicopy_fams.py $agrodata/nucgb_alns $agrodata/unicopy_nucgb_aln_list

# prunier detection of transfers
mkdir $agrodata/unicopy_prunier_out
# [TO ADAPT] computation in parallel - here using a PBS/Torque queue submission system
# parallel jobs are cast with predefined walltime and as many as possible tasks are executed by each job:
# a stack of tasks is distributed to parallel workers; involve call to task manager utils/getnextbig.
$scripts/prunier.qsub $agrodata/unicopy_nucgb_aln_list $agrodata/reftree/$subdbname.reftree.prunier $agrodata/unicopy_prunier_out $agrodata/logs/prunier q1day 1 400 dna 0.9 $agrodata/phyml_trees

# rooting of trees as by Prunier
mkdir $agrodata/unicopy_rooted_trees
python $scripts/root_as_prunier.py $agrodata/unicopy_gene_tree_list $reftree $agrodata/unicopy_prunier_out $agrodata/unicopy_rooted_trees
cp -f $agrodata/unicopy_rooted_trees/* $agrodata/rooted_trees

################################################################
# II.2 # find putative duplications and extract unicopy subtrees
################################################################


# find duplications in gene trees - computation in parallel (20/06/12)
setenv currentrec $agrodata/duplications/maxlossrate02
mkdir -p $currentrec
mkdir $agrodata/logs/duplications
scripts/lsfullpath.py $agrodata/rooted_trees > rooted_tree_list
# launch with parameter maxlossrate = 0.2 (hard-coded in finddup.csh)
# [TO ADAPT] computation in parallel - here using a PBS/Torque queue submission system
# parallel jobs are cast with predefined walltime and as many as possible tasks are executed by each job:
# a stack of tasks is distributed to parallel workers; involve call to task manager utils/getnextbig.
scripts/finddup.qsub $agrodata/rooted_tree_list $agrodata/nucgb_alns $currentrec $reftree $agrodata/logs/duplications q1day 1 300

# check coherrence of result files
$scripts/lsfullpath.py $currentrec/modified_trees > $currentrec/modified_tree_list
$scripts/lsfullpath.py $currentrec/normed_trees > $currentrec/normed_tree_list
$scripts/lsfullpath.py $currentrec/subtrees > $currentrec/subtree_list
$scripts/lsfullpath.py $currentrec/subfams > $currentrec/subfam_list
$scripts/lsfullpath.py $currentrec/treeshelves > $currentrec/shelf_list
cut $currentrec/subfam_list -d'.' -f1 | awk -F'/' '{print $NF}' | sort -u > $currentrec/fams_in_subfam

# verify completion of jobs
$scripts/whichfam.py $currentrec/treeshelves > $currentrec/famsin_treeshelves
$scripts/whichfam.py $currentrec/modified_trees > $currentrec/famsin_modified_trees
diff $currentrec/famsin_modified_trees $currentrec/famsin_treeshelves | grep "<" | cut -d' ' -f 2 > $currentrec/fams_toredo_finddups
rm -f $currentrec/trees_toredo_finddups && touch $currentrec/trees_toredo_finddups
foreach fam (`cat $currentrec/fams_toredo_finddups`)
echo "$agrodata/rooted_trees/$fam.phb" >> $currentrec/trees_toredo_finddups
end
# following families did not complete:
# 49RHIZOB_1203 : family of EF-Tu has a pair of copy in each genome, but those are always closer from their paralog than from any ortholog, certainly due to continuous homogeneization by intra-genomic gene conversion >> a nightmare for the duplication-minimization algorithm, that exceeds recursion limits. Are missing in a lot of genomes (due to frequent rejection of duplicated sequences in recent assemblies)
# 49RHIZOB_84_1 : just a big family (888 members) of ABC transporter inner membrane proteins ; but same prolem as above (recursion limit exceeded).
$scripts/finddup_simple.qsub $currentrec/trees_toredo_finddups $agrodata/nucgb_alns $currentrec $reftree $agrodata/logs/duplications q1week
# in both cases topology amelioration run indefinetely + cannot find combinations of unicopy subtrees in an appropriate time

# finally (12/07/12):
# 10774 trees and full families : OK
# 28343 sub-families of orthologs
# 38186 non-redundant unicopy subtrees

########################################################################
# II.3 # detect HGT based on statistically suported topological conflict
########################################################################

# fetch unicopy sub-family alignments from previous family alignments
mkdir $currentrec/subalns_fasta
python $scripts/extract_subaln.py $currentrec/subtree2leaves $agrodata/nucgb_alns_fasta codon-gb $currentrec/subalns_fasta

## perform horizontal transfer detection by Prunier on unicopy subtrees (06/07/12)
mkdir $currentrec/prunier_out
mkdir logs/prunier
$scripts/lsfullpath.py $currentrec/subalns_fasta > $currentrec/subalns_fasta_list
$scripts/prunier.qsub $currentrec/subalns_fasta_list $agrodata/reftree/$subdbname.reftree.prunier $currentrec/prunier_out $agrodata/logs/prunier q1day 1 400 dna 0.9 $currentrec/subtrees

###########################################################################
# II.4 # Completion and integration + db upload of reconciliation scenarios
###########################################################################

# within the same script rec_to_db.py, several key steps:
# - Prunier reconciliation of many subtrees from a same gene family tree are interpreted (inference of HGT donor and recipients)
#   and integrated into one single gene tree reconciliation;
# - for the need of the db upload procedure (see below), reconciliation scenarios need to be complete, i.e. that every node of the gene tree has an avent associated to it.
#	For this reason, these reconciliations are completed (i.e. remaining unexplained topological conflict between the gene tree and the species tree)
#   by infering additional HGT based on minimization of subtrees taxonomic depth (using taxonomic incongruence scoring algorithm from TPMS (Bigot et al. 2013) as re-implemented in this script).
#   Note however that this inference will be re-made at step II.7, based on the final set of inferred HGT events. This step of the procedure thus leads to infering
#	events that have no final value, and for this reason, it is not described in the pipeline summary
# - load detailed (e.g. source of HGT inference, mapping of unicopy subtrees) and global (list of all events per node) reconciliation 
#   information in SQL database.

# (22/01/13)
# [TO ADAPT] computation in parallel - here using a PBS/Torque queue submission system
# parallel jobs are cast with predefined walltime and as many as possible tasks are executed by each job:
# a stack of tasks is distributed to parallel workers; involve dynamic query of the database.
cd $currentrec
scripts/lsfullpath.py $currentrec/modified_trees > $currentrec/modified_tree_list
# jobs query a task table in the SQL db to fetch a new task; conccurent access of jobs to task table is avoided with explicit locks to prevent creating duplicate
psql -h $yourservername -U $yourusername -d $yourdbname << EOF
CREATE TABLE IF NOT EXISTS reconciliation_to_db (gene_tree_file_path varchar(500) PRIMARY KEY, current_status int DEFAULT 0, date timestamp with time zone DEFAULT now(), job_id DEFAULT 0);
\copy reconciliation_to_db (gene_tree_file_path) from $currentrec/modified_tree_list
UPDATE reconciliation_to_db SET current_status=0, date=now(), job_id=0;
EOF
$scripts/rec_to_db.qsub public.reconciliation_to_db $currentrec/subtrees $currentrec/subtree2leaves $currentrec/prunier_out $agrodata/reftree/47RHIZOB.reftree.phb X $currentrec/rec_to_db_220113 0.9 1 40 q1day $agrodata/logs/rec_to_db
# database is not updated in live (too inefficient due to locks on tables during concurrent access)
# a data dump is generated that will be imported to the db at the next step.


# import of data (with interactive validation of commits)
# copy dumps and delete duplicates created by simultaneous jobs working on the same families (legacy requirement from before instauring the table lock system, cf. above)
python $scripts/import_sql_dumps.py -c -i -d --job.dumps.dir=$currentrec/rec_to_db_220113 --concatenated.dump.dir=$currentrec/rec_to_db_220113/concat_sql_dump

# check redundancy of duplicate rows for the same unicopy subtree create by different jobs working on the same task (again, should not happen)
SELECT subtree_name, unicopy_subtree_id FROM phylogeny.unicopy_subtree INNER JOIN (SELECT subtree_name FROM phylogeny.unicopy_subtree GROUP BY subtree_name HAVING count(unicopy_subtree_id)>1) AS duplicates USING (subtree_name);
SELECT * FROM phylogeny.reconciled_gene_tree_node INNER JOIN (SELECT rec_gi_id, node_id FROM phylogeny.reconciled_gene_tree_node WHERE count(nr_node_id)>1 GROUP BY rec_gi_id, node_id) AS duplicates USING (rec_gi_id, node_id) ;

# Because several gene subtree may overlap in a same gene tree (in presence of in-paralogs, different sets of co-orthologs were extracted),
# and were independently reconciled by Prunier, it is possible to infer different events at the same gene tree node
# all inferred events are recorded in db, and are given different unique identifiers (node_id) ;
# the correspoding gene tree nodes are given a non-redundant identifier (nr_node_id) that is shared by potentially several events
# check that every non-redundant node has at least one event affected to it, and that when there are several events, that they are of one single type:
psql -h $yourservername -U $yourusername -d $yourdbname
select count(*) from phylogeny.reconciled_gene_tree_node;
select count(*) from (select distinct rec_gi_id, rec_gi_start_node from phylogeny.event ) as foo ;
select count(*) from (select distinct rec_gi_id, rec_gi_start_node, event_type from phylogeny.event ) as foo ;
# each of the three return values are 467528 : OK.

# make a backup the database schema and data at this step
pg_dump -U $yourusername -h $yourservername -f $agrodata/database/agrogenom_070113.dump -F c agrogenom


####################################################
# II.5 # Build genomic blocks of evolutionary events
####################################################

# This script getblockevents.py finds block events and create the corresponding objects in the database.
# It first explore the contemporary genomes (replicon by replicon), recognizing events on their lineage that are similar between neighbouring genes and aggregate them into "leaf block events";
# In asecond step, leaf block events from different genomes referring to the same evolutionary (duplication, transfer) events are aggregated into an "ancestral block event".
# 
# Ancestral blocks are thus the units of actual evotutionary events that gather homologous single gene events set in the past
# (i.e. located at internal gene tree nodes) - or rather, gathering independent observations in different gene families of the same evolutionary event.
# Leaf blocks are the descendent of the ancestral blocks, i.e. their realization in contemporary genomes.
#
# [TO ADAPT] computation in parallel - here using a PBS/Torque queue submission system
# parallel jobs are cast with predefined walltime and as many as possible tasks are executed by each job:
# a stack of tasks is distributed to parallel workers; involve dynamic query of the database.
# jobs query a task table in the SQL db to fetch a new task; conccurent access of jobs to task table is avoided with explicit locks to prevent creating duplicate.
psql -h $yourservername -U $yourusername -d $yourdbname <<EOF
CREATE TABLE public.buildingblocks AS (SELECT chromosome, count(*) as nb_genes, 0 as job_id, 0 as current_status, now() as date FROM genome.gene where chromosome!='LIBAP_1' AND chromosome!='LIBSC_1' GROUP BY chromosome);
EOF

# leaf blocks reconstruction
$scripts/getblockevents.qsub $currentrec/rec_to_db_220113/reconciled_tree_pickles $reftree $currentrec/blockevents_240113 $agrodata/logs/getblockevents q1day 80 public.buildingblocks l speciation,duplication:old.N1.N2.N3
# ancestral blocks reconstruction 
$scripts/getblockevents.qsub $currentrec/rec_to_db_220113/reconciled_tree_pickles $reftree $currentrec/blockevents_240113 $agrodata/logs/getblockevents q1day 1 public.buildingblocks a speciation,duplication:old.N1.N2.N3
setenv blocks $currentrec/blockevents_240113
setenv ancblocks $blocks/getblockevent_out.a.81
cd $blocks/getblockevent_out.a.81/
# modify table formats (enum() types must not be quoted)
foreach file (`ls blocks_db_dump/*`)
sed -e "s/'//g" $file > $file.2
rm $file ; mv $file.2 $file
end
# load in db
psql -h $yourservername -U $yourusername -d $yourdbname < $ancblocks/blocks_db_dump/allblocktables_sqldump.sql
# make a backup the database schema and data at this step
pg_dump -U $yourusername -h $yourservername -f $agrodata/database/agrogenom_dump_210213.tar -F t agrogenom


##################################
# II.6 # Reconciliation refinement
##################################

# When multiple transfer events were proposed by independent Prunier inferences on overlapping unicopy subtrees ( = replicates)
# for a same non-redundant gene tree node, choose the best reconciliation using the the block event information,
# i.e. favouring the evnts that are found in the largest blocks of genes. When this criterion is not discriminating
# (most events involve single gens), the frequency of equivalent transfer events among the independent replicates.
# The choice of the optimal event may create the cohabitation of events from different replicates, which does not guaratee
# the consistency of the scenario; the potential conflict is resolved by a posteriori integration and completion of the reconciliation:
# remaining unexplained topological conflict between the gene tree and the species tree is resolved by infering additional HGT
# based on minimization of subtrees taxonomic depth (using XD algorithm from TPMS (Bigot et al. 2013) as re-implemented in this script)
# (call to same functions as in step II.5).

qsub -q q1day -N refine_transfer -e $agrodata/logs/refine_transfer.stderr -o $agrodata/logs/refine_transfer.stdout << EOF
python $scripts/refine_transfer_annotation.py -d phylariane -p $dbpwd $currentrec/modified_trees $agrodata/reftree/47RHIZOB.reftree.phb $currentrec/refined_trees_260313
EOF

# update database
psql -h $yourservername -U $yourusername -d $yourdbname
\copy phylogeny.event from '/pandata/lassalle/agrogenom/duplications/maxlossrate02/refined_trees_260313/event_dump.tab'
\copy phylogeny.event_possible_species_node from '/pandata/lassalle/agrogenom/duplications/maxlossrate02/refined_trees_260313/event_possible_species_node_dump.tab'
INSERT INTO phylogeny.reconciliation_collection VALUES ('2013-03-26', '/pandata/lassalle/agrogenom/duplications/maxlossrate02/refined_trees_260313');
\copy phylogeny.represented_event FROM '/pandata/lassalle/agrogenom/duplications/maxlossrate02/refined_trees_260313/represented_event_dump.tab'

# make a backup the database schema and data at this step
pg_dump -U $yourusername -h $yourservername -f $agrodata/database/agrogenom_dump_260313.tar -F t agrogenom


#########################################
# II.7 # Ancestral content reconstruction
#########################################

# separation of families into orthologous subfamilies and reconstruction of ancestral states with Wagner pasimony (call to COUNT program (Csuros 2008))
# this leads to the inference of the last set of additional HGT events, based on heterogeneous representation of species within othologuous gene subfamilies.

mkdir $currentrec/ancestral_content_230413
python $scripts/ancestral_content.py -d phylariane -p $dbpwd --infer='Count.AsymmetricWagner' -o '-gain 1.5' --tmp-dir=$currentrec/ancestral_content_230413/tmp15 --output-subfamily-trees=true --output-pickled-family-histories=true --output-genome-synthesis=true $agrodata/nucfam_fasta_list $currentrec/refined_trees_260313/reconciled_tree_pickles $reftree $currentrec/ancestral_content_230413 >& $currentrec/ancestral_content_230413.AsymetricWagner.g_15.oe &
setenv anccont $currentrec/ancestral_content_230413/Count.AsymmetricWagner.gain_15/

# new trees (with additional transfers) have been generated. complete the annotated tree collection with previously calculated trees
# complete with the unchanged ones
mkdir $anccont/reconciliation_collection
cp -r $anccont/modified_trees/* $anccont/reconciliation_collection/
cp -np $currentrec/refined_trees_260313/normed_phyloxml_trees/*.xml $anccont/reconciliation_collection/normed_phyloxml_trees/
cp -np $currentrec/refined_trees_260313/reconciled_tree_pickles/*.pickle $anccont/reconciliation_collection/reconciled_tree_pickles/
python $scripts/dump_reconciliation_collection_tables.py -c $anccont/reconciliation_collection/ -r $reftree -i '2013-04-23' -t -a

# store the collection into database
psql -h $yourservername -U $yourusername -d $yourdbname
\copy phylogeny.event from '/pandata/lassalle/agrogenom/duplications/maxlossrate02/ancestral_content_230413/Count.AsymmetricWagner.gain_15/reconciliation_collection/event_dump.tab'
\copy phylogeny.event_possible_species_node from '/pandata/lassalle/agrogenom/duplications/maxlossrate02/ancestral_content_230413/Count.AsymmetricWagner.gain_15/reconciliation_collection/event_possible_species_node_dump.tab'
\copy phylogeny.reconciliation_collection FROM '/pandata/lassalle/agrogenom/duplications/maxlossrate02/ancestral_content_230413/Count.AsymmetricWagner.gain_15/reconciliation_collection/reconciliation_collection_dump.tab'
\copy phylogeny.represented_event FROM '/pandata/lassalle/agrogenom/duplications/maxlossrate02/ancestral_content_230413/Count.AsymmetricWagner.gain_15/reconciliation_collection/represented_event_dump.tab'
-- superfluous reconciliation collections
delete from phylogeny.reconciliation_collection where rec_col_id NOT IN ('2013-03-26', '2013-04-23');
-- to delete the burden of redundent and unused events
create temporary table repr_evt as (select * from phylogeny.event left join phylogeny.represented_event using (event_id));
-- superfluous speciations and duplications
delete from phylogeny.event where event_id in (select distinct(event_id) from repr_evt WHERE event_type IN ('S', 'D') AND rec_col_id IS NULL);
-- superfluous transfers
delete from phylogeny.event where event_id in (select distinct(event_id) from repr_evt LEFT JOIN analysis.prunier_inference USING (event_id) where event_type='T' AND prinf_id IS NULL AND rec_col_id IS NULL);


#################################
# II.8 # Finalizing and synthesis
#################################

## regenrates blocks with new events documented
# clean database
psql -h $yourservername -U $yourusername -d $yourdbname
CREATE TABLE public.buildingblocks AS (SELECT chromosome, count(*) as nb_genes, 0 as job_id, 0 as current_status, now() as date FROM genome.gene where chromosome!='LIBAP_1' AND chromosome!='LIBSC_1' GROUP BY chromosome);
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
# make a backup the database schema and data at this step
pg_dump -U $yourusername -f agrogenomdb_230413.dump.tar -F t agrogenom

mkdir $currentrec/blockevents_240413

# sequential run (could be parallel, but not as many events to consider now)
# leaf blocks reconstruction
python $scripts/getblockevents.py $anccont/reconciliation_collection/reconciled_tree_pickles $reftree $currentrec/blockevents_240413 99 -c public.buildingblocks -e speciation,duplication:old.N1.N2.N3,gain:N1.N2.N3 -d phylariane -p $dbpwd >& $currentrec/blockevents_240413.oe &
python $scripts/getblockevents.py $anccont/reconciliation_collection/reconciled_tree_pickles $reftree $currentrec/blockevents_240413 98 -c public.buildingblocks -e speciation,duplication:old.N1.N2.N3,gain:N1.N2.N3 -d phylariane -p $dbpwd >& $currentrec/blockevents_240413.oe2 &
# ancestral blocks reconstruction !!! high memory use
qsub -q q1tera -l mem=60g -N getblockevents.a.1 -e $agrodata/logs/getblockevents.a.1.stderr -o $agrodata/logs/getblockevents.a.1.stdout << EOF
python $scripts/getblockevents.py $anccont/reconciliation_collection/reconciled_tree_pickles $reftree $currentrec/blockevents_240413 1 -a -v -c public.buildingblocks -e speciation,duplication:old.N1.N2.N3,gain:N1.N2.N3 -d phylariane -p $dbpwd
EOF
setenv blocks $currentrec/blockevents_240413
setenv ancblocks $blocks/getblockevent_out.a.1
psql -h $yourservername -U $yourusername -d $yourdbname < $ancblocks/blocks_db_dump/allblocktables_sqldump.sql

# make a backup the database schema and data at this step
pg_dump -U $yourusername -f agrogenomdb_240413.dump.tar -F t agrogenom

# generate synthesis files: genome history and clade-specific gene tables
python $scripts/ancestral_content.py -d phylariane -p $dbpwd --input-family-history-pickles=$anccont/family_history_pickles --input-subfamily-trees=$anccont/ortho_subtrees --output-all --output-subfamily-trees=false --output-pickled-family-histories=false $agrodata/nucfam_fasta_list $currentrec/refined_trees_260313/reconciled_tree_pickles $reftree $anccont/synthesis >& $currentrec/ancestral_content_230413.synthesis.oe


######################################
# III # Ananalysis of genome histories
######################################

#########################################
# III.1 # Load gene sets into database
#########################################

# load subfamily informations into DB
cd $anccont/synthesis/
psql -h $yourservername -U $yourusername -d $yourdbname
\copy genome.gene2subfam from 'ortho_subfam_leaves.tab'
\copy phylogeny.event2subfam from 'ortho_subfam_events.tab'

# dump phylogenetic profile and clade-specific gene sets to DB
# create dump data files
python $scripts/spegenes_to_db.py $anccont/synthesis/annottables $reftree $anccont/synthesis/specific_gene_dump.tab
python $scripts/phyloprofile_to_db.py $anccont/synthesis/phylogenetic_profiles/subfam_node_counts.mat $anccont/synthesis/phylogenetic_profile_dump.tab

# create tables and load data
psql -h $yourservername -U $yourusername -d $yourdbname
CREATE TABLE genome.phylogenetic_profile (
	subfam_id varchar(50) NOT NULL,
	number varchar(50) NOT NULL,
	present boolean,
	count smallint,
	PRIMARY KEY (subfam_id, number)
);
CREATE INDEX phylogetic_profile_subfam_id_key ON genome.phylogenetic_profile (subfam_id);
CREATE INDEX phylogetic_profile_number_key ON genome.phylogenetic_profile (number);
\copy genome.phylogenetic_profile from 'phylogenetic_profile_dump.tab'

CREATE TABLE genome.specific_gene (
	subfam_id varchar(50) NOT NULL,
	number varchar(50) NOT NULL,
	specificity varchar(20) NOT NULL,
	relaxed smallint,
	PRIMARY KEY (subfam_id, number, specificity)
);
CREATE INDEX specific_genes_subfam_id_key ON genome.specific_gene (subfam_id);
CREATE INDEX specific_genes_number_key ON genome.specific_gene (number);
CREATE INDEX specific_genes_specificity_key ON genome.specific_gene (specificity);
\copy genome.specific_gene from 'specific_gene_dump.tab' 


###############################################
# III.2 # Functional homogeneity of gene blocks
###############################################

# [TO DO YOURSELF] requires two table:
# - genome.go: contain the reference Gene Ontology term id and their functional annotation; was obtained from Gene Ontology consortium website
# - genome.gene_go: a mapping of genes to GO terms; was obtained from running InterProScan on all CDSs and then using UniProt IP2GO mapping.
# annotate subfamilies with GO terms shared by at least 50% of its members
psql -h $yourservername -U $yourusername -d $yourdbname
create table genome.subfam_go as (
	select subfam_id, go.go_id, count(gene_id) 
		from genome.gene inner join genome.gene_go using (gene_id) 
		inner join genome.go using (go_id) 
		inner join genome.gene2subfam using (hogenom_gene_id) 
	group by subfam_id, go.go_id
);


mkdir $agrodata/logs/funsim
# [TO ADAPT] computation in parallel - here using a PBS/Torque queue submission system
psql -h $yourservername -U $yourusername -d $yourdbname -A -t -c "select distinct(number) from phylogeny.species_node where left(number, 3) IN ('AGR', 'ATU');" > $agrodata/GeneOntology/Agrobacterium.taxids
foreach spe (`cat $agrodata/GeneOntology/Agrobacterium.taxids`)
qsub -N funsim.$spe -q q1day  -e $agrodata/logs/funsim/funsim.$spe.stderr -o $agrodata/logs/funsim/funsim.$spe.stdout << EOF
python $scripts/score_genegroup_funsim.py -t $spe -r $agrodata/GeneOntology/rand_gene_pair_sim -o $agrodata/GeneOntology/sggf_out -g $agrodata/GeneOntology/go_daily-termdb.obo-xml -w 2:10 -b -p $dbpwd
EOF
end
# summary of results, generate graphics
R --vanilla --silent < $scripts/score_genegroup_funsim.r
psql -h $yourservername -U $yourusername -d $yourdbname 
# store in database
\copy blocks.functional_homogeneity from '/pandata/lassalle/agrogenom/GeneOntology/sggf_out/funsim_dump.tab'
