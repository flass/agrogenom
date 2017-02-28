create temporary table distr_uncertainty as (
  select event_type, characteristic, nb_rec, count(event_id) 
  from ( 
    select event_id, event_type, characteristic, count(sp_node_id) as nb_rec 
      from phylogeny.event 
      inner join phylogeny.event_possible_species_node using (event_id) 
      inner join  phylogeny.represented_event using (event_id) 
    where rec_col_id='2013-08-14' 
    group by event_id, event_type, characteristic
  ) as a 
  group by event_type, characteristic, nb_rec 
  order by event_type, characteristic, nb_rec
);

\copy distr_uncertainty to 'distr_uncertainty.tab'

create temporary table diffd as (
 select * from (
  select anc_block_id, pe.event_id, nb_epsn, nb_bpsn, count(evc.event_id) as nb_ev from (
    select e.event_id, count(epsn.sp_node_id) as nb_epsn, e.anc_block_id
      from phylogeny.event as e 
      inner join phylogeny.event_possible_species_node as epsn using (event_id)
    where e.event_type='D'
    group by e.event_id, e.anc_block_id
  ) as pe inner join (
    select ab.anc_block_id, count(bpsn.sp_node_id) as nb_bpsn
      from blocks.ancestral_block as ab 
      inner join blocks.block_possible_species_node as bpsn on bpsn.block_id=ab.anc_block_id
    group by ab.anc_block_id
  ) as pb using (anc_block_id) 
  inner join (select event_id, anc_block_id from phylogeny.event) as evc using (anc_block_id)
  group by pe.event_id, nb_epsn, anc_block_id, nb_bpsn
 ) as gth where nb_ev>1
);

\copy diffd to 'compare_duplication_event_block_location_precision.tab'

create temporary table diffg as (
 select * from (
  select anc_block_id, pe.event_id, nb_epsn, nb_bpsn, count(evc.event_id) as nb_ev from (
    select e.event_id, count(epsn.sp_node_id) as nb_epsn, e.anc_block_id
      from phylogeny.event as e 
      inner join phylogeny.event_possible_species_node as epsn using (event_id)
    where e.event_type='G'
    group by e.event_id, e.anc_block_id
  ) as pe inner join (
    select ab.anc_block_id, count(bpsn.sp_node_id) as nb_bpsn
      from blocks.ancestral_block as ab 
      inner join blocks.block_possible_species_node as bpsn on bpsn.block_id=ab.anc_block_id
    group by ab.anc_block_id
  ) as pb using (anc_block_id) 
  inner join (select event_id, anc_block_id from phylogeny.event) as evc using (anc_block_id)
  group by pe.event_id, nb_epsn, anc_block_id, nb_bpsn
 ) as gth where nb_ev>1
);

\copy diffg to 'compare_gain_event_block_location_precision.tab'


create temporary table difft as (
 select * from (
  select anc_block_id, pe.event_id, nb_repsn, nb_depsn, nb_rbpsn, nb_dbpsn, count(evc.event_id) as nb_ev
  from (
    select e.event_id, nb_repsn, nb_depsn, e.anc_block_id
      from phylogeny.event as e 
      inner join (
        select event_id, count(sp_node_id) as nb_repsn from phylogeny.event_possible_species_node where characteristic='rec' group by event_id
      ) as repsn using (event_id)
      inner join (
        select event_id, count(sp_node_id) as nb_depsn from phylogeny.event_possible_species_node where characteristic='don' group by event_id
      ) as depsn using (event_id)
    where e.event_type='T'
  ) as pe inner join (
    select ab.anc_block_id, nb_rbpsn, nb_dbpsn
      from blocks.ancestral_block as ab 
      inner join (
        select block_id, count(sp_node_id) as nb_rbpsn from blocks.block_possible_species_node where block_type='ancestral' and characteristic='rec' group by block_id
      ) as rbpsn on rbpsn.block_id=ab.anc_block_id
      inner join (
        select block_id, count(sp_node_id) as nb_dbpsn from blocks.block_possible_species_node where block_type='ancestral' and characteristic='don' group by block_id
      ) as dbpsn on dbpsn.block_id=ab.anc_block_id
  ) as pb using (anc_block_id) 
  inner join (select event_id, anc_block_id from phylogeny.event) as evc using (anc_block_id)
  group by pe.event_id, nb_repsn, nb_depsn, anc_block_id, nb_rbpsn, nb_dbpsn
 ) as gth where nb_ev>1
);

\copy difft to 'compare_transfer_event_block_location_precision.tab'

-- distribution of sizes of transferred leaf block according to their (ancestral) recipient (in projection on current genomes)
create temporary table distlblreftree as (
  select number, currentgenome, lblen, count(distinct(leaf_block_id)) 
    from (
      select leaf_block_id, number, currentgenome, count(distinct(gene_id)) as lblen 
          from (
            select gene_id, split_part(chromosome, '_', 1) as currentgenome from genome.gene
          ) as g 
          inner join blocks.gene2block using (gene_id) 
          inner join blocks.leaf_block using (leaf_block_id) 
          inner join blocks.ancestral_block using (anc_block_id) 
          inner join blocks.block_possible_species_node on block_id=leaf_block_id 
          inner join phylogeny.species_node using (sp_node_id) 
          where characteristic='rec' 
          and event_type IN ('T', 'G') 
          and reference_node is true 
      group by leaf_block_id, number, currentgenome
    ) as lbl 
  group by number, currentgenome, lblen
);
\copy distlblreftree to 'distrib_size_leafblocks_by_rec_by_currgenome.tab'


-- distribution of of functional homogeneity of transferred leaf block according to their (ancestral) recipient (in projection on current genomes)
create temporary table distlbfhreftree as (
  select * from( 
    select leaf_block_id, number, currentgenome, count(distinct(gene_id)) as lblen , funsimP, funsimF, ((funsimP^2)+(funsimF^2))/2 as funsim
      from (
        select gene_id, split_part(chromosome, '_', 1) as currentgenome from genome.gene
      ) as g 
      inner join blocks.gene2block using (gene_id) 
      inner join blocks.leaf_block using (leaf_block_id) 
      inner join (
        select leaf_block_id, funsim as funsimP from blocks.functional_homogeneity
        where pairwise_metric='funSimMax'
        and group_metric='meanMaxSim'
        and aspect ='P'
      ) as fsp using (leaf_block_id)
      inner join (
        select leaf_block_id, funsim as funsimF from blocks.functional_homogeneity
        where pairwise_metric='funSimMax'
        and group_metric='meanMaxSim'
        and aspect ='F'
      ) as fsf using (leaf_block_id)
      inner join blocks.ancestral_block using (anc_block_id) 
      inner join blocks.block_possible_species_node on block_id=leaf_block_id 
      inner join phylogeny.species_node using (sp_node_id) 
      where characteristic='rec' 
      and event_type IN ('T', 'G') 
      and reference_node is true 
    group by leaf_block_id, number, currentgenome, funsimF, funsimP
  ) as fblb where funsim is not null
);
\copy distlbfhreftree to 'distrib_funcsim_leafblocks_by_rec_by_currgenome.tab'




-- number of events
select event_type, count(*) from phylogeny.event inner join phylogeny.represented_event using (event_id) where rec_col_id='2013-08-14' group by event.event_type;
 --~ event_type | count  
--~ ------------+--------
 --~ T          |  43233
 --~ G          |   5189
 --~ S          | 411766
 --~ D          |   7340

-- number of events eligible for block event reconstruction
select event_type, count(*) from phylogeny.event inner join phylogeny.represented_event using (event_id) where rec_col_id='2013-08-14'and anc_block_id is not null group by event.event_type;
 --~ event_type | count 
--~ ------------+-------
 --~ T          | 43131
 --~ G          |  2603
 --~ S          |   299
 --~ D          |  4406

-- number of block events
select event_type, count(*) from blocks.ancestral_block inner join (select distinct(anc_block_id) from phylogeny.event inner join phylogeny.represented_event using (event_id) where rec_col_id='2013-08-14') as r using (anc_block_id) group by event_type;
 --~ event_type | count 
--~ ------------+-------
 --~ D          |  3035
 --~ T          | 32139
 --~ G          |  1681

-- idem  specifically for  replacing transfers
-- number of events
select count(*) from phylogeny.event inner join phylogeny.represented_event using (event_id) inner join phylogeny.event2subfam using (event_id) left join (select rec_col_id as rec_col_id_old, event_id from phylogeny.represented_event where rec_col_id!='2013-08-14') as oldrec using (event_id) where rec_col_id='2013-08-14' and event.event_type='T' and origin='f' and rec_col_id_old is not null ;
 --~ count 
--~ -------
 --~ 9271
 
-- number of events eligible for block event reconstruction
select count(*) 
  from phylogeny.event 
  inner join phylogeny.represented_event using (event_id) 
  inner join phylogeny.event2subfam using (event_id) 
  left join (
    select rec_col_id as rec_col_id_old, event_id 
      from phylogeny.represented_event 
    where rec_col_id!='2013-08-14'
  ) as oldrec using (event_id) 
where rec_col_id='2013-08-14' 
and event.event_type='T' 
and origin='f' 
and rec_col_id_old is not null 
and anc_block_id is not null;
 --~ count 
--~ -------
 --~ 9249
 
select event_type, count(*) 
  from blocks.ancestral_block 
  inner join (
    select distinct(anc_block_id) 
      from phylogeny.event 
      inner join phylogeny.represented_event using (event_id) 
      inner join phylogeny.event2subfam using (event_id) 
      left join (
        select rec_col_id as rec_col_id_old, event_id 
          from phylogeny.represented_event 
        where rec_col_id!='2013-08-14'
      ) as oldrec using (event_id) 
    where rec_col_id='2013-08-14' 
    and event.event_type='T' 
    and origin='f' 
    and rec_col_id_old is not null
  ) as oldrec using (anc_block_id) 
group by event_type;
 --~ count 
--~ -------
 --~ 8288

-- leaf block size distribution 
-- (duplications)
--~ create temp table duplblocsizedistrib as (
  --~ select max, count(*) from (
    --~ select anc_block_id, max(count) from (
      --~ select leaf_block_id, anc_block_id, count(distinct(gene_id)) 
        --~ from blocks.leaf_block 
        --~ inner join blocks.gene2block using (leaf_block_id) 
        --~ inner join blocks.ancestral_block using (anc_block_id) 
      --~ where event_type='D' 
      --~ group by leaf_block_id, anc_block_id
    --~ ) as lb group by anc_block_id
  --~ ) as ab
  --~ group by max order by max
--~ );
--~ \copy duplblocsizedistrib to 'dup_leafblock_size_distrib.tab'

--~ -- (transfers)
--~ create temp table tralblocsizedistrib as (
  --~ select max, count(*) from (
    --~ select anc_block_id, max(count) from (
      --~ select leaf_block_id, anc_block_id, count(distinct(gene_id)) 
        --~ from blocks.leaf_block 
        --~ inner join blocks.gene2block using (leaf_block_id) 
        --~ inner join blocks.ancestral_block using (anc_block_id) 
      --~ where event_type='T' 
      --~ group by leaf_block_id, anc_block_id
    --~ ) as lb group by anc_block_id
  --~ ) as ab
  --~ group by max order by max
--~ );
--~ \copy tralblocsizedistrib to 'trans_leafblock_size_distrib.tab'
--~ 

create temp table duplblocsizedistrib as (
  select max, count(*) from (
    select anc_block_id, max(count) from (
      select leaf_block_id, anc_block_id, count(distinct(gene_id)) 
        from blocks.leaf_block 
        inner join blocks.gene2block using (leaf_block_id) 
        inner join  phylogeny.event using (anc_block_id)
        inner join phylogeny.represented_event using (event_id) 
      where event_type='D' 
      and rec_col_id='2013-08-14' 
      group by leaf_block_id, anc_block_id
    ) as lb group by anc_block_id
  ) as ab
  group by max order by max
);
\copy duplblocsizedistrib to 'dup_leafblock_size_distrib.tab'

-- (transfers)
create temp table tralblocsizedistrib as (
  select max, count(*) from (
    select anc_block_id, max(count) from (
      select leaf_block_id, anc_block_id, count(distinct(gene_id)) 
        from blocks.leaf_block 
        inner join blocks.gene2block using (leaf_block_id) 
        inner join  phylogeny.event using (anc_block_id)
        inner join phylogeny.represented_event using (event_id) 
      where event_type='T' 
      and rec_col_id='2013-08-14' 
      group by leaf_block_id, anc_block_id
    ) as lb group by anc_block_id
  ) as ab
  group by max order by max
);
\copy tralblocsizedistrib to 'trans_leafblock_size_distrib.tab'

-- (originations)
create temp table orilblocsizedistrib as (
  select max, count(*) from (
    select anc_block_id, max(count) from (
      select leaf_block_id, anc_block_id, count(distinct(gene_id)) 
        from blocks.leaf_block 
        inner join blocks.gene2block using (leaf_block_id) 
        inner join  phylogeny.event using (anc_block_id)
        inner join phylogeny.represented_event using (event_id) 
      where event_type='G' 
      and rec_col_id='2013-08-14' 
      group by leaf_block_id, anc_block_id
    ) as lb group by anc_block_id
  ) as ab
  group by max order by max
);
\copy orilblocsizedistrib to 'ori_leafblock_size_distrib.tab'

-- ancestral block size distribution 
--~ -- (duplications)
--~ create temp table dupablocsizedistrib as (
  --~ select nbfam, count(*) from (
    --~ select anc_block_id, count(family_accession) as nbfam from (
      --~ select distinct(family_accession), anc_block_id
        --~ from blocks.ancestral_block
        --~ inner join blocks.fam2block using (anc_block_id) 
      --~ where event_type='D' 
    --~ ) as lb group by anc_block_id
  --~ ) as ab
  --~ group by nbfam order by nbfam
--~ );
--~ \copy dupablocsizedistrib to 'dup_ancblock_size_distrib.tab'
--~ -- (transfers)
--~ create temp table traablocsizedistrib as (
  --~ select nbfam, count(*) from (
    --~ select anc_block_id, count(family_accession) as nbfam from (
      --~ select distinct(family_accession), anc_block_id
        --~ from blocks.ancestral_block
        --~ inner join blocks.fam2block using (anc_block_id) 
      --~ where event_type='T' 
    --~ ) as lb group by anc_block_id
  --~ ) as ab
  --~ group by nbfam order by nbfam
--~ );
--~ \copy traablocsizedistrib to 'trans_ancblock_size_distrib.tab'
--~ -- (originations)
--~ create temp table oriablocsizedistrib as (
  --~ select nbfam, count(*) from (
    --~ select anc_block_id, count(family_accession) as nbfam from (
      --~ select distinct(family_accession), anc_block_id
        --~ from blocks.ancestral_block
        --~ inner join blocks.fam2block using (anc_block_id) 
      --~ where event_type='G' 
    --~ ) as lb group by anc_block_id
  --~ ) as ab
  --~ group by nbfam order by nbfam
--~ );
--~ \copy oriablocsizedistrib to 'ori_ancblock_size_distrib.tab'
--~ 
-- (duplications)
create temp table dupablocsizedistrib as (
  select nbfam, count(*) from (
    select anc_block_id, count(family_accession) as nbfam from (
      select distinct(family_accession), anc_block_id
        from phylogeny.event 
        inner join phylogeny.represented_event using (event_id) 
        inner join blocks.fam2block using (anc_block_id) 
      where event_type='D' 
      and rec_col_id='2013-08-14' 
    ) as lb group by anc_block_id
  ) as ab
  group by nbfam order by nbfam
);
\copy dupablocsizedistrib to 'dup_ancblock_size_distrib.tab'
-- (transfers)
create temp table traablocsizedistrib as (
  select nbfam, count(*) from (
    select anc_block_id, count(family_accession) as nbfam from (
      select distinct(family_accession), anc_block_id
        from phylogeny.event 
        inner join phylogeny.represented_event using (event_id) 
        inner join blocks.fam2block using (anc_block_id) 
      where event_type='T' 
      and rec_col_id='2013-08-14' 
    ) as lb group by anc_block_id
  ) as ab
  group by nbfam order by nbfam
);
\copy traablocsizedistrib to 'trans_ancblock_size_distrib.tab'
-- (originations)
create temp table oriablocsizedistrib as (
  select nbfam, count(*) from (
    select anc_block_id, count(family_accession) as nbfam from (
      select distinct(family_accession), anc_block_id
        from phylogeny.event 
        inner join phylogeny.represented_event using (event_id) 
        inner join blocks.fam2block using (anc_block_id) 
      where event_type='G' 
      and rec_col_id='2013-08-14' 
    ) as lb group by anc_block_id
  ) as ab
  group by nbfam order by nbfam
);
\copy oriablocsizedistrib to 'ori_ancblock_size_distrib.tab'

-- blocks non investigated for block amalgamation
create temp table evtnoninvestigatedforblocks as (
  select event_type, count(*) from 
    phylogeny.event 
    inner join phylogeny.represented_event using (event_id) 
  where rec_col_id='2013-08-14' 
  and anc_block_id is null 
  group by event_type
);
\copy evtnoninvestigatedforblocks to 'event_non_investigated_for_blocks.tab'

-- large transfer events ?
select distinct anc_block_id, hogenom_gene_id, locus_tag, new_locus_tag, number, characteristic from blocks.leaf_block
        inner join blocks.ancestral_block using (anc_block_id) 
        inner join blocks.block_possible_species_node on block_id=anc_block_id 
        inner join phylogeny.species_node using (sp_node_id) 
        inner join blocks.gene2block using (leaf_block_id) 
        inner join genome.gene using (gene_id)
        left join genome.old2new_labels on locus_tag=old_locus_tag
where anc_block_id in (select anc_block_id from (
   select anc_block_id, count(family_accession) as nbfam from (
      select distinct(family_accession), anc_block_id
        from blocks.ancestral_block
        inner join blocks.fam2block using (anc_block_id) 
      where event_type='T' 
    ) as lb group by anc_block_id
  ) as ab where nbfam > 20) 
 order by anc_block_id, locus_tag;

-- large transfer events with how much prunier support?
select * from (
  select anc_block_id, rbpsn.number as rec, dbpsn.number as don, count(distinct(family_accession)) as nb_fam, nb_prunier 
    from (
      select * 
        from blocks.block_possible_species_node 
        inner join phylogeny.species_node using (sp_node_id) 
      where characteristic='rec'
    ) as rbpsn 
    inner join (
      select * 
        from blocks.block_possible_species_node 
        inner join phylogeny.species_node using (sp_node_id) 
      where characteristic='don'
    ) as dbpsn using (block_id) 
    inner join (
      select anc_block_id, ancestral_block.event_type, count(distinct(event_id)) as nb_prunier 
        from blocks.ancestral_block 
        left join ( select anc_block_id, event_id from
          phylogeny.event 
          inner join analysis.prunier_inference using (event_id) 
        ) as q using (anc_block_id)
      group by anc_block_id, ancestral_block.event_type
    ) as a on block_id=anc_block_id 
    inner join blocks.fam2block using (anc_block_id) 
  where rbpsn.block_type='ancestral' 
  and event_type='T' 
  group by anc_block_id, rbpsn.number, dbpsn.number, nb_prunier
) as n 
where nb_fam > 20 
and nb_prunier >= nb_fam/3
order by nb_fam desc;


-- frequency of allelic replcement by replicon?
create temporary table repltransloc as(
select location, number, count(*), totlocnum
  from phylogeny.event 
  inner join phylogeny.represented_event using (event_id) 
  inner join phylogeny.event2subfam using (event_id) 
  left join (
    select rec_col_id as rec_col_id_old, event_id 
      from phylogeny.represented_event 
    where rec_col_id!='2013-08-14'
  ) as oldrec using (event_id) 
  inner join phylogeny.event_possible_species_node using (event_id) 
  inner join phylogeny.species_node using (sp_node_id) 
  inner join genome.location_profile using (subfam_id, number) 
  inner join (
    select location, number, count(*) as totlocnum
    from genome.location_profile 
    group by location, number 
  ) as totloc using (location, number)
where rec_col_id='2013-08-14' 
and event.event_type='T' 
and origin='f' 
and rec_col_id_old is not null 
and characteristic='rec' 
and reference_node='t' 
group by location, number, totlocnum
order by location, number
);
\copy repltransloc to 'replacing_transfer_locations.tab'


-- frequency of clade-specific genes in clusters?
--~ select distinct specific_gene.number, nbfam, count(distinct(hogenom_gene_id)) as nbspegene
  --~ from genome.gene
    --~ inner join genome.gene2subfam using (hogenom_gene_id)
    --~ inner join genome.specific_gene using (subfam_id)
    --~ inner join blocks.gene2block using (gene_id)
    --~ inner join blocks.leaf_block using (leaf_block_id) 
    --~ inner join blocks.ancestral_block using (anc_block_id) 
    --~ inner join blocks.block_possible_species_node as bpsn on block_id=anc_block_id 
    --~ inner join phylogeny.species_node using (sp_node_id) 
    --~ inner join (
      --~ select anc_block_id, count(distinct(family_accession)) as nbfam
        --~ from blocks.ancestral_block
        --~ inner join blocks.fam2block using (anc_block_id) 
      --~ group by anc_block_id
    --~ ) as lab using (anc_block_id) 
--~ where specificity='specific_presence'
--~ and bpsn.reference_node is true
--~ and characteristic='rec'
--~ and specific_gene.number=species_node.number
--~ group by specific_gene.number, nbfam
--~ order by specific_gene.number, nbfam
--~ ;

-- coordinates of clade-specifc genes. !!! relxed criterion, presence elsewhere to be filtered
create temporary table coord_spegenes as (
select distinct species_node.number as current_genome, specific_gene.number as clade, chromosome, hogenom_gene_id, genomic_beg
  from genome.gene
    inner join phylogeny.species_node using (tax_id)
    inner join genome.gene2subfam using (hogenom_gene_id)
    inner join genome.specific_gene using (subfam_id)
where specificity='specific_presence'
and relaxed <=2 
order by species_node.number, specific_gene.number, chromosome, genomic_beg
);
\copy coord_spegenes to 'coordinates_clade-specifc_genes_relaxed2'

\copy genome.gene (chromosome, hogenom_gene_id, genomic_beg)  to 'coordinates_all_genes'
