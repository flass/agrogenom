
create temporary table diffd as (
  select * from (
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
);

\copy diffd to 'compare_duplication_event_block_location_presision.tab'

create temporary table diffg as (
  select * from (
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
);

\copy diffg to 'compare_gain_event_block_location_presision.tab'


create temporary table difft as (
  select * from (
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
);

\copy difft to 'compare_transfer_event_block_location_presision.tab'

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
select count(*) from phylogeny.event inner join phylogeny.represented_event using (event_id) inner join phylogeny.event2subfam using (event_id) left join (select rec_col_id as rec_col_id_old, event_id from phylogeny.represented_event where rec_col_id!='2013-08-14') as oldrec using (event_id) where rec_col_id='2013-08-14' and event.event_type='T' and origin='f' and rec_col_id_old is not null and anc_block_id is not null;
 --~ count 
--~ -------
 --~ 9249
 
select event_type, count(*) from blocks.ancestral_block inner join (select distinct(anc_block_id) from phylogeny.event inner join phylogeny.represented_event using (event_id) inner join phylogeny.event2subfam using (event_id) left join (select rec_col_id as rec_col_id_old, event_id from phylogeny.represented_event where rec_col_id!='2013-08-14') as oldrec using (event_id) where rec_col_id='2013-08-14' and event.event_type='T' and origin='f' and rec_col_id_old is not null) as oldrec using (anc_block_id) group by event_type;
 --~ count 
--~ -------
 --~ 8288

-- ancestral block size distribution (transfers)
create temp table dupblocsizeditrib as (
  select max, count(*) from (
    select anc_block_id, max(count) from (
      select leaf_block_id, anc_block_id, count(distinct(gene_id)) 
        from blocks.leaf_block 
        inner join blocks.gene2block using (leaf_block_id) 
        inner join blocks.ancestral_block using (anc_block_id) 
      where event_type='D' 
      group by leaf_block_id, anc_block_id
    ) as lb group by anc_block_id
  ) as ab
  group by max order by max
);

\copy dupblocsizeditrib to 'dup_block_size_ditrib.tab'

create temp table trablocsizeditrib as (
  select max, count(*) from (
    select anc_block_id, max(count) from (
      select leaf_block_id, anc_block_id, count(distinct(gene_id)) 
        from blocks.leaf_block 
        inner join blocks.gene2block using (leaf_block_id) 
        inner join blocks.ancestral_block using (anc_block_id) 
      where event_type='T' 
      group by leaf_block_id, anc_block_id
    ) as lb group by anc_block_id
  ) as ab
  group by max order by max
);
\copy trablocsizeditrib to 'trans_block_size_ditrib.tab'
