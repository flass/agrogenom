--
-- This script present queries for exploring the structure of gene flow by HGT and gene lineage creation (by origination, duplication and HGT) and loss across the species tree
--

-- create handy view
set search_path = phylogeny,public;
create view tree_events as 
  select family_accession, rec_gi_id, reconciled_gene_tree.filepath, event_id, event_type, rec_gi_start_node, rec_gi_end_node
  from gene_tree
    inner join reconciled_gene_tree on (input_gene_tree_id=gene_tree_id)
    inner join event using (rec_gi_id)
;

-- count of all recorded events; may include alternative events, i.e. more than one event / gene tree node
select count(*) as tevents from tree_events
	where event_type='T'
), (
select count(*) as devents from tree_events
	where event_type='D'
);

-- list of all reference (i.e. optimal, final set of) transfer events, with detail of donor and recipient species tree (ancestor) nodes
select distinct on ( rec_gi_id, rec_gi_start_node ) 
family_accession, event_id, rec_gi_id, rec_gi_start_node, recnode.number as receptor, donnode.number as donor
from gene_tree
inner join phylogeny.reconciled_gene_tree on (gene_tree_id=input_gene_tree_id)
inner join phylogeny.event  using (rec_gi_id)
inner join (
	select * from phylogeny.event_possible_species_node
	inner join species_node using (sp_node_id)
		where characteristic='rec'
		and reference_node=TRUE
	) as recnode using (event_id)
inner join (
	select * from phylogeny.event_possible_species_node
	inner join species_node using (sp_node_id)
		where characteristic='don'
		and reference_node=TRUE
	) as donnode using (event_id)
order by rec_gi_id, rec_gi_start_node;


-- count of all transfer events by donor and recipient pairs
select family_accession, count( distinct event_id ), rec_gi_id, rec_gi_start_node, recnode.number as receptor, donnode.number as donor, subtree_file_path
from gene_tree
inner join phylogeny.reconciled_gene_tree on (gene_tree_id=input_gene_tree_id)
inner join phylogeny.event using (rec_gi_id)
inner join analysis.prunier_inference using (event_id)
inner join phylogeny.unicopy_subtree using (unicopy_subtree_id)
inner join (
	select * from phylogeny.event_possible_species_node
	inner join species_node using (sp_node_id)
		where characteristic='rec'
		and reference_node=TRUE
	) as recnode using (event_id)
inner join (
	select * from phylogeny.event_possible_species_node
	inner join species_node using (sp_node_id)
		where characteristic='don'
		and reference_node=TRUE
	) as donnode using (event_id)
group by family_accession, rec_gi_id, rec_gi_start_node, recnode.number, donnode.number;

select count(*) from phylogeny.event;

-- count of all transfer and duplication events by gene family
set search_path = phylogeny,public;
select devents.family_accession, tcounts, dcounts from (
	select family_accession, count(*) as tcounts from tree_events
		where event_type='T'
	group by tree_events.family_accession
) as tevents, (
	select family_accession, count(*) as dcounts from tree_events
		where event_type='D'
	group by tree_events.family_accession
	) as devents
where tevents.family_accession=devents.family_accession
group by tree_events.family_accession;


-- Example search of every DTL event occuring in a particular gene family '49RHIZOB_3475'
-- and display their potential location in the species tree (recipient for Ts)
set search_path = phylogeny,public,pg_catalog;
select * from tree_events
inner join event_possible_species_node using (event_id)
inner join species_node using (sp_node_id)
where family_accession='49RHIZOB_3475'
and characteristic='rec'
and reference_node=TRUE;
