--
-- Agrgenom PostgreSQL database dump
--

SET statement_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = off;
SET check_function_bodies = false;
SET client_min_messages = warning;
SET escape_string_warning = off;

--
-- as postgres super_user
--

CREATE EXTENSION ltree;

--
-- Name: analysis; Type: SCHEMA; Schema: -; Owner: agrogenomadmin
--

CREATE SCHEMA analysis;


ALTER SCHEMA analysis OWNER TO agrogenomadmin;

--
-- Name: SCHEMA analysis; Type: COMMENT; Schema: -; Owner: agrogenomadmin
--

COMMENT ON SCHEMA analysis IS 'Store information about analysis . TO COMPLETE';


--
-- Name: genome; Type: SCHEMA; Schema: -; Owner: agrogenomadmin
--

CREATE SCHEMA genome;


ALTER SCHEMA genome OWNER TO agrogenomadmin;

--
-- Name: SCHEMA genome; Type: COMMENT; Schema: -; Owner: agrogenomadmin
--

COMMENT ON SCHEMA genome IS 'Store information about genome. TO COMPLETE';


--
-- Name: phylogeny; Type: SCHEMA; Schema: -; Owner: agrogenomadmin
--

CREATE SCHEMA phylogeny;


ALTER SCHEMA phylogeny OWNER TO agrogenomadmin;

--
-- Name: SCHEMA phylogeny; Type: COMMENT; Schema: -; Owner: agrogenomadmin
--

COMMENT ON SCHEMA phylogeny IS 'Store information about phylogeny. TO COMPLETE';


--
-- Name: search; Type: SCHEMA; Schema: -; Owner: agrogenomadmin
--

CREATE SCHEMA search;


ALTER SCHEMA search OWNER TO agrogenomadmin;

--
-- Name: SCHEMA search; Type: COMMENT; Schema: -; Owner: agrogenomadmin
--

COMMENT ON SCHEMA search IS 'Store information that will be requested a lot.';


--
-- Name: taxonomy; Type: SCHEMA; Schema: -; Owner: agrogenomadmin
--

CREATE SCHEMA taxonomy;


ALTER SCHEMA taxonomy OWNER TO agrogenomadmin;

--
-- Name: SCHEMA taxonomy; Type: COMMENT; Schema: -; Owner: agrogenomadmin
--

COMMENT ON SCHEMA taxonomy IS 'Store information about taxonomy. The taxonomy schema is based on the information provided by the NCBI in readme.txt and in ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz TO COMPLETE';


--
-- Name: blocks; Type: SCHEMA; Schema: -; Owner: agrogenomadmin
--

CREATE SCHEMA blocks;


ALTER SCHEMA blocks OWNER TO agrogenomadmin;

--
-- Name: SCHEMA taxonomy; Type: COMMENT; Schema: -; Owner: agrogenomadmin
--

COMMENT ON SCHEMA blocks IS 'Store information relative to block of neighbouring genes in genomes sharing common events in their respecitive phylogenetic trees';


--
-- Name: plpgsql; Type: PROCEDURAL LANGUAGE; Schema: -; Owner: postgres
--

CREATE PROCEDURAL LANGUAGE plpgsql;


ALTER PROCEDURAL LANGUAGE plpgsql OWNER TO postgres;

SET search_path = public, pg_catalog;

--
-- Name: analysis; Type: TYPE; Schema: public; Owner: agrogenomadmin
--

CREATE TYPE analysis AS ENUM (
	'REC',
	'G_T',
	'Sp_T',
	'ROOT'
);


ALTER TYPE public.analysis OWNER TO agrogenomadmin;

--
-- Name: eventnodecharacteristic; Type: TYPE; Schema: public; Owner: agrogenomadmin
--
--~ 
--~ CREATE TYPE eventnodecharacteristic AS ENUM (
	--~ 'start',
	--~ 'end',
	--~ 'null'
--~ );

CREATE TYPE eventnodecharacteristic AS ENUM (
	'rec',
	'don',
	'location'
);


ALTER TYPE public.eventnodecharacteristic OWNER TO agrogenomadmin;

--
-- Name: evolutionaryevent; Type: TYPE; Schema: public; Owner: agrogenomadmin
--

CREATE TYPE evolutionaryevent AS ENUM (
	'D',
	'L',
	'T',
	'S',
	'G'
);


ALTER TYPE public.evolutionaryevent OWNER TO agrogenomadmin;

--
-- Name: treestatus; Type: TYPE; Schema: public; Owner: agrogenomadmin
--

CREATE TYPE treestatus AS ENUM (
	'I',
	'M'
);


ALTER TYPE public.treestatus OWNER TO agrogenomadmin;

--
-- Name: blocktype; Type: TYPE; Schema: public; Owner: agrogenomadmin
--

CREATE TYPE blocktype AS ENUM (
	'ancestral',
	'leaf'
);


ALTER TYPE public.blocktype OWNER TO agrogenomadmin;

--
-- Name: prunier_status; Type: TYPE; Schema: public; Owner: agrogenomadmin
--

--~ CREATE TYPE prunier_status AS ENUM (
	--~ 'not_done',
	--~ 'interupted',
	--~ 'done',
	--~ 'failed'
--~ );
--~ 
--~ 
--~ ALTER TYPE public.prunier_status OWNER TO agrogenomadmin;

SET search_path = public, pg_catalog;


CREATE FUNCTION get_rec_gi_id(text) RETURNS integer
	LANGUAGE plpgsql
	AS $_$
declare
		num alias for $1;
		rec_gi_id integer;
begin
		select into rec_gi_id rec_gi_id from phylogeny.reconciled_gene_tree where species_tree_id=num;
return rec_gi_id;
end;
$_$;


ALTER FUNCTION public.get_rec_gi_id(text) OWNER TO agrogenomadmin;

--
-- Name: get_rec_gi_id(integer); Type: FUNCTION; Schema: public; Owner: agrogenomadmin
--

CREATE FUNCTION get_rec_gi_id(integer) RETURNS integer
	LANGUAGE plpgsql
	AS $_$
declare
		num alias for $1;
		rec_gi_id integer;
begin
		select into rec_gi_id rec_gi_id from phylogeny.reconciled_gene_tree where species_tree_id=num;
return rec_gi_id;
end;
$_$;


ALTER FUNCTION public.get_rec_gi_id(integer) OWNER TO agrogenomadmin;

--
-- Name: get_species_tree_id(text); Type: FUNCTION; Schema: public; Owner: agrogenomadmin
--

CREATE FUNCTION get_species_tree_id(text) RETURNS integer
	LANGUAGE plpgsql
	AS $_$
declare
	name_file_path alias for $1;
	species_tree_id integer;
begin
	select into species_tree_id species_tree_id from phylogeny.species_tree where file_path=name_file_path;
return species_tree_id;
end;
$_$;


ALTER FUNCTION public.get_species_tree_id(text) OWNER TO agrogenomadmin;


SET search_path = taxonomy, public, pg_catalog;

--
-- Name: findtaxon(text); Type: FUNCTION; Schema: taxonomy; Owner: agrogenomadmin
--

CREATE FUNCTION findtaxon(text) RETURNS integer
	LANGUAGE sql STABLE STRICT
	AS $_$
select distinct taxonomy.nodes.tax_id
from taxonomy.nodes,
taxonomy.names
where taxonomy.nodes.tax_id = taxonomy.names.tax_id
and taxonomy.names.name_txt = $1;
$_$;


ALTER FUNCTION taxonomy.findtaxon(text) OWNER TO agrogenomadmin;

--
-- Name: getspecies(integer); Type: FUNCTION; Schema: taxonomy; Owner: agrogenomadmin
--

CREATE FUNCTION getspecies(integer) RETURNS integer
	LANGUAGE sql STABLE STRICT
	AS $_$
select answer.tax_id
from taxonomy.nodes as query , 
taxonomy.nodes as answer 
where answer.rank='species' 
and answer.path @> query.path 
and query.tax_id= $1
order by answer.path;
$_$;


ALTER FUNCTION taxonomy.getspecies(integer) OWNER TO agrogenomadmin;

--
-- Name: getsuperkingdom(integer); Type: FUNCTION; Schema: taxonomy; Owner: agrogenomadmin
--

CREATE FUNCTION getsuperkingdom(integer) RETURNS integer
	LANGUAGE sql STABLE STRICT
	AS $_$
select answer.tax_id
from taxonomy.nodes as query , 
taxonomy.nodes as answer 
where answer.rank='superkingdom' 
and answer.path @> query.path 
and query.tax_id= $1
$_$;


ALTER FUNCTION taxonomy.getsuperkingdom(integer) OWNER TO agrogenomadmin;

--
-- Name: getthirdlevel(integer); Type: FUNCTION; Schema: taxonomy; Owner: agrogenomadmin
--

CREATE FUNCTION getthirdlevel(id integer) RETURNS integer
	LANGUAGE sql
	AS $_$
select answer.tax_id
from taxonomy.nodes as query , 
taxonomy.nodes as answer 
where
answer.path @> query.path 
and nlevel(answer.path) = 3
and query.tax_id= $1
$_$;


ALTER FUNCTION taxonomy.getthirdlevel(id integer) OWNER TO agrogenomadmin;

--
-- Name: isunderof(integer, integer); Type: FUNCTION; Schema: taxonomy; Owner: agrogenomadmin
--

CREATE FUNCTION isunderof(integer, integer) RETURNS boolean
	LANGUAGE sql STABLE STRICT
	AS $_$
select ancetre.path @> noeud.path 
from taxonomy.nodes as ancetre,
taxonomy.nodes as noeud
where ancetre.tax_id = $1 
and noeud.tax_id = $2;
$_$;


ALTER FUNCTION taxonomy.isunderof(integer, integer) OWNER TO agrogenomadmin;

--
-- Name: path(integer); Type: FUNCTION; Schema: taxonomy; Owner: agrogenomadmin
--

CREATE FUNCTION path(integer) RETURNS text
	LANGUAGE plpgsql STABLE STRICT
	AS $_$
	 
declare
   taxa alias for $1;
   rep  text;
   rec  record;
begin
   rep = '';
   for rec in select taxonomy.names.name_txt 
				from taxonomy.names,
					 taxonomy.nodes as query, 
					 taxonomy.nodes as answer 
			   where answer.path @> query.path 
				 and answer.genbank_hidden_flag!=1
				 and taxonomy.names.tax_id = answer.tax_id 
				 and taxonomy.names.name_class = 'scientific name' 
				 and answer.tax_id > 1 
				 and query.tax_id != answer.tax_id 
				 and query.tax_id = taxa 
			order by nlevel(answer.path)
   loop
	 rep := rep || ':' || rec.name_txt;
   end loop;
   return substring(rep from 2);
   end;
 $_$;


ALTER FUNCTION taxonomy.path(integer) OWNER TO agrogenomadmin;

--
-- Name: scientificname(integer); Type: FUNCTION; Schema: taxonomy; Owner: agrogenomadmin
--

CREATE FUNCTION scientificname(integer) RETURNS text
	LANGUAGE sql STABLE STRICT
	AS $_$
select taxonomy.names.name_txt
from taxonomy.names
where taxonomy.names.tax_id = $1
and taxonomy.names.name_class ='scientific name';
$_$;


ALTER FUNCTION taxonomy.scientificname(integer) OWNER TO agrogenomadmin;

--
-- Name: synonymousname(integer); Type: FUNCTION; Schema: taxonomy; Owner: agrogenomadmin
--

CREATE FUNCTION synonymousname(integer) RETURNS text
	LANGUAGE plpgsql STABLE STRICT
	AS $_$
declare
   taxa alias for $1;
   rep  text;
   rec  record;
begin
   rep = '';
   for rec in select taxonomy.names.name_txt
		  from taxonomy.names
		  where taxonomy.names.tax_id = taxa
			and taxonomy.names.name_class in ('synonym',
										'common name',
										'equivalent name')
   loop
	 rep := rep || rec.name_txt || ';' ;
   end loop;
   return trim(trailing ';' from rep);
end;
$_$;


ALTER FUNCTION taxonomy.synonymousname(integer) OWNER TO agrogenomadmin;


SET default_tablespace = '';

SET default_with_oids = false;

SET search_path = analysis, pg_catalog;

--
-- Name: analysis; Type: TABLE; Schema: analysis; Owner: agrogenomadmin; Tablespace: 
--


CREATE TABLE analysis (
	analysis_id serial PRIMARY KEY,
	parameter_id integer NOT NULL,
	prog_id integer NOT NULL,
	type public.analysis,
	comments text
);


ALTER TABLE analysis.analysis OWNER TO agrogenomadmin;

--
-- Name: TABLE analysis; Type: COMMENT; Schema: analysis; Owner: agrogenomadmin
--

COMMENT ON TABLE analysis IS 'TO COMPLETE';


--
-- Name: cost; Type: TABLE; Schema: analysis; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE cost (
	cost_id serial PRIMARY KEY,
	cost_value double precision NOT NULL,
	eventtype public.evolutionaryevent NOT NULL,
	analysis_id integer NOT NULL
);

ALTER TABLE analysis.cost OWNER TO agrogenomadmin;


--
-- Name: inference_parameter; Type: TABLE; Schema: analysis; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE inference_parameter (
	parameter_id serial PRIMARY KEY,
	parameter_string text
);

ALTER TABLE analysis.inference_parameter OWNER TO agrogenomadmin;

--
-- Name: prunier_inference; Type: TABLE; Schema: analysis; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE prunier_inference (
	prinf_id serial PRIMARY KEY,
	unicopy_subtree_id integer NOT NULL,
	prunier_round smallint NOT NULL,
	event_id integer,
	analysis_id integer NOT NULL,
	inverted_rec_don bool
);

ALTER TABLE analysis.prunier_inference OWNER TO agrogenomadmin;

--
-- Name: prunier_inference_support; Type: TABLE; Schema: analysis; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE prunier_inference_support (
	prinf_id integer NOT NULL,
	branch_support double precision NOT NULL
);


ALTER TABLE analysis.prunier_inference_support OWNER TO agrogenomadmin;

--
-- Name: pruned_leaf; Type: TABLE; Schema: analysis; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE pruned_leaf (
	prinf_id integer NOT NULL,
	gene_id integer NOT NULL
);


ALTER TABLE analysis.pruned_leaf OWNER TO agrogenomadmin;

--
-- Name: program; Type: TABLE; Schema: analysis; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE program (
	prog_id serial PRIMARY KEY,
	name varchar(50) NOT NULL,
	version real NOT NULL
);

ALTER TABLE analysis.program OWNER TO agrogenomadmin;


SET search_path = genome, pg_catalog;

--
-- Name: gene; Type: TABLE; Schema: genome; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE gene (
	gene_id serial PRIMARY KEY,
	tax_id integer,
	previous_5_gene_id integer,
	following_3_gene_id integer,
	family_accession varchar(50),
	family_type_id smallint,
	name varchar(500) default NULL,
	locus_tag varchar(500) default NULL,
	chromosome varchar(500),
	genomic_beg integer,
	genomic_end integer,
	strand smallint NOT NULL,
	description text,
	hogenom_gene_id varchar(50) NOT NULL
);


ALTER TABLE genome.gene OWNER TO agrogenomadmin;


CREATE TABLE chromosome (
	chromosome varchar(500) PRIMARY KEY,
	tax_id integer NOT NULL,
	topology varchar(500)
);


ALTER TABLE genome.chromosome OWNER TO agrogenomadmin;


CREATE TABLE genome.phylogenetic_profile (
subfam_id varchar(50) NOT NULL,
number varchar(50) NOT NULL,
present boolean,
count smallint,
PRIMARY KEY (subfam_id, number)
);
ALTER TABLE genome.phylogenetic_profile OWNER TO agrogenomadmin;
CREATE INDEX phylogetic_profile_subfam_id_key ON genome.phylogenetic_profile (subfam_id);
CREATE INDEX phylogetic_profile_number_key ON genome.phylogenetic_profile (number);

CREATE TABLE genome.specific_gene (
subfam_id varchar(50) NOT NULL,
number varchar(50) NOT NULL,
specificity varchar(20) NOT NULL,
relaxed smallint,
PRIMARY KEY (subfam_id, number, specificity)
);
ALTER TABLE genome.specific_gene OWNER TO agrogenomadmin;
CREATE INDEX specific_genes_subfam_id_key ON genome.specific_gene (subfam_id);
CREATE INDEX specific_genes_number_key ON genome.specific_gene (number);
CREATE INDEX specific_genes_specificity_key ON genome.specific_gene (specificity);

CREATE TABLE genome.location_profile (
subfam_id varchar(50) NOT NULL,
number varchar(50) NOT NULL,
location varchar(10),
PRIMARY KEY (subfam_id, number)
);
ALTER TABLE genome.location_profile OWNER TO agrogenomadmin;
CREATE INDEX location_profile_subfam_id_key ON genome.location_profile (subfam_id);
CREATE INDEX location_profile_number_key ON genome.location_profile (number);


-- phylogeny schema tables

SET search_path = phylogeny, pg_catalog;

--
-- Name: event; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE event (
	event_id serial PRIMARY KEY,
	event_type public.evolutionaryevent NOT NULL,
	rec_gi_id integer NOT NULL,
	inf_event_date integer,
	sup_event_date integer,
	support double precision,
	anc_block_id bigint NOT NULL,
	rec_gi_start_node integer NOT NULL,
	rec_gi_end_node integer
);

ALTER TABLE phylogeny.event OWNER TO agrogenomadmin;

--
-- Name: represented_event; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE represented_event (
	event_id INTEGER NOT NULL,
	nr_node_id INTEGER NOT NULL,
	rec_col_id DATE NOT NULL,
	PRIMARY KEY (nr_node_id, rec_col_id)
);

ALTER TABLE phylogeny.represented_event OWNER TO agrogenomadmin;


--
-- Name: reconcilations; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE reconciliation_collection (
	rec_col_id DATE PRIMARY KEY,
	rec_tree_directory_path VARCHAR(500)
);

ALTER TABLE phylogeny.reconciliation_collection OWNER TO agrogenomadmin;

--
-- Name: event_per_node; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE event_possible_species_node (
	event_id integer NOT NULL,
	sp_node_id integer NOT NULL,
	characteristic public.eventnodecharacteristic,
	reference_node bool NOT NULL,
	UNIQUE (event_id, sp_node_id, characteristic)
);

ALTER TABLE phylogeny.event_possible_species_node OWNER TO agrogenomadmin;

--
-- Name: family; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE family (
	accession varchar(50) NOT NULL,
	family_type_id integer NOT NULL
);

ALTER TABLE phylogeny.family OWNER TO agrogenomadmin;


--
-- Name: family_type; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE family_type (
	family_type_id serial NOT NULL,
	name varchar(256) NOT NULL,
	version smallint NOT NULL
);

ALTER TABLE phylogeny.family_type OWNER TO agrogenomadmin;

--
-- Name: gene_tree; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE gene_tree (
	gene_tree_id serial PRIMARY KEY,
	file_path text NOT NULL,
	family_accession varchar(50) NOT NULL,
	family_type_id smallint,
	type public.treestatus DEFAULT 'I',
	analysis_id integer,
	source_gene_tree_id integer
);

ALTER TABLE phylogeny.gene_tree OWNER TO agrogenomadmin;


--
-- Name: reconciled_gene_tree; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

--~ CREATE TABLE reconciled_gene_tree (
	--~ rec_gi_id integer NOT NULL,
	--~ filepath text NOT NULL,
	--~ rec_cost double precision NOT NULL,
	--~ input_gene_tree_id integer NOT NULL,
	--~ species_tree_id integer NOT NULL,
	--~ analysis_id integer NOT NULL,
	--~ tree_num integer
--~ );
CREATE TABLE reconciled_gene_tree (
	rec_gi_id serial PRIMARY KEY,
	filepath text,
	rec_cost double precision,
	input_gene_tree_id integer NOT NULL,
	species_tree_id integer NOT NULL,
	analysis_id integer NOT NULL,
	tree_num integer
);

ALTER TABLE phylogeny.reconciled_gene_tree OWNER TO agrogenomadmin;
--
-- Name: reconciled_gene_tree_node; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE reconciled_gene_tree_node (
	nr_node_id serial PRIMARY KEY,
	rec_gi_id integer NOT NULL,
	node_id integer NOT NULL,
	branch_support double precision,
	branch_length double precision
);

ALTER TABLE phylogeny.reconciled_gene_tree_node OWNER TO agrogenomadmin;

--
-- Name: unicopy_subtree; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE unicopy_subtree (
	unicopy_subtree_id serial PRIMARY KEY,
	subtree_name varchar(30) NOT NULL,
	subtree_file_path text NOT NULL,
	prunier_outfile_path text
);

ALTER TABLE phylogeny.unicopy_subtree OWNER TO agrogenomadmin;

--
-- Name: node2inference; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE node2inference (
	nr_node_id integer NOT NULL,
	prinf_id integer NOT NULL
);

ALTER TABLE phylogeny.node2inference OWNER TO agrogenomadmin;

--
-- Name: node2subtree; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE node2subtree (
	nr_node_id integer NOT NULL,
	unicopy_subtree_id integer NOT NULL
);

ALTER TABLE phylogeny.node2subtree OWNER TO agrogenomadmin;

--
-- Name: gene2subtree; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE gene2subtree (
	gene_id integer NOT NULL,
	unicopy_subtree_id integer NOT NULL
);

ALTER TABLE phylogeny.gene2subtree OWNER TO agrogenomadmin;

--
-- Name: subtree_implantation; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE subtree_implantation (
	unicopy_subtree_id integer NOT NULL,
	prunier_congruent boolean
);

ALTER TABLE phylogeny.subtree_implantation OWNER TO agrogenomadmin;

--
-- Name: species_node; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE species_node (
	sp_node_id serial PRIMARY KEY,
	number varchar(50),
	parent_sp_node_id integer,
	species_tree_id integer,
	tax_id integer,
	support double precision,
	uniprot_id varchar(50),
	inf_date integer,
	sup_date integer,
	left_num integer,
	right_num integer,
	branch_length double precision,
	UNIQUE(species_tree_id, right_num, left_num)
);


ALTER TABLE phylogeny.species_node OWNER TO agrogenomadmin;

--
-- Name: species_tree; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE species_tree (
	species_tree_id serial NOT NULL,
	name text NOT NULL,
	file_path text NOT NULL,
	analysis_id integer NOT NULL,
	source_species_tree_id integer
);


ALTER TABLE phylogeny.species_tree OWNER TO agrogenomadmin;

--
-- Name: substitution; Type: TABLE; Schema: phylogeny; Owner: agrogenomadmin; Tablespace: 
--


CREATE TABLE substitution (
  nr_node_id BIGINT NOT NULL, 
  dN float NOT NULL, 
  dS float NOT NULL, 
  analysis_id DATE, 
  PRIMARY KEY (nr_node_id, analysis_id)
);

ALTER TABLE phylogeny.substitution OWNER TO agrogenomadmin;




-- search schema tables 

SET search_path = search, pg_catalog;

--
-- Name: lineage; Type: TABLE; Schema: search; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE lineage (
	tax_id integer,
	superkingdom_l text,
	phylum_l text,
	class_l text,
	order_l text,
	family_l text,
	genus_l text,
	species_l text
);


ALTER TABLE search.lineage OWNER TO agrogenomadmin;

-- taxonomy schema tables

SET search_path = taxonomy, pg_catalog;

--
-- Name: delnodes; Type: TABLE; Schema: taxonomy; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE delnodes (
	tax_id integer NOT NULL
);


ALTER TABLE taxonomy.delnodes OWNER TO agrogenomadmin;

--
-- Name: division; Type: TABLE; Schema: taxonomy; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE division (
	division_id integer NOT NULL,
	division_cde character(3),
	division_name text,
	comments text
);


ALTER TABLE taxonomy.division OWNER TO agrogenomadmin;

--
-- Name: gencode; Type: TABLE; Schema: taxonomy; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE gencode (
	genetic_code_id varchar(3) NOT NULL,
	abbreviation text,
	name text,
	cde text,
	starts text
);


ALTER TABLE taxonomy.gencode OWNER TO agrogenomadmin;

--
-- Name: merged; Type: TABLE; Schema: taxonomy; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE merged (
	old_tax_id integer NOT NULL,
	new_tax_id integer
);


ALTER TABLE taxonomy.merged OWNER TO agrogenomadmin;

--
-- Name: names; Type: TABLE; Schema: taxonomy; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE names (
	tax_id integer,
	name_txt text,
	unique_name text,
	name_class text
);


ALTER TABLE taxonomy.names OWNER TO agrogenomadmin;

--
-- Name: nodes; Type: TABLE; Schema: taxonomy; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE nodes (
	tax_id integer NOT NULL,
	parent_tax_id integer NOT NULL,
	rank varchar(50) NOT NULL,
	embl_code character(2),
	division_id integer,
	inherited_div_flag boolean,
	genetic_code_id varchar(3),
	inherited_gc_flag boolean,
	mitochondrial_genetic_code_id varchar(3),
	inherited_mgc_flag boolean,
	genbank_hidden_flag boolean,
	hidden_subtree_root_flag boolean,
	comments text,
	path public.ltree NOT NULL
);


ALTER TABLE taxonomy.nodes OWNER TO agrogenomadmin;

-- blocks schema tables

SET search_path = blocks, pg_catalog;

--
-- Name: ancestral_block; Type: TABLE; Schema: blocks; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE ancestral_block (
	anc_block_id bigint PRIMARY KEY,
	event_type public.evolutionaryevent NOT NULL
);

ALTER TABLE blocks.ancestral_block OWNER TO agrogenomadmin;

--
-- Name: fam2block; Type: TABLE; Schema: blocks; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE fam2block (
	anc_block_id bigint NOT NULL,
	family_accession varchar(50) NOT NULL
);

ALTER TABLE blocks.fam2block OWNER TO agrogenomadmin;

--
-- Name: leaf_blocks; Type: TABLE; Schema: blocks; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE leaf_block (
	leaf_block_id bigint PRIMARY KEY,
	anc_block_id bigint
);

ALTER TABLE blocks.leaf_block OWNER TO agrogenomadmin;

--
-- Name: gene2block; Type: TABLE; Schema: blocks; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE gene2block (
	leaf_block_id bigint NOT NULL,
	gene_id bigint NOT NULL
);

ALTER TABLE blocks.gene2block OWNER TO agrogenomadmin;


--
-- Name: block_possible_species_node; Type: TABLE; Schema: blocks; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE block_possible_species_node (
	block_id bigint NOT NULL,
	sp_node_id integer NOT NULL,
	characteristic public.eventnodecharacteristic,
	block_type public.blocktype,
	reference_node bool NOT NULL,
	PRIMARY KEY (block_id, block_type, sp_node_id, characteristic)
);

ALTER TABLE blocks.block_possible_species_node OWNER TO agrogenomadmin;


--
-- Name: functional_homogeneity; Type: TABLE; Schema: blocks; Owner: agrogenomadmin; Tablespace: 
--

CREATE TABLE functional_homogeneity (
	leaf_block_id bigint NOT NULL,
	aspect varchar(30) NOT NULL, 
	funsim float NOT NULL,
	paiwise_metric varchar(30) NOT NULL,
	group_metric varchar(30) NOT NULL,
	PRIMARY KEY (leaf_block_id, aspect, paiwise_metric, group_metric)
	);

ALTER TABLE blocks.functional_homogeneity OWNER TO agrogenomadmin;

--
-- CREATE TABLE INDEXES
--

CREATE UNIQUE INDEX gene_tree_family_accession_key ON phylogeny.gene_tree (family_accession);

CREATE INDEX event_event_type_key ON phylogeny.event (event_type);
CREATE INDEX event_rec_gi_id_key ON phylogeny.event (rec_gi_id);
CREATE INDEX event_rec_gi_start_node_key ON phylogeny.event (rec_gi_start_node);
CREATE INDEX event_rec_gi_end_node_key ON phylogeny.event (rec_gi_end_node);

CREATE INDEX epsn_sp_node_id_key ON phylogeny.event_possible_species_node (sp_node_id);
CREATE INDEX epsn_event_id_key ON phylogeny.event_possible_species_node (event_id);
CREATE INDEX epsn_characteristic_key ON phylogeny.event_possible_species_node (characteristic);
CREATE INDEX epsn_reference_node_key ON phylogeny.event_possible_species_node (reference_node);

CREATE UNIQUE INDEX species_node_number_key ON phylogeny.species_node (number);

CREATE INDEX g2st_gene_id_key ON phylogeny.gene2subtree (gene_id);
CREATE INDEX g2st_unicopy_subtree_id_key ON phylogeny.gene2subtree (unicopy_subtree_id);

CREATE INDEX n2st_nr_node_id_key ON phylogeny.node2subtree (nr_node_id);
CREATE INDEX n2st_unicopy_subtree_id_key ON phylogeny.node2subtree (unicopy_subtree_id);

CREATE INDEX n2i_nr_node_id_key ON phylogeny.node2inference (nr_node_id);
CREATE INDEX n2i_prinf_id_key ON phylogeny.node2inference (prinf_id);

CREATE INDEX subtree_implantation_unicopy_subtree_id_key ON phylogeny.subtree_implantation (unicopy_subtree_id);

CREATE INDEX pruned_leaf_gene_id_key ON analysis.pruned_leaf (gene_id);
CREATE INDEX pruned_leaf_prinf_id_key ON analysis.pruned_leaf (prinf_id);

CREATE INDEX prunier_inference_event_id_key ON analysis.prunier_inference (event_id);
CREATE INDEX prunier_inference_unicopy_subtree_id_key ON analysis.prunier_inference (unicopy_subtree_id);
CREATE UNIQUE INDEX prunier_inference_inference_key ON analysis.prunier_inference (unicopy_subtree_id, prunier_round);

CREATE INDEX pis_prinf_id_key ON analysis.prunier_inference_support (prinf_id);

ALTER TABLE genome.gene ADD CONSTRAINT gene_hogenom_gene_id_key UNIQUE (hogenom_gene_id);
CREATE INDEX gene_locus_tag_key ON genome.gene (locus_tag);
CREATE INDEX gene_name_key ON genome.gene (name);
CREATE INDEX gene_beg_key ON genome.gene (genomic_beg);
CREATE INDEX gene_end_key ON genome.gene (genomic_end);
CREATE INDEX gene_chromosome_beg_key ON genome.gene (chromosome, genomic_beg);
CLUSTER genome.gene USING gene_chromosome_beg_key;

ALTER TABLE phylogeny.reconciled_gene_tree_node ADD CONSTRAINT reconciled_gene_tree_node_rec_gi_id_node_id_key UNIQUE (rec_gi_id, node_id);

--
-- CREATE FOREING KEY CONSTRAINTS
--

-- depending on phylogeny.unicopy_subtree
ALTER TABLE phylogeny.subtree_implantation ADD FOREIGN KEY (unicopy_subtree_id) REFERENCES phylogeny.unicopy_subtree (unicopy_subtree_id) ON DELETE CASCADE ;
ALTER TABLE phylogeny.gene2subtree ADD FOREIGN KEY (unicopy_subtree_id) REFERENCES phylogeny.unicopy_subtree (unicopy_subtree_id) ON DELETE CASCADE ;
ALTER TABLE phylogeny.node2subtree ADD FOREIGN KEY (unicopy_subtree_id) REFERENCES phylogeny.unicopy_subtree (unicopy_subtree_id) ON DELETE CASCADE ;
ALTER TABLE analysis.prunier_inference ADD FOREIGN KEY (unicopy_subtree_id) REFERENCES phylogeny.unicopy_subtree (unicopy_subtree_id) ON DELETE CASCADE ;

-- depending on phylogeny.gene_tree
ALTER TABLE phylogeny.reconciled_gene_tree ADD FOREIGN KEY (input_gene_tree_id) REFERENCES phylogeny.gene_tree (gene_tree_id) ON DELETE CASCADE ;

-- depending on phylogeny.reconciled_gene_tree
ALTER TABLE phylogeny.reconciled_gene_tree_node ADD FOREIGN KEY (rec_gi_id) REFERENCES phylogeny.reconciled_gene_tree (rec_gi_id) ON DELETE CASCADE ;

-- depending on phylogeny.reconciled_gene_tree_node
ALTER TABLE phylogeny.event ADD FOREIGN KEY (rec_gi_id, rec_gi_start_node) REFERENCES phylogeny.reconciled_gene_tree_node (rec_gi_id, node_id) ON DELETE RESTRICT ;
--~ ALTER TABLE phylogeny.node2gene ADD FOREIGN KEY (nr_node_id) REFERENCES phylogeny.reconciled_gene_tree_node (nr_node_id) ON DELETE CASCADE ;
ALTER TABLE phylogeny.node2inference ADD FOREIGN KEY (nr_node_id) REFERENCES phylogeny.reconciled_gene_tree_node (nr_node_id) ON DELETE CASCADE ;
ALTER TABLE phylogeny.node2subtree ADD FOREIGN KEY (nr_node_id) REFERENCES phylogeny.reconciled_gene_tree_node (nr_node_id) ON DELETE CASCADE ;
ALTER TABLE phylogeny.substitution ADD FOREIGN KEY (nr_node_id) REFERENCES phylogeny.reconciled_gene_tree_node (nr_node_id) ON DELETE RESTRICT ;

-- depending on phylogeny.event
ALTER TABLE phylogeny.event_possible_species_node ADD FOREIGN KEY (event_id) REFERENCES phylogeny.event (event_id) ON DELETE CASCADE ;
ALTER TABLE  phylogeny.represented_event ADD FOREIGN KEY (event_id) REFERENCES phylogeny.event (event_id) ON DELETE CASCADE ;

-- depending on phylogeny.reconciliation_collection
ALTER TABLE phylogeny.represented_event  ADD FOREIGN KEY (rec_col_id) REFERENCES phylogeny.reconciliation_collection (rec_col_id) ON DELETE CASCADE ;

-- depending on analysis.prunier_inference
ALTER TABLE analysis.prunier_inference_support ADD FOREIGN KEY (prinf_id) REFERENCES analysis.prunier_inference (prinf_id) ON DELETE CASCADE ;
ALTER TABLE analysis.pruned_leaf ADD FOREIGN KEY (prinf_id) REFERENCES analysis.prunier_inference (prinf_id) ON DELETE CASCADE ;
ALTER TABLE phylogeny.node2inference ADD FOREIGN KEY (prinf_id) REFERENCES analysis.prunier_inference (prinf_id) ON DELETE CASCADE ;

-- depending on phylogeny.species_node
ALTER TABLE phylogeny.event_possible_species_node ADD FOREIGN KEY (sp_node_id) REFERENCES phylogeny.species_node (sp_node_id) ON DELETE RESTRICT ;
ALTER TABLE blocks.block_possible_species_node ADD FOREIGN KEY (sp_node_id) REFERENCES phylogeny.species_node (sp_node_id) ON DELETE RESTRICT ;
ALTER TABLE genome.phylogenetic_profile ADD FOREIGN KEY (number) REFERENCES phylogeny.species_node (number) ON DELETE RESTRICT ;
ALTER TABLE genome.location_profile ADD FOREIGN KEY (number) REFERENCES phylogeny.species_node (number) ON DELETE RESTRICT ;
ALTER TABLE genome.specific_gene ADD FOREIGN KEY (number) REFERENCES phylogeny.species_node (number) ON DELETE RESTRICT ;

-- depending on blocks.ancestral_block
ALTER TABLE blocks.leaf_block ADD FOREIGN KEY (anc_block_id) REFERENCES blocks.ancestral_block (anc_block_id) ON DELETE RESTRICT ;
--~ ALTER TABLE blocks.block_possible_species_node ADD FOREIGN KEY (block_id, block_type) REFERENCES blocks.ancestral_block (anc_block_id, 'ancestral') ON DELETE CASCADE ;
ALTER TABLE blocks.fam2block ADD FOREIGN KEY (anc_block_id) REFERENCES blocks.ancestral_block (anc_block_id) ON DELETE CASCADE ;
ALTER TABLE phylogeny.event ADD FOREIGN KEY (anc_block_id) REFERENCES blocks.ancestral_block (anc_block_id) ON DELETE RESTRICT ;

-- depending on blocks.leaf_block
--~ ALTER TABLE blocks.block_possible_species_node ADD FOREIGN KEY (block_id, block_type) REFERENCES blocks.leaf_block (leaf_block_id, 'leaf') ON DELETE CASCADE ;
ALTER TABLE blocks.gene2block ADD FOREIGN KEY (leaf_block_id) REFERENCES blocks.leaf_block (leaf_block_id) ON DELETE CASCADE ;
ALTER TABLE blocks.functional_homogeneity ADD FOREIGN KEY (leaf_block_id) REFERENCES blocks.leaf_block (leaf_block_id) ON DELETE CASCADE ;


--
-- CREATE VIEWS
--

-- event-related informations

SET search_path = phylogeny, pg_catalog;

-- (reconciled) gene tree 

CREATE VIEW rec_gene_tree AS (
	SELECT gene_tree.family_accession, reconciled_gene_tree.rec_gi_id
		FROM phylogeny.gene_tree
		INNER JOIN phylogeny.reconciled_gene_tree ON reconciled_gene_tree.input_gene_tree_id = gene_tree.gene_tree_id
);


-- event information query/display

CREATE VIEW tree_events AS (
	SELECT gene_tree.family_accession, reconciled_gene_tree.rec_gi_id, reconciled_gene_tree.filepath, event.event_id, event.event_type, event.rec_gi_start_node, event.rec_gi_end_node
		FROM phylogeny.gene_tree
		INNER JOIN phylogeny.reconciled_gene_tree ON reconciled_gene_tree.input_gene_tree_id = gene_tree.gene_tree_id
		INNER JOIN phylogeny.event USING (rec_gi_id)
);

CREATE VIEW event_locs AS (
	SELECT event_id, sp_node_id, characteristic, number 
		FROM event
		INNER JOIN phylogeny.event_possible_species_node USING (event_id)
		INNER JOIN phylogeny.species_node USING (sp_node_id)
);

CREATE OR REPLACE VIEW tree_event_reflocs AS (
	SELECT tree_events.family_accession, tree_events.rec_gi_id, tree_events.event_id, tree_events.event_type, tree_events.rec_gi_start_node, tree_events.rec_gi_end_node, epsnrec.sp_node_id AS rec_sp_node_id, epsndon.sp_node_id AS don_sp_node_id
		FROM phylogeny.tree_events
		INNER JOIN (SELECT * FROM phylogeny.event_possible_species_node INNER JOIN phylogeny.species_node USING (sp_node_id) WHERE characteristic IN ('rec', 'location') AND reference_node=TRUE) AS epsnrec ON tree_events.event_id=epsnrec.event_id
		LEFT JOIN (SELECT * FROM phylogeny.event_possible_species_node INNER JOIN phylogeny.species_node USING (sp_node_id) WHERE characteristic='don' AND reference_node=TRUE) AS epsndon ON tree_events.event_id=epsndon.event_id
);

CREATE OR REPLACE VIEW tree_event_refloc_names AS (
	SELECT tree_events.family_accession, tree_events.rec_gi_id, tree_events.event_id, tree_events.event_type, tree_events.rec_gi_start_node, tree_events.rec_gi_end_node, epsnrec.sp_node_id AS rec_sp_node_id, epsndon.sp_node_id AS don_sp_node_id, epsnrec.number AS rec_spe, epsndon.number AS don_spe
		FROM phylogeny.tree_events
		INNER JOIN (SELECT * FROM phylogeny.event_possible_species_node INNER JOIN phylogeny.species_node USING (sp_node_id) WHERE characteristic IN ('rec', 'location') AND reference_node=TRUE) AS epsnrec ON tree_events.event_id=epsnrec.event_id
		LEFT JOIN (SELECT * FROM phylogeny.event_possible_species_node INNER JOIN phylogeny.species_node USING (sp_node_id) WHERE characteristic='don' AND reference_node=TRUE) AS epsndon ON tree_events.event_id=epsndon.event_id
);

-- event inference confidence

CREATE VIEW node_observing_replicates AS (
	SELECT nr_node_id, unicopy_subtree_id, subtree_name
		FROM phylogeny.reconciled_gene_tree_node
		INNER JOIN phylogeny.node2subtree using (nr_node_id)
		INNER JOIN phylogeny.unicopy_subtree USING (unicopy_subtree_id)
);
	
CREATE VIEW transfer_supporting_replicates AS (
	SELECT event_id, nr_node_id, unicopy_subtree_id, subtree_name
		FROM phylogeny.reconciled_gene_tree_node
		INNER JOIN phylogeny.node2inference USING (nr_node_id)
		INNER JOIN analysis.prunier_inference USING (prinf_id)
		INNER JOIN phylogeny.unicopy_subtree USING (unicopy_subtree_id)
);

CREATE VIEW node_observation AS (
	SELECT nr_node_id, count(unicopy_subtree_id) AS tested
		FROM phylogeny.reconciled_gene_tree_node
		INNER JOIN phylogeny.node2subtree using (nr_node_id)
		INNER JOIN phylogeny.unicopy_subtree USING (unicopy_subtree_id)
	GROUP BY nr_node_id
);

CREATE VIEW transfer_observation AS (
	SELECT event_id, nr_node_id, count(unicopy_subtree_id) AS infered
		FROM phylogeny.reconciled_gene_tree_node
		INNER JOIN phylogeny.node2inference USING (nr_node_id)
		INNER JOIN analysis.prunier_inference USING (prinf_id)
	GROUP BY event_id, nr_node_id
);

CREATE VIEW transfer_inference_confidence AS (
	SELECT nr_node_id, family_accession, node_id, event_id, infered, tested 
		FROM phylogeny.tree_events as te
		INNER JOIN phylogeny.reconciled_gene_tree_node AS rgtn ON te.rec_gi_id=rgtn.rec_gi_id AND te.rec_gi_start_node=rgtn.node_id
		INNER JOIN phylogeny.node_observation USING (nr_node_id)
		INNER JOIN phylogeny.transfer_observation USING (event_id, nr_node_id)
);

CREATE VIEW event_inference_confidence AS (
	SELECT nr_node_id, family_accession, node_id, event_id, infered, tested 
		FROM phylogeny.tree_events as te
		INNER JOIN phylogeny.reconciled_gene_tree_node AS rgtn ON te.rec_gi_id=rgtn.rec_gi_id AND te.rec_gi_start_node=rgtn.node_id
		INNER JOIN phylogeny.node_observation USING (nr_node_id)
		LEFT JOIN phylogeny.transfer_observation USING (event_id, nr_node_id)
);

-- block event-related informations

SET search_path = blocks, pg_catalog;

CREATE VIEW gene_block_events AS (
	SELECT anc_block_id, leaf_block_id, gene_id, event_type
	FROM genome.gene 
		INNER JOIN blocks.gene2block USING (gene_id) 
		INNER JOIN blocks.leaf_block USING (leaf_block_id) 
		INNER JOIN blocks.ancestral_block USING (anc_block_id) 
	ORDER BY anc_block_id
);
	
CREATE VIEW detailed_gene_block_events AS (
	SELECT anc_block_id, leaf_block_id, gene_id, event_type, hogenom_gene_id, family_accession, locus_tag, description
	FROM genome.gene 
		INNER JOIN blocks.gene2block USING (gene_id) 
		INNER JOIN blocks.leaf_block USING (leaf_block_id) 
		INNER JOIN blocks.ancestral_block USING (anc_block_id) 
	ORDER BY anc_block_id
);

CREATE VIEW block_event_reflocs AS (
	SELECT block_id, bpsnrec.sp_node_id AS rec_sp_node_id, bpsndon.sp_node_id AS don_sp_node_id, bpsnrec.number AS rec, bpsndon.number AS don
		FROM (
			SELECT * 
				FROM blocks.block_possible_species_node 
				INNER JOIN phylogeny.species_node USING (sp_node_id)
			WHERE characteristic IN ('rec', 'location') AND reference_node=TRUE
		) AS bpsnrec 
		LEFT JOIN (
			SELECT * 
				FROM blocks.block_possible_species_node 
				INNER JOIN phylogeny.species_node USING (sp_node_id)
			WHERE characteristic='don' AND reference_node=TRUE
		) AS bpsndon USING (block_id)
);


CREATE VIEW gene_block_event_reflocs AS (
	SELECT anc_block_id, leaf_block_id, gene_id, event_type, rec_sp_node_id, don_sp_node_id
		FROM blocks.gene_block_events
		INNER JOIN blocks.block_event_reflocs ON gene_block_events.anc_block_id=block_event_reflocs.block_id
);
--~ CREATE VIEW gene_block_event_reflocs AS (
	--~ SELECT anc_block_id, leaf_block_id, gene_id, event_type, bpsnrec.sp_node_id AS rec_sp_node_id, bpsndon.sp_node_id AS don_sp_node_id
		--~ FROM blocks.gene_block_events
		--~ INNER JOIN (SELECT * FROM blocks.block_possible_species_node WHERE characteristic IN ('rec', 'location') AND reference_node=TRUE) AS bpsnrec ON gene_block_events.anc_block_id=bpsnrec.block_id
		--~ LEFT JOIN (SELECT * FROM blocks.block_possible_species_node WHERE characteristic='don' AND reference_node=TRUE) AS bpsndon ON gene_block_events.anc_block_id=bpsndon.block_id
--~ );
	
	
CREATE VIEW detailed_gene_block_event_reflocs AS (
	SELECT anc_block_id, leaf_block_id, hogenom_gene_id, family_accession, locus_tag, description, event_type, rec, don
		FROM blocks.detailed_gene_block_events
		INNER JOIN blocks.block_event_reflocs ON detailed_gene_block_events.anc_block_id=block_event_reflocs.block_id
);

	
	
