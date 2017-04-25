--
-- PostgreSQL database dump
--

SET statement_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = off;
SET check_function_bodies = false;
SET client_min_messages = warning;
SET escape_string_warning = off;

--
-- Name: analysis; Type: SCHEMA; Schema: -; Owner: lassalle
--

CREATE SCHEMA analysis;


ALTER SCHEMA analysis OWNER TO lassalle;

--
-- Name: SCHEMA analysis; Type: COMMENT; Schema: -; Owner: lassalle
--

COMMENT ON SCHEMA analysis IS 'Store information about analysis . TO COMPLETE';


--
-- Name: genome; Type: SCHEMA; Schema: -; Owner: lassalle
--

CREATE SCHEMA genome;


ALTER SCHEMA genome OWNER TO lassalle;

--
-- Name: SCHEMA genome; Type: COMMENT; Schema: -; Owner: lassalle
--

COMMENT ON SCHEMA genome IS 'Store information about genome. TO COMPLETE';


--
-- Name: phylogeny; Type: SCHEMA; Schema: -; Owner: lassalle
--

CREATE SCHEMA phylogeny;


ALTER SCHEMA phylogeny OWNER TO lassalle;

--
-- Name: SCHEMA phylogeny; Type: COMMENT; Schema: -; Owner: lassalle
--

COMMENT ON SCHEMA phylogeny IS 'Store information about phylogeny. TO COMPLETE';


--
-- Name: search; Type: SCHEMA; Schema: -; Owner: lassalle
--

CREATE SCHEMA search;


ALTER SCHEMA search OWNER TO lassalle;

--
-- Name: SCHEMA search; Type: COMMENT; Schema: -; Owner: lassalle
--

COMMENT ON SCHEMA search IS 'Store information that will be requested a lot.';


--
-- Name: taxonomy; Type: SCHEMA; Schema: -; Owner: lassalle
--

CREATE SCHEMA taxonomy;


ALTER SCHEMA taxonomy OWNER TO lassalle;

--
-- Name: SCHEMA taxonomy; Type: COMMENT; Schema: -; Owner: lassalle
--

COMMENT ON SCHEMA taxonomy IS 'Store information about taxonomy. The taxonomy schema is based on the information provided by the NCBI in readme.txt and in ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz TO COMPLETE';


--
-- Name: plpgsql; Type: PROCEDURAL LANGUAGE; Schema: -; Owner: postgres
--

CREATE PROCEDURAL LANGUAGE plpgsql;


ALTER PROCEDURAL LANGUAGE plpgsql OWNER TO postgres;

SET search_path = public, pg_catalog;

--
-- Name: analysis; Type: TYPE; Schema: public; Owner: lassalle
--

CREATE TYPE analysis AS ENUM (
    'REC',
    'G_T',
    'Sp_T',
    'ROOT'
);


ALTER TYPE public.analysis OWNER TO lassalle;

--
-- Name: eventnodecharacteristic; Type: TYPE; Schema: public; Owner: lassalle
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


ALTER TYPE public.eventnodecharacteristic OWNER TO lassalle;

--
-- Name: evolutionaryevent; Type: TYPE; Schema: public; Owner: lassalle
--

CREATE TYPE evolutionaryevent AS ENUM (
    'D',
    'L',
    'T',
    'S',
    'G'
);


ALTER TYPE public.evolutionaryevent OWNER TO lassalle;

--
-- Name: prunier_status; Type: TYPE; Schema: public; Owner: lassalle
--

CREATE TYPE prunier_status AS ENUM (
    'not_done',
    'interupted',
    'done',
    'failed'
);


ALTER TYPE public.prunier_status OWNER TO lassalle;

--
-- Name: analysis; Type: TABLE; Schema: analysis; Owner: lassalle; Tablespace: 
--

SET search_path = analysis, pg_catalog;

CREATE TABLE analysis (
    analysis_id integer NOT NULL,
    parameter_id integer NOT NULL,
    prog_id integer NOT NULL,
    type public.analysis,
    comments text
);


ALTER TABLE analysis.analysis OWNER TO lassalle;

--
-- Name: TABLE analysis; Type: COMMENT; Schema: analysis; Owner: lassalle
--

COMMENT ON TABLE analysis IS 'TO COMPLETE';


--
-- Name: analysis_analysis_id_seq; Type: SEQUENCE; Schema: analysis; Owner: lassalle
--

CREATE SEQUENCE analysis_analysis_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE analysis.analysis_analysis_id_seq OWNER TO lassalle;

--
-- Name: analysis_analysis_id_seq; Type: SEQUENCE OWNED BY; Schema: analysis; Owner: lassalle
--

ALTER SEQUENCE analysis_analysis_id_seq OWNED BY analysis.analysis_id;


--
-- Name: cost; Type: TABLE; Schema: analysis; Owner: lassalle; Tablespace: 
--

CREATE TABLE cost (
    cost_id integer NOT NULL,
    cost_value double precision NOT NULL,
    eventtype public.evolutionaryevent NOT NULL,
    analysis_id integer NOT NULL
);


ALTER TABLE analysis.cost OWNER TO lassalle;

--
-- Name: TABLE cost; Type: COMMENT; Schema: analysis; Owner: lassalle
--

COMMENT ON TABLE cost IS 'TO COMPLETE';


--
-- Name: cost_cost_id_seq; Type: SEQUENCE; Schema: analysis; Owner: lassalle
--

CREATE SEQUENCE cost_cost_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE analysis.cost_cost_id_seq OWNER TO lassalle;

--
-- Name: cost_cost_id_seq; Type: SEQUENCE OWNED BY; Schema: analysis; Owner: lassalle
--

ALTER SEQUENCE cost_cost_id_seq OWNED BY cost.cost_id;


--
-- Name: inference_parameter; Type: TABLE; Schema: analysis; Owner: lassalle; Tablespace: 
--

CREATE TABLE inference_parameter (
    parameter_id integer NOT NULL,
    parameter_string text
);


ALTER TABLE analysis.inference_parameter OWNER TO lassalle;

--
-- Name: TABLE inference_parameter; Type: COMMENT; Schema: analysis; Owner: lassalle
--

COMMENT ON TABLE inference_parameter IS 'TO COMPLETE';


--
-- Name: inference_parameter_parameter_id_seq; Type: SEQUENCE; Schema: analysis; Owner: lassalle
--

CREATE SEQUENCE inference_parameter_parameter_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE analysis.inference_parameter_parameter_id_seq OWNER TO lassalle;

--
-- Name: inference_parameter_parameter_id_seq; Type: SEQUENCE OWNED BY; Schema: analysis; Owner: lassalle
--

ALTER SEQUENCE inference_parameter_parameter_id_seq OWNED BY inference_parameter.parameter_id;


--
-- Name: prunier_inference; Type: TABLE; Schema: analysis; Owner: lassalle; Tablespace: 
--

CREATE TABLE prunier_inference (
    prinf_id integer PRIMARY KEY,
    unicopy_subtree_id integer NOT NULL,
    prunier_round smallint NOT NULL,
    event_id integer,
    analysis_id integer NOT NULL,
    file_path text NOT NULL,
    inverted_rec_don bool
);


ALTER TABLE analysis.prunier_inference OWNER TO lassalle;

--
-- Name: TABLE prunier_inference; Type: COMMENT; Schema: analysis; Owner: lassalle
--

COMMENT ON TABLE prunier_inference IS 'TO COMPLETE';


--
-- Name: inference_parameter_parameter_id_seq; Type: SEQUENCE; Schema: analysis; Owner: lassalle
--

CREATE SEQUENCE prunier_inference_prinf_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE analysis.prunier_inference_prinf_id_seq OWNER TO lassalle;

--
-- Name: prunier_inference_prinf_id_seq; Type: SEQUENCE OWNED BY; Schema: analysis; Owner: lassalle
--

ALTER SEQUENCE prunier_inference_prinf_id_seq OWNED BY prunier_inference.prinf_id;

--
-- Name: prunier_inference_support; Type: TABLE; Schema: analysis; Owner: lassalle; Tablespace: 
--

CREATE TABLE prunier_inference_support (
    prinf_id integer NOT NULL,
    branch_support double precision NOT NULL
);


ALTER TABLE analysis.prunier_inference_support OWNER TO lassalle;

--
-- Name: TABLE prunier_inference_support; Type: COMMENT; Schema: analysis; Owner: lassalle
--

COMMENT ON TABLE prunier_inference_support IS 'TO COMPLETE';



--
-- Name: program; Type: TABLE; Schema: analysis; Owner: lassalle; Tablespace: 
--

CREATE TABLE program (
    prog_id integer NOT NULL,
    name character varying(50) NOT NULL,
    version real NOT NULL
);


ALTER TABLE analysis.program OWNER TO lassalle;

--
-- Name: TABLE program; Type: COMMENT; Schema: analysis; Owner: lassalle
--

COMMENT ON TABLE program IS 'TO COMPLETE';


--
-- Name: program_prog_id_seq; Type: SEQUENCE; Schema: analysis; Owner: lassalle
--

CREATE SEQUENCE program_prog_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE analysis.program_prog_id_seq OWNER TO lassalle;

--
-- Name: program_prog_id_seq; Type: SEQUENCE OWNED BY; Schema: analysis; Owner: lassalle
--

ALTER SEQUENCE program_prog_id_seq OWNED BY program.prog_id;


SET search_path = genome, pg_catalog;

--
-- Name: gene; Type: TABLE; Schema: genome; Owner: lassalle; Tablespace: 
--

CREATE TABLE gene (
    gene_id integer NOT NULL,
    tax_id integer,
    previous_5_gene_id integer,
    following_3_gene_id integer,
    family_accession character varying(50),
    family_type_id smallint,
    name character varying(500),
    chromosome character varying(500),
    genomic_beg integer,
    genomic_end integer,
    strand smallint NOT NULL,
    description text,
    hogenom_gene_id character varying(50) NOT NULL
);


ALTER TABLE genome.gene OWNER TO lassalle;

--
-- Name: TABLE gene; Type: COMMENT; Schema: genome; Owner: lassalle
--

COMMENT ON TABLE gene IS 'TO COMPLETE';


--
-- Name: gene_gene_id_seq; Type: SEQUENCE; Schema: genome; Owner: lassalle
--

CREATE SEQUENCE gene_gene_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE genome.gene_gene_id_seq OWNER TO lassalle;

--
-- Name: gene_gene_id_seq; Type: SEQUENCE OWNED BY; Schema: genome; Owner: lassalle
--

ALTER SEQUENCE gene_gene_id_seq OWNED BY gene.gene_id;


SET search_path = phylogeny, pg_catalog;

--
-- Name: event; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

--~ CREATE TABLE event (
    --~ event_id integer NOT NULL,
    --~ event_type public.evolutionaryevent NOT NULL,
    --~ rec_gi_id integer NOT NULL,
    --~ inf_event_date integer,
    --~ sup_event_date integer,
    --~ support double precision,
    --~ anc_block_id integer NOT NULL,
    --~ rec_gi_start_node character varying(50) NOT NULL,
    --~ rec_gi_end_node character varying(50) NOT NULL
--~ );
CREATE TABLE event (
    event_id integer NOT NULL,
    event_type public.evolutionaryevent NOT NULL,
    rec_gi_id integer NOT NULL,
    inf_event_date integer,
    sup_event_date integer,
    support double precision,
    anc_block_id integer NOT NULL,
    rec_gi_start_node character varying(50) NOT NULL,
    rec_gi_end_node character varying(50)
);


ALTER TABLE phylogeny.event OWNER TO lassalle;

--
-- Name: TABLE event; Type: COMMENT; Schema: phylogeny; Owner: lassalle
--

COMMENT ON TABLE event IS 'TO COMPLETE';


--
-- Name: event_event_id_seq; Type: SEQUENCE; Schema: phylogeny; Owner: lassalle
--

CREATE SEQUENCE event_event_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE phylogeny.event_event_id_seq OWNER TO lassalle;

--
-- Name: event_event_id_seq; Type: SEQUENCE OWNED BY; Schema: phylogeny; Owner: lassalle
--

ALTER SEQUENCE event_event_id_seq OWNED BY event.event_id;

--~ 
--~ --
--~ -- Name: event_per_node; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--~ --
--~ 
--~ CREATE TABLE event_per_node (
    --~ event_id integer NOT NULL,
    --~ sp_node_id integer NOT NULL,
    --~ characteristic public.eventnodecharacteristic
--~ );
--~ 
--~ 
--~ ALTER TABLE phylogeny.event_per_node OWNER TO lassalle;

--
-- Name: event_per_node; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE event_possible_species_node (
    event_id integer NOT NULL,
    sp_node_id integer NOT NULL,
    characteristic public.eventnodecharacteristic,
    reference_node bool NOT NULL
);


ALTER TABLE phylogeny.event_possible_species_node OWNER TO lassalle;

--
-- Name: family; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE family (
    accession character varying(50) NOT NULL,
    family_type_id integer NOT NULL
);


ALTER TABLE phylogeny.family OWNER TO lassalle;

--
-- Name: TABLE family; Type: COMMENT; Schema: phylogeny; Owner: lassalle
--

COMMENT ON TABLE family IS 'TO COMPLETE';


--
-- Name: family_type; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE family_type (
    family_type_id integer NOT NULL,
    name character varying(256) NOT NULL,
    version smallint NOT NULL
);


ALTER TABLE phylogeny.family_type OWNER TO lassalle;

--
-- Name: family_type_family_type_id_seq; Type: SEQUENCE; Schema: phylogeny; Owner: lassalle
--

CREATE SEQUENCE family_type_family_type_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE phylogeny.family_type_family_type_id_seq OWNER TO lassalle;

--
-- Name: family_type_family_type_id_seq; Type: SEQUENCE OWNED BY; Schema: phylogeny; Owner: lassalle
--

ALTER SEQUENCE family_type_family_type_id_seq OWNED BY family_type.family_type_id;


--
-- Name: gene_tree; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE gene_tree (
    gene_tree_id integer NOT NULL,
    file_path text NOT NULL,
    family_accession character varying(50) NOT NULL,
    family_type_id smallint NOT NULL,
    type public.treestatus NOT NULL,
    analysis_id integer,
    source_gene_tree_id integer
);


ALTER TABLE phylogeny.gene_tree OWNER TO lassalle;

--
-- Name: TABLE gene_tree; Type: COMMENT; Schema: phylogeny; Owner: lassalle
--

COMMENT ON TABLE gene_tree IS 'TO COMPLETE';


--
-- Name: gene_tree_gene_tree_id_seq; Type: SEQUENCE; Schema: phylogeny; Owner: lassalle
--

CREATE SEQUENCE gene_tree_gene_tree_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE phylogeny.gene_tree_gene_tree_id_seq OWNER TO lassalle;

--
-- Name: gene_tree_gene_tree_id_seq; Type: SEQUENCE OWNED BY; Schema: phylogeny; Owner: lassalle
--

ALTER SEQUENCE gene_tree_gene_tree_id_seq OWNED BY gene_tree.gene_tree_id;


--
-- Name: reconciled_gene_tree; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
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
    rec_gi_id integer NOT NULL,
    filepath text,
    rec_cost double precision,
    input_gene_tree_id integer NOT NULL,
    species_tree_id integer NOT NULL,
    analysis_id integer NOT NULL,
    tree_num integer
);


ALTER TABLE phylogeny.reconciled_gene_tree OWNER TO lassalle;

--
-- Name: TABLE reconciled_gene_tree; Type: COMMENT; Schema: phylogeny; Owner: lassalle
--

COMMENT ON TABLE reconciled_gene_tree IS 'TO COMPLETE';


--
-- Name: reconciled_gene_tree_rec_gi_id_seq; Type: SEQUENCE; Schema: phylogeny; Owner: lassalle
--

CREATE SEQUENCE reconciled_gene_tree_rec_gi_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE phylogeny.reconciled_gene_tree_rec_gi_id_seq OWNER TO lassalle;

--
-- Name: reconciled_gene_tree_rec_gi_id_seq; Type: SEQUENCE OWNED BY; Schema: phylogeny; Owner: lassalle
--

ALTER SEQUENCE reconciled_gene_tree_rec_gi_id_seq OWNED BY reconciled_gene_tree.rec_gi_id;

--
-- Name: reconciled_gene_tree_node; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE reconciled_gene_tree_node (
    nr_node_id integer PRIMARY KEY,
    node_id integer NOT NULL,
    branch_support double precision,
    branch_length double precision
);


ALTER TABLE phylogeny.reconciled_gene_tree_node OWNER TO lassalle;

--
-- Name: TABLE reconciled_gene_tree_node; Type: COMMENT; Schema: phylogeny; Owner: lassalle
--

COMMENT ON TABLE reconciled_gene_tree_node IS 'TO COMPLETE';


--
-- Name: reconciled_gene_tree_node_nr_node_id_seq; Type: SEQUENCE; Schema: phylogeny; Owner: lassalle
--

CREATE SEQUENCE reconciled_gene_tree_node_nr_node_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE phylogeny.reconciled_gene_tree_node_nr_node_id_seq OWNER TO lassalle;

--
-- Name: reconciled_gene_tree_node_nr_node_id_seq; Type: SEQUENCE OWNED BY; Schema: phylogeny; Owner: lassalle
--

ALTER SEQUENCE reconciled_gene_tree_node_nr_node_id_seq OWNED BY reconciled_gene_tree_node.nr_node_id;

--
-- Name: node2inference; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

--
-- Name: unicopy_subtree; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE unicopy_subtree (
    unicopy_subtree_id integer PRIMARY KEY,
    subtree_name character varying(30) NOT NULL,
    prunier_status public.prunier_status NOT NULL,
    file_path text NOT NULL
);


ALTER TABLE phylogeny.unicopy_subtree OWNER TO lassalle;

--
-- Name: TABLE unicopy_subtree; Type: COMMENT; Schema: phylogeny; Owner: lassalle
--

COMMENT ON TABLE unicopy_subtree IS 'TO COMPLETE';


--
-- Name: unicopy_subtree_unicopy_subtree_id_seq; Type: SEQUENCE; Schema: phylogeny; Owner: lassalle
--

CREATE SEQUENCE unicopy_subtree_unicopy_subtree_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE phylogeny.unicopy_subtree_unicopy_subtree_id_seq OWNER TO lassalle;

--
-- Name: reconciled_gene_tree_node_nr_node_id_seq; Type: SEQUENCE OWNED BY; Schema: phylogeny; Owner: lassalle
--

ALTER SEQUENCE unicopy_subtree_unicopy_subtree_id_seq OWNED BY unicopy_subtree.unicopy_subtree_id;

--
-- Name: node2inference; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE node2inference (
    nr_node_id integer NOT NULL,
    prinf_id integer NOT NULL
);

ALTER TABLE phylogeny.node2inference OWNER TO lassalle;

--
-- Name: node2subtree; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE node2subtree (
    nr_node_id integer NOT NULL,
    unicopy_subtree_id integer NOT NULL
);

ALTER TABLE phylogeny.node2subtree OWNER TO lassalle;

--
-- Name: gene2subtree; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE gene2subtree (
    gene_id integer NOT NULL,
    unicopy_subtree_id integer NOT NULL
);

ALTER TABLE phylogeny.gene2subtree OWNER TO lassalle;

--
-- Name: species_node; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE species_node (
    sp_node_id integer NOT NULL,
    number character varying(50),
    parent_sp_node_id integer,
    species_tree_id integer,
    tax_id integer,
    support double precision,
    uniprot_id character varying(50),
    inf_date integer,
    sup_date integer,
    left_num integer,
    right_num integer,
    branch_length double precision
);


ALTER TABLE phylogeny.species_node OWNER TO lassalle;

--
-- Name: TABLE species_node; Type: COMMENT; Schema: phylogeny; Owner: lassalle
--

COMMENT ON TABLE species_node IS 'TO COMPLETE';


--
-- Name: species_node_sp_node_id_seq; Type: SEQUENCE; Schema: phylogeny; Owner: lassalle
--

CREATE SEQUENCE species_node_sp_node_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE phylogeny.species_node_sp_node_id_seq OWNER TO lassalle;

--
-- Name: species_node_sp_node_id_seq; Type: SEQUENCE OWNED BY; Schema: phylogeny; Owner: lassalle
--

ALTER SEQUENCE species_node_sp_node_id_seq OWNED BY species_node.sp_node_id;


--
-- Name: species_tree; Type: TABLE; Schema: phylogeny; Owner: lassalle; Tablespace: 
--

CREATE TABLE species_tree (
    species_tree_id integer NOT NULL,
    name text NOT NULL,
    file_path text NOT NULL,
    analysis_id integer NOT NULL,
    source_species_tree_id integer
);


ALTER TABLE phylogeny.species_tree OWNER TO lassalle;

--
-- Name: TABLE species_tree; Type: COMMENT; Schema: phylogeny; Owner: lassalle
--

COMMENT ON TABLE species_tree IS 'TO COMPLETE';


--
-- Name: species_tree_species_tree_id_seq; Type: SEQUENCE; Schema: phylogeny; Owner: lassalle
--

CREATE SEQUENCE species_tree_species_tree_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MAXVALUE
    NO MINVALUE
    CACHE 1;


ALTER TABLE phylogeny.species_tree_species_tree_id_seq OWNER TO lassalle;

--
-- Name: species_tree_species_tree_id_seq; Type: SEQUENCE OWNED BY; Schema: phylogeny; Owner: lassalle
--

ALTER SEQUENCE species_tree_species_tree_id_seq OWNED BY species_tree.species_tree_id;


SET search_path = search, pg_catalog;

--
-- Name: lineage; Type: TABLE; Schema: search; Owner: lassalle; Tablespace: 
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


ALTER TABLE search.lineage OWNER TO lassalle;

SET search_path = taxonomy, pg_catalog;

--
-- Name: delnodes; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE delnodes (
    tax_id integer NOT NULL
);


ALTER TABLE taxonomy.delnodes OWNER TO lassalle;

--
-- Name: division; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE division (
    division_id integer NOT NULL,
    division_cde character(3),
    division_name text,
    comments text
);


ALTER TABLE taxonomy.division OWNER TO lassalle;

--
-- Name: gencode; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE gencode (
    genetic_code_id character varying(3) NOT NULL,
    abbreviation text,
    name text,
    cde text,
    starts text
);


ALTER TABLE taxonomy.gencode OWNER TO lassalle;

--
-- Name: merged; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE merged (
    old_tax_id integer NOT NULL,
    new_tax_id integer
);


ALTER TABLE taxonomy.merged OWNER TO lassalle;

--
-- Name: names; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE names (
    tax_id integer,
    name_txt text,
    unique_name text,
    name_class text
);


ALTER TABLE taxonomy.names OWNER TO lassalle;

--
-- Name: nodes; Type: TABLE; Schema: taxonomy; Owner: lassalle; Tablespace: 
--

CREATE TABLE nodes (
    tax_id integer NOT NULL,
    parent_tax_id integer NOT NULL,
    rank character varying(50) NOT NULL,
    embl_code character(2),
    division_id integer,
    inherited_div_flag boolean,
    genetic_code_id character varying(3),
    inherited_gc_flag boolean,
    mitochondrial_genetic_code_id character varying(3),
    inherited_mgc_flag boolean,
    genbank_hidden_flag boolean,
    hidden_subtree_root_flag boolean,
    comments text,
    path public.ltree NOT NULL
);


ALTER TABLE taxonomy.nodes OWNER TO lassalle;

SET search_path = analysis, pg_catalog;

--
-- Name: analysis_id; Type: DEFAULT; Schema: analysis; Owner: lassalle
--

ALTER TABLE analysis ALTER COLUMN analysis_id SET DEFAULT nextval('analysis_analysis_id_seq'::regclass);


--
-- Name: cost_id; Type: DEFAULT; Schema: analysis; Owner: lassalle
--

ALTER TABLE cost ALTER COLUMN cost_id SET DEFAULT nextval('cost_cost_id_seq'::regclass);


--
-- Name: parameter_id; Type: DEFAULT; Schema: analysis; Owner: lassalle
--

ALTER TABLE inference_parameter ALTER COLUMN parameter_id SET DEFAULT nextval('inference_parameter_parameter_id_seq'::regclass);


--
-- Name: prog_id; Type: DEFAULT; Schema: analysis; Owner: lassalle
--

ALTER TABLE program ALTER COLUMN prog_id SET DEFAULT nextval('program_prog_id_seq'::regclass);


--
-- Name: prinf_id; Type: DEFAULT; Schema: analysis; Owner: lassalle
--

ALTER TABLE prunier_inference ALTER COLUMN prinf_id SET DEFAULT nextval('prunier_inference_prinf_id_seq'::regclass);


SET search_path = genome, pg_catalog;

--
-- Name: gene_id; Type: DEFAULT; Schema: genome; Owner: lassalle
--

ALTER TABLE gene ALTER COLUMN gene_id SET DEFAULT nextval('gene_gene_id_seq'::regclass);


SET search_path = phylogeny, pg_catalog;

--
-- Name: event_id; Type: DEFAULT; Schema: phylogeny; Owner: lassalle
--

ALTER TABLE event ALTER COLUMN event_id SET DEFAULT nextval('event_event_id_seq'::regclass);


--
-- Name: family_type_id; Type: DEFAULT; Schema: phylogeny; Owner: lassalle
--

ALTER TABLE family_type ALTER COLUMN family_type_id SET DEFAULT nextval('family_type_family_type_id_seq'::regclass);


--
-- Name: gene_tree_id; Type: DEFAULT; Schema: phylogeny; Owner: lassalle
--

ALTER TABLE gene_tree ALTER COLUMN gene_tree_id SET DEFAULT nextval('gene_tree_gene_tree_id_seq'::regclass);


--
-- Name: rec_gi_id; Type: DEFAULT; Schema: phylogeny; Owner: lassalle
--

ALTER TABLE reconciled_gene_tree ALTER COLUMN rec_gi_id SET DEFAULT nextval('reconciled_gene_tree_rec_gi_id_seq'::regclass);


--
-- Name: nr_node_id; Type: DEFAULT; Schema: phylogeny; Owner: lassalle
--

ALTER TABLE reconciled_gene_tree_node ALTER COLUMN nr_node_id SET DEFAULT nextval('reconciled_gene_tree_node_nr_node_id_seq'::regclass);


--
-- Name: unicopy_subtree_id; Type: DEFAULT; Schema: phylogeny; Owner: lassalle
--

ALTER TABLE unicopy_subtree ALTER COLUMN unicopy_subtree_id SET DEFAULT nextval('unicopy_subtree_unicopy_subtree_id_seq'::regclass);


--
-- Name: sp_node_id; Type: DEFAULT; Schema: phylogeny; Owner: lassalle
--

ALTER TABLE species_node ALTER COLUMN sp_node_id SET DEFAULT nextval('species_node_sp_node_id_seq'::regclass);


--
-- Name: species_tree_id; Type: DEFAULT; Schema: phylogeny; Owner: lassalle
--

ALTER TABLE species_tree ALTER COLUMN species_tree_id SET DEFAULT nextval('species_tree_species_tree_id_seq'::regclass);


SET search_path = analysis, pg_catalog;

