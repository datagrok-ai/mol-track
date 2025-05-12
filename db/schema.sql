CREATE SCHEMA moltrack;

-- Unique chemical structures
CREATE TABLE moltrack.compounds (
    id SERIAL PRIMARY KEY,
    canonical_smiles TEXT NOT NULL,  -- RDKit canonical SMILES
    original_molfile TEXT,           -- as sketched by the chemist
    inchi TEXT NOT NULL,             -- IUPAC InChI
    inchikey TEXT NOT NULL UNIQUE,   -- IUPAC InChIKey
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    is_archived BOOLEAN DEFAULT FALSE
);

-- Batch is a physical sample of the compound (for instance synthesized or procured).
CREATE TABLE moltrack.batches (
    id SERIAL PRIMARY KEY,
    compound_id INTEGER NOT NULL REFERENCES moltrack.compounds(id),
    batch_number TEXT NOT NULL,
    amount REAL,
    amount_unit TEXT,
    purity REAL,
    notes TEXT,
    expiry_date DATE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by INTEGER
);

-- Explains the meaning of a scalar property. 
CREATE TABLE moltrack.semantic_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,    -- e.g., "Molecule", "Cell", "Tissue", "Organism", "Treatment", "Drug", "Image"
    description TEXT
);

-- Properties table - for calculated and measured properties
CREATE TABLE moltrack.properties (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,

    -- value_type defines the colummn in the batch_details, assay_details and assay_results 
    -- tables that store the property value:
    -- * [value_num] for "int" and "double", 
    -- * [value_datetime] for "datetime", 
    -- * [value_uuid] for "uuid", 
    -- * [value_string] for "string"
    value_type TEXT CHECK (value_type IN ('int', 'double', 'bool', 'datetime', 'uuid', 'string')),
    semantic_type_id INTEGER REFERENCES moltrack.semantic_types(id),
    property_class TEXT CHECK (property_class IN ('CALCULATED', 'MEASURED', 'PREDICTED')),
    unit TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- User-defined batch details
create table moltrack.batch_details (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER NOT NULL REFERENCES moltrack.batches(id),
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),

    value_qualifier SMALLINT NOT NULL DEFAULT 0,  -- 0 for "=", 1 for "<", 2 for ">". Only used for numeric properties.
    value_datetime TIMESTAMP WITH TIME ZONE,
    value_uuid uuid,
    value_num REAL,
    value_string TEXT
);

-- Assay type defines what properties are measured in an assay.
-- All of the required properties (see assay_type_properties) must be present in the assay results.
CREATE TABLE moltrack.assay_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT,
    --author_id uuid references moltrack.users(id),
    created_on timestamp without time zone,
    updated_on timestamp without time zone
);

-- Assay type metadata
CREATE TABLE moltrack.assay_type_details (
    assay_type_id INTEGER NOT NULL REFERENCES moltrack.assay_types(id),
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),

    value_datetime TIMESTAMP WITH TIME ZONE,
    value_uuid uuid,
    value_num REAL,
    value_string TEXT,

    PRIMARY KEY (assay_type_id, property_id)
);

-- An experiment where one or more properties were measured for one or more batches.
CREATE TABLE moltrack.assays (
    id SERIAL PRIMARY KEY,
    assay_type_id INTEGER NOT NULL REFERENCES moltrack.assay_types(id),
    name TEXT NOT NULL,      -- e.g., "Kinase Inhibition Assay"
    description TEXT,        -- Detailed description of the assay
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- User-defined assay details 
-- (like when the assay was performed, where, by whom, etc.)
CREATE TABLE moltrack.assay_details (
    assay_id INTEGER NOT NULL REFERENCES moltrack.assays(id),
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),

    value_datetime TIMESTAMP WITH TIME ZONE,
    value_uuid uuid,
    value_num REAL,
    value_string TEXT,

    PRIMARY KEY (assay_id, property_id)
);

-- A set of measurement types for a given assay type.
-- This will be used to validate the data submitted for an assay.
CREATE TABLE moltrack.assay_type_properties (
    assay_type_id INTEGER NOT NULL REFERENCES moltrack.assay_types(id),
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),
    required BOOLEAN NOT NULL DEFAULT FALSE,   -- if true, used for validation of the submitted data
    PRIMARY KEY (assay_type_id, property_id)
);

-- Actual measurements.
CREATE TABLE moltrack.assay_results (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER NOT NULL REFERENCES moltrack.batches(id),
    assay_id INTEGER NOT NULL REFERENCES moltrack.assays(id),  
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),  -- should also be in the corresponding assay_type
    
    -- For performance reasons, only numbers, strings, and booleans are supported for assay results (no datetime or uuid).
    value_qualifier SMALLINT NOT NULL DEFAULT 0,  -- 0 for "=", 1 for "<", 2 for ">". Only used for numeric properties.
    value_num REAL,
    value_string TEXT,
    value_bool BOOLEAN
);

-- Synonym types
CREATE TABLE moltrack.synonym_types (
    id SERIAL PRIMARY KEY,
    synonym_level TEXT NOT NULL,
    name TEXT NOT NULL,
    pattern TEXT NOT NULL,
    description TEXT NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Compound synonyms
CREATE TABLE moltrack.compound_synonyms (
    id SERIAL PRIMARY KEY,
    compound_id INTEGER NOT NULL REFERENCES moltrack.compounds(id),
    synonym_type_id INTEGER NOT NULL REFERENCES moltrack.synonym_types(id),
    synonym_value TEXT NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Batch synonyms
CREATE TABLE moltrack.batch_synonyms (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER NOT NULL moltrack.batches(id),
    synonym_type_id INTEGER NOT NULL REFERENCES moltrack.synonym_types(id),
    synonym_value TEXT NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

GRANT ALL PRIVILEGES ON SCHEMA moltrack TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA moltrack TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA moltrack TO CURRENT_USER;


-- Maintaining RDKit cartrige data in the rdk schema
create extension if not exists rdkit;
drop schema if exists rdk cascade;
create schema rdk;
drop table if exists rdk.mols;

-- rdk.mols contains ids and molecules in RDKit Mol format
select * into rdk.mols 
from (select id, mol_from_smiles(canonical_smiles::cstring) m from moltrack.compounds) tmp;

create index molidx on rdk.mols using gist(m);

-- trigger to insert into rdk.mols after a new compound is inserted into moltrack.compounds
-- Note that we have the "^" symbol in the body below instead of the semicolon because we
-- split the file by semicolon when executing commands (and convert back after splitting)
create or replace function insert_into_rdk_mols()
returns trigger as $$
begin
  insert into rdk.mols(id, m)
  values (
    new.id,
    mol_from_smiles(new.canonical_smiles::cstring)
  )^
  return null^
end^
$$ language plpgsql;

-- trigger to insert into rdk.mols after a new compound is inserted into moltrack.compounds
create trigger compounds_after_insert
after insert on moltrack.compounds
for each row
execute function insert_into_rdk_mols();