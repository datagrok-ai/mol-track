CREATE SCHEMA moltrack;

-- -- Unique chemical structures
-- CREATE TABLE moltrack.compounds (
--     id SERIAL PRIMARY KEY,
--     canonical_smiles TEXT NOT NULL,  -- RDKit canonical SMILES
--     original_molfile TEXT,           -- as sketched by the chemist
--     inchi TEXT NOT NULL,             -- IUPAC InChI
--     molhash TEXT NOT NULL,
--     formula TEXT NOT NULL,
--     tautomer TEXT NOT NULL,
--     no_stereo_smiles TEXT NOT NULL,
--     no_stereo_tautomer TEXT NOT NULL,
--     sgroup_data TEXT NOT NULL,
--     inchikey TEXT NOT NULL UNIQUE,   -- IUPAC InChIKey
--     created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
--     updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
--     is_archived BOOLEAN DEFAULT FALSE
-- );

-- Table for Compound
CREATE TABLE moltrack.compounds (
    id SERIAL PRIMARY KEY,
    canonical_smiles VARCHAR NULL,
    original_molfile TEXT NULL,
    inchi TEXT NULL,
    inchikey TEXT NULL,
    molhash TEXT NULL,
    formula VARCHAR NULL,
    tautomer VARCHAR NULL,
    no_stereo_smiles VARCHAR NULL,
    no_stereo_tautomer VARCHAR NULL,
    sgroup_data TEXT NULL,
    is_archived BOOLEAN NOT NULL DEFAULT FALSE,
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMP NOT NULL DEFAULT NOW()
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
    value_type TEXT CHECK (value_type IN ('int', 'double', 'bool', 'datetime', 'string')),
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
    result_value REAL
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

-- Assay type properties define what properties are measured in an assay (see assay_type_properties).
CREATE TABLE moltrack.assay_type_properties (
    assay_type_id INTEGER NOT NULL REFERENCES moltrack.assay_types(id),
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),
    required BOOLEAN NOT NULL DEFAULT FALSE,   -- if true, used for validation of the submitted data
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

-- Assay Properties table - for linking assays to properties
CREATE TABLE moltrack.assay_properties (
    assay_id INTEGER NOT NULL REFERENCES moltrack.assays(id),
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),
    PRIMARY KEY (assay_id, property_id)
);

-- Actual measurements.
CREATE TABLE moltrack.assay_results (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER NOT NULL REFERENCES moltrack.batches(id),
    assay_id INTEGER NOT NULL REFERENCES moltrack.assays(id),  
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),  -- should also be in the corresponding assay_type
    result_value REAL
);

GRANT ALL PRIVILEGES ON SCHEMA moltrack TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA moltrack TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA moltrack TO CURRENT_USER;

-- Commenting below for now --

-- -- Maintaining RDKit cartrige data in the rdk schema
-- create extension if not exists rdkit;
-- drop schema if exists rdk cascade;
-- create schema rdk;
-- drop table if exists rdk.mols;

-- -- rdk.mols contains ids and molecules in RDKit Mol format
-- select * into rdk.mols 
-- from (select id, mol_from_smiles(canonical_smiles::cstring) m from moltrack.compounds) tmp;

-- create index molidx on rdk.mols using gist(m);

-- -- trigger to insert into rdk.mols after a new compound is inserted into moltrack.compounds
-- -- Note that we have the "^" symbol in the body below instead of the semicolon because we
-- -- split the file by semicolon when executing commands (and convert back after splitting)
-- create or replace function insert_into_rdk_mols()
-- returns trigger as $$
-- begin
--   insert into rdk.mols(id, m)
--   values (
--     new.id,
--     mol_from_smiles(new.canonical_smiles::cstring)
--   )^
--   return null^
-- end^
-- $$ language plpgsql;

-- -- trigger to insert into rdk.mols after a new compound is inserted into moltrack.compounds
-- create trigger compounds_after_insert
-- after insert on moltrack.compounds
-- for each row
-- execute function insert_into_rdk_mols();