CREATE SCHEMA moltrack;

-- Compounds table - core table for chemical structures
CREATE TABLE moltrack.compounds (
    id SERIAL PRIMARY KEY,
    canonical_smiles TEXT NOT NULL,
    inchi TEXT NOT NULL,
    inchikey TEXT NOT NULL UNIQUE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    is_archived BOOLEAN DEFAULT FALSE
);

-- Batches table - for tracking physical samples
CREATE TABLE moltrack.batches (
    id SERIAL PRIMARY KEY,
    compound_id INTEGER NOT NULL REFERENCES moltrack.compounds(id),
    batch_number TEXT NOT NULL,
    amount REAL,
    amount_unit TEXT,
    purity REAL,
    vendor TEXT,
    catalog_id TEXT,
    acquisition_date DATE,
    expiry_date DATE,
    storage_location TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by INTEGER,
    notes TEXT
);

-- Properties table - for calculated and measured properties
CREATE TABLE moltrack.properties (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    value_type TEXT CHECK (value_type IN ('NUMBER', 'STRING', 'BOOL')),
    property_class TEXT CHECK (property_class IN ('CALCULATED', 'MEASURED', 'PREDICTED')),
    unit TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Assays table - for assay types
CREATE TABLE moltrack.assays (
    id SERIAL PRIMARY KEY,
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

-- Assay Results table - for assay results
CREATE TABLE moltrack.assay_results (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER NOT NULL REFERENCES moltrack.batches(id),
    assay_id INTEGER NOT NULL REFERENCES moltrack.assays(id),
    property_id INTEGER NOT NULL REFERENCES moltrack.properties(id),
    result_value REAL
);

GRANT ALL PRIVILEGES ON SCHEMA moltrack TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA moltrack TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA moltrack TO CURRENT_USER;