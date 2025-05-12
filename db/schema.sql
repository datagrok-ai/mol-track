-- SQL dump generated using DBML (dbml.dbdiagram.io)
-- Database: PostgreSQL
-- Generated at: 2025-05-12T12:19:53.665Z

CREATE SCHEMA "moltrack";

CREATE TABLE "moltrack"."users" (
  "id" uuid PRIMARY KEY NOT NULL,
  "email" text UNIQUE NOT NULL,
  "first_name" text NOT NULL,
  "last_name" text NOT NULL,
  "has_password" boolean NOT NULL,
  "is_active" boolean NOT NULL DEFAULT false,
  "is_service_account" boolean NOT NULL DEFAULT false,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL
);

CREATE TABLE "moltrack"."settings" (
  "id" serial PRIMARY KEY,
  "name" text NOT NULL,
  "value" text NOT NULL,
  "description" text NOT NULL
);

CREATE TABLE "moltrack"."compounds" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "canonical_smiles" text NOT NULL,
  "original_molfile" text,
  "molregno" int UNIQUE NOT NULL,
  "inchi" text NOT NULL,
  "inchikey" text UNIQUE NOT NULL,
  "formula" text NOT NULL,
  "hash_mol" uuid UNIQUE NOT NULL,
  "hash_tautomer" uuid UNIQUE NOT NULL,
  "hash_canonical_smiles" uuid UNIQUE NOT NULL,
  "hash_no_stereo_smiles" uuid UNIQUE NOT NULL,
  "hash_no_stereo_tautomer" uuid UNIQUE NOT NULL,
  "is_archived" boolean DEFAULT false,
  "deleted_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "deleted_by" uuid NOT NULL
);

CREATE TABLE "moltrack"."batches" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "compound_id" int NOT NULL,
  "batch_number" text NOT NULL,
  "amount" real,
  "amount_unit" text,
  "purity" real,
  "notes" text,
  "expiry_date" date
);

CREATE TABLE "moltrack"."additions" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "name" text NOT NULL,
  "description" text,
  "display_name" text,
  "display_order" int,
  "is_active" boolean DEFAULT true,
  "formula" text,
  "molecular_weight" real,
  "smiles" text,
  "molfile" text,
  "role" text,
  "is_archived" boolean DEFAULT false,
  "deleted_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "deleted_by" uuid NOT NULL
);

CREATE TABLE "moltrack"."batch_additions" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "batch_id" int NOT NULL,
  "addition_id" int NOT NULL,
  "addition_equivalent" real NOT NULL DEFAULT 1,
  PRIMARY KEY ("batch_id", "addition_id")
);

CREATE TABLE "moltrack"."semantic_types" (
  "id" serial PRIMARY KEY,
  "name" text NOT NULL,
  "description" text
);

CREATE TABLE "moltrack"."properties" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "name" text NOT NULL,
  "value_type" text,
  "semantic_type_id" int,
  "property_class" text,
  "unit" text
);

CREATE TABLE "moltrack"."batch_details" (
  "id" serial PRIMARY KEY,
  "batch_id" int NOT NULL,
  "property_id" int NOT NULL,
  "value_datetime" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "value_uuid" uuid,
  "value_num" float,
  "value_string" text
);

CREATE TABLE "moltrack"."synonym_types" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "synonym_level" text NOT NULL,
  "name" text NOT NULL,
  "pattern" text NOT NULL,
  "description" text NOT NULL
);

CREATE TABLE "moltrack"."compound_synonyms" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "compound_id" int NOT NULL,
  "synonym_type_id" int NOT NULL,
  "synonym_value" text NOT NULL
);

CREATE TABLE "moltrack"."batch_synonyms" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "batch_id" int NOT NULL,
  "synonym_type_id" int NOT NULL,
  "synonym_value" text NOT NULL
);

CREATE TABLE "moltrack"."compound_details" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "compound_id" int NOT NULL,
  "property_id" int NOT NULL,
  "value_datetime" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "value_uuid" uuid,
  "value_num" float,
  "value_string" text
);

CREATE TABLE "moltrack"."assay_types" (
  "id" serial PRIMARY KEY,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL,
  "name" text NOT NULL,
  "description" text
);

CREATE TABLE "moltrack"."assay_type_details" (
  "assay_type_id" int NOT NULL,
  "property_id" int NOT NULL,
  "required" boolean NOT NULL DEFAULT false,
  "value_datetime" timestamp,
  "value_uuid" uuid,
  "value_num" float,
  "value_string" text
);

CREATE TABLE "moltrack"."assays" (
  "id" serial PRIMARY KEY,
  "assay_type_id" int NOT NULL,
  "name" text NOT NULL,
  "description" text,
  "created_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "updated_at" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "created_by" uuid NOT NULL,
  "updated_by" uuid NOT NULL
);

CREATE TABLE "moltrack"."assay_details" (
  "assay_id" int NOT NULL,
  "property_id" int NOT NULL,
  "value_datetime" timestamp DEFAULT (CURRENT_TIMESTAMP),
  "value_uuid" uuid,
  "value_num" float,
  "value_string" text
);

CREATE TABLE "moltrack"."assay_type_properties" (
  "assay_type_id" int NOT NULL,
  "property_id" int NOT NULL,
  "required" bool NOT NULL DEFAULT false,
  PRIMARY KEY ("assay_type_id", "property_id")
);

CREATE TABLE "moltrack"."assay_results" (
  "id" serial PRIMARY KEY,
  "batch_id" int NOT NULL,
  "assay_id" int NOT NULL,
  "property_id" int NOT NULL,
  "value_qualifier" smallint NOT NULL DEFAULT 0,
  "value_num" float,
  "value_string" text,
  "value_bool" boolean
);

CREATE UNIQUE INDEX ON "moltrack"."users" ("email");

CREATE INDEX ON "moltrack"."compounds" ("canonical_smiles");

CREATE INDEX ON "moltrack"."compounds" ("inchi");

CREATE INDEX ON "moltrack"."compounds" ("inchikey");

CREATE INDEX ON "moltrack"."compounds" ("formula");

CREATE UNIQUE INDEX ON "moltrack"."compounds" ("hash_mol");

CREATE UNIQUE INDEX ON "moltrack"."compounds" ("hash_tautomer");

CREATE UNIQUE INDEX ON "moltrack"."compounds" ("hash_canonical_smiles");

CREATE UNIQUE INDEX ON "moltrack"."compounds" ("hash_no_stereo_smiles");

CREATE UNIQUE INDEX ON "moltrack"."compounds" ("hash_no_stereo_tautomer");

CREATE UNIQUE INDEX ON "moltrack"."compounds" ("molregno");

CREATE INDEX ON "moltrack"."batches" ("compound_id");

CREATE INDEX ON "moltrack"."batches" ("batch_number");

CREATE INDEX ON "moltrack"."additions" ("name");

CREATE INDEX ON "moltrack"."additions" ("display_name");

CREATE INDEX ON "moltrack"."additions" ("display_order");

CREATE INDEX ON "moltrack"."additions" ("formula");

CREATE INDEX ON "moltrack"."additions" ("role");

CREATE INDEX ON "moltrack"."properties" ("name");

CREATE INDEX ON "moltrack"."properties" ("value_type");

CREATE INDEX ON "moltrack"."properties" ("semantic_type_id");

CREATE INDEX ON "moltrack"."properties" ("property_class");

CREATE INDEX ON "moltrack"."batch_details" ("batch_id");

CREATE INDEX ON "moltrack"."batch_details" ("property_id");

CREATE UNIQUE INDEX ON "moltrack"."synonym_types" ("synonym_level", "name");

CREATE INDEX ON "moltrack"."compound_synonyms" ("compound_id");

CREATE INDEX ON "moltrack"."compound_synonyms" ("synonym_type_id");

CREATE INDEX ON "moltrack"."batch_synonyms" ("batch_id");

CREATE INDEX ON "moltrack"."batch_synonyms" ("synonym_type_id");

CREATE INDEX ON "moltrack"."compound_details" ("compound_id");

CREATE INDEX ON "moltrack"."compound_details" ("property_id");

CREATE UNIQUE INDEX ON "moltrack"."assay_type_details" ("assay_type_id", "property_id");

CREATE UNIQUE INDEX ON "moltrack"."assay_details" ("assay_id", "property_id");

CREATE INDEX ON "moltrack"."assay_results" ("batch_id");

CREATE INDEX ON "moltrack"."assay_results" ("assay_id");

CREATE INDEX ON "moltrack"."assay_results" ("property_id");

COMMENT ON TABLE "moltrack"."users" IS 'for MVP, just a single global user will be configured
';

COMMENT ON TABLE "moltrack"."settings" IS 'This table is used to store configuration settings for the application.
examples include:
- compound standardization rules
- compound uniqueness rules
- compound identification rules and synonym generation rules
- batch identification rules and synonym generation rules
I suspect that there should be some more structure to this table, but for MVP, this will be used to store the rules for the global user.
';

COMMENT ON COLUMN "moltrack"."compounds"."canonical_smiles" IS 'RDKit canonical SMILES';

COMMENT ON COLUMN "moltrack"."compounds"."original_molfile" IS 'as sketched by the chemist';

COMMENT ON COLUMN "moltrack"."compounds"."molregno" IS 'generated by the system based on business rules';

COMMENT ON COLUMN "moltrack"."compounds"."inchi" IS 'IUPAC InChI';

COMMENT ON COLUMN "moltrack"."compounds"."inchikey" IS 'IUPAC InChIKey';

COMMENT ON TABLE "moltrack"."batches" IS 'need to refactor to identify what should be in details and synonyms
I believe that amount, amount_unit, purity shoujld be properties
';

COMMENT ON TABLE "moltrack"."additions" IS 'Questions:
- how are name, display_name, and descrition used
- should we have UI related fields in this table, e.g., display_name, display_order
- should we do check constraints for role?
';

COMMENT ON COLUMN "moltrack"."additions"."role" IS '(''SALT'', ''SOLVATE'', ''REAGENT'', ''SOLVENT'', ''OTHER'')';

COMMENT ON TABLE "moltrack"."batch_additions" IS 'This allows us to have multiple salt and solvates for a given batch''
';

COMMENT ON COLUMN "moltrack"."batch_additions"."addition_equivalent" IS 'equivalent relative to parent batch';

COMMENT ON COLUMN "moltrack"."semantic_types"."name" IS 'e.g., "Molecule", "Cell", ...';

COMMENT ON TABLE "moltrack"."properties" IS 'for validity checking should we have field to specify what entity this should be relevant
future: property should have ontology reference
';

COMMENT ON COLUMN "moltrack"."properties"."value_type" IS '(''int'', ''double'', ''bool'', ''datetime'', ''string'')';

COMMENT ON COLUMN "moltrack"."properties"."property_class" IS '(''CALCULATED'', ''MEASURED'', ''PREDICTED'')';

COMMENT ON COLUMN "moltrack"."properties"."unit" IS 'future: unit should point to unit ontology';

COMMENT ON TABLE "moltrack"."batch_details" IS 'future: batch_details.value_string should have optional ontology reference';

COMMENT ON COLUMN "moltrack"."batch_details"."value_datetime" IS 'with time zone';

COMMENT ON TABLE "moltrack"."synonym_types" IS 'how to address combined uniqueness constraint';

COMMENT ON COLUMN "moltrack"."synonym_types"."synonym_level" IS '(''BATCH'',''COMPOUND'')';

COMMENT ON COLUMN "moltrack"."synonym_types"."name" IS 'e.g., "CAS", "USAN", "INN", "tradename"';

COMMENT ON COLUMN "moltrack"."synonym_types"."pattern" IS 'regex for identifier: CHEMBL.*';

COMMENT ON TABLE "moltrack"."compound_details" IS 'future: compound_details.value_string should have optional ontology reference';

COMMENT ON COLUMN "moltrack"."compound_details"."value_datetime" IS 'with time zone';

COMMENT ON TABLE "moltrack"."assay_type_details" IS 'PK: (assay_type_id, property_id)';

COMMENT ON COLUMN "moltrack"."assay_type_details"."value_datetime" IS 'timestamp with time zone';

COMMENT ON COLUMN "moltrack"."assays"."name" IS 'e.g., "Kinase Inhibition Assay"';

COMMENT ON TABLE "moltrack"."assay_details" IS 'PK: (assay_id, property_id)';

COMMENT ON COLUMN "moltrack"."assay_details"."value_datetime" IS 'with time zone';

COMMENT ON TABLE "moltrack"."assay_type_properties" IS 'A set of measurement types for a given assay type.
This will be used to validate the data submitted for an assay.
';

COMMENT ON COLUMN "moltrack"."assay_type_properties"."required" IS 'if true, must exist in submitted data for validation';

COMMENT ON TABLE "moltrack"."assay_results" IS 'For performance reasons, only numbers, strings, and booleans are supported for assay results (no datetime or uuid).';

COMMENT ON COLUMN "moltrack"."assay_results"."value_qualifier" IS 'SMALLINT, 0 for "=", 1 for "<", 2 for ">". Only used for numeric properties.';

ALTER TABLE "moltrack"."batches" ADD FOREIGN KEY ("compound_id") REFERENCES "moltrack"."compounds" ("id");

ALTER TABLE "moltrack"."batch_additions" ADD FOREIGN KEY ("batch_id") REFERENCES "moltrack"."batches" ("id");

ALTER TABLE "moltrack"."batch_additions" ADD FOREIGN KEY ("addition_id") REFERENCES "moltrack"."additions" ("id");

ALTER TABLE "moltrack"."properties" ADD FOREIGN KEY ("semantic_type_id") REFERENCES "moltrack"."semantic_types" ("id");

ALTER TABLE "moltrack"."batch_details" ADD FOREIGN KEY ("batch_id") REFERENCES "moltrack"."batches" ("id");

ALTER TABLE "moltrack"."batch_details" ADD FOREIGN KEY ("property_id") REFERENCES "moltrack"."properties" ("id");

ALTER TABLE "moltrack"."compound_synonyms" ADD FOREIGN KEY ("compound_id") REFERENCES "moltrack"."compounds" ("id");

ALTER TABLE "moltrack"."compound_synonyms" ADD FOREIGN KEY ("synonym_type_id") REFERENCES "moltrack"."synonym_types" ("id");

ALTER TABLE "moltrack"."batch_synonyms" ADD FOREIGN KEY ("batch_id") REFERENCES "moltrack"."batches" ("id");

ALTER TABLE "moltrack"."batch_synonyms" ADD FOREIGN KEY ("synonym_type_id") REFERENCES "moltrack"."synonym_types" ("id");

ALTER TABLE "moltrack"."compound_details" ADD FOREIGN KEY ("compound_id") REFERENCES "moltrack"."compounds" ("id");

ALTER TABLE "moltrack"."compound_details" ADD FOREIGN KEY ("property_id") REFERENCES "moltrack"."properties" ("id");

ALTER TABLE "moltrack"."assay_type_details" ADD FOREIGN KEY ("assay_type_id") REFERENCES "moltrack"."assay_types" ("id");

ALTER TABLE "moltrack"."assay_type_details" ADD FOREIGN KEY ("property_id") REFERENCES "moltrack"."properties" ("id");

ALTER TABLE "moltrack"."assays" ADD FOREIGN KEY ("assay_type_id") REFERENCES "moltrack"."assay_types" ("id");

ALTER TABLE "moltrack"."assay_details" ADD FOREIGN KEY ("assay_id") REFERENCES "moltrack"."assays" ("id");

ALTER TABLE "moltrack"."assay_details" ADD FOREIGN KEY ("property_id") REFERENCES "moltrack"."properties" ("id");

ALTER TABLE "moltrack"."assay_type_properties" ADD FOREIGN KEY ("assay_type_id") REFERENCES "moltrack"."assay_types" ("id");

ALTER TABLE "moltrack"."assay_type_properties" ADD FOREIGN KEY ("property_id") REFERENCES "moltrack"."properties" ("id");

ALTER TABLE "moltrack"."assay_results" ADD FOREIGN KEY ("batch_id") REFERENCES "moltrack"."batches" ("id");

ALTER TABLE "moltrack"."assay_results" ADD FOREIGN KEY ("assay_id") REFERENCES "moltrack"."assays" ("id");

ALTER TABLE "moltrack"."assay_results" ADD FOREIGN KEY ("property_id") REFERENCES "moltrack"."properties" ("id");

ALTER TABLE "moltrack"."users" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."users" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."compounds" ADD FOREIGN KEY ("deleted_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."compounds" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."compounds" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."batches" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."batches" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."additions" ADD FOREIGN KEY ("deleted_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."additions" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."additions" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."batch_additions" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."batch_additions" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."properties" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."properties" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."synonym_types" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."synonym_types" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."compound_synonyms" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."compound_synonyms" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."batch_synonyms" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."batch_synonyms" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."compound_details" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."compound_details" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."assay_types" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."assay_types" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."assays" ADD FOREIGN KEY ("created_by") REFERENCES "moltrack"."users" ("id");

ALTER TABLE "moltrack"."assays" ADD FOREIGN KEY ("updated_by") REFERENCES "moltrack"."users" ("id");
