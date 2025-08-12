# Example 2 Dataset #

This is a synthetic dataset with a reasonably coherent set of molecule and assay data at both single point and curve level and multiple assays

## Compounds ##

- [schema](compounds_schema.json)
- [mapping](compounds_mapping.json)
- [data](compounds.csv)

## Batches ##

- [schema](batches_schema.json)
- [mapping](batches_mapping.json)
- [data](batches.csv)

## Assay Data ##

- Assays
  - [schema](assay_types_schema.json)
  - [data](./assays.json)
- Assay Runs
  - [schema](assay_runs_schema.json)
  - [mapping](assay_runs_mapping.json)
  - [data](assay_runs.csv)
- Assay Results
  - [schema](./assay_results_schema.json)
  - [mapping](./assay_results_mapping.json)
  - [data](./assay_results.csv)

[assay_results.csv](assay_results.csv) is the original data file.  I created [cols.csv](./cols.csv) to characterize which columns should be mapped to which entity_type in the assay domain. 

[load.sh](load.csv) is used in conjunction with the client/client.py to perform the loading.  Mapping files are required for the loading to accommodate the original column names.
