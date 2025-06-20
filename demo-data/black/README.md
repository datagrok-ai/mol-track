# Black Demo Data Set #

The 'black' demo data set is derived from [this paper](https://doi.org/10.1016/j.tox.2021.152819).  It is intended to for modest testing of the MolTrack capability by adding compounds, batches, assays (types), assays (assay_runs), and assay_results for Hepatocyte stability of 54 compounds tested against multiple species on multiple days.

To boot strap MolTrack, you should

- [ ] Load the schema for compounds using [compounds_schema.json](./compounds_schema.json)
- [ ] Load the schema for batches using [batches_schema.json](./batches_schema.json)
- [ ] Load the schema for assay data domain using [assay_data_schema.json](assay_data_schema.json)

## Compounds ##

To test compounds functionality, you should register [compounds.csv](./compounds.csv).  A mapping file, [compounds_mapping](compounds_mapping.json) is provided for your convenience.

## Batches ##

To test the batches functionality, you should register [batches.csv](./batches.csv) which should register both the batches and the dependent compounds.  A mapping file [batches_mapping](./batches_mapping.json) is provided for your convenience.

## Assay related ##

Assays and the corresponding results profiles (assay_properties) will be created using [assays_instances.json](assays_instances.json)

Assays Runs will be created using [assay_runs.csv](assay_runs.csv).  The mapping of csv columns to fixed and dynamic assay_run attributes will be done by case insensitve mapping or a mapping file can be provided. A potential mapping file is [assay_runs_mapping.json](assay_runs_mapping.json)

Assay results will be shown in [assay_results.csv](assay_results.csv).  Processing this file you should be creating assays instances with *assay_details* property values and a *assay_results* instances with property/values that correspond to those declared in the *assay_type_properties* table. A potential mapping file is [assay_results_mapping.json](./assay_results_mapping.json).  It is critical that one column of the csv can be used to uniquely lookup the batch.

Assay raw data will be shown in assay_raw_results.csv.  This will be used for assay results calculation, not MolTrack registration.