# Black Demo Data Set #

The 'black' demo data set is derived from [this paper](https://doi.org/10.1016/j.tox.2021.152819).  It is intended to do some modest testing of the MolTrack capability by adding compounds, batches, assays (types), assays (assay_runs), and assay_results for Hepatocyte stability of 54 compounds tested against multiple species on multiple days.

To boot strap MolTrack, you should

- [ ] Load the schema for compounds using [compounds_schema.json](./compounds_schema.json)
- [ ] Load the schema for batches using [batches_schema.json](./batches_schema.json)
- [ ] Load the schema for assay_types using [assay_types_schema.json](./assay_type_schema.json)
- [ ] Load the schema for assays using [assays_schema.json](./assays_schema.json)

## Compounds ##

To test compounds functionality, you should register [compounds.csv](./compounds.csv)

## Batches ##

To test the batches functionality, you should register [batches.csv](./batches.csv) which should register both the batches and the dependent compounds.

## Assay related ##

To be described.