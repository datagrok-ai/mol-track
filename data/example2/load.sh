#! /bin/bash

working_dir=./data/example2

# clean up any existing data
python -m client.client database clean -f 

# # load additions
# python -m client.client additions load ${working_dir}/additions.csv

# load compound schema
python -m client.client schema load ${working_dir}/compounds_schema.json

# load batches schema
python -m client.client schema load ${working_dir}/batches_schema.json

# load compounds
python -m client.client compounds load ${working_dir}/compounds.csv -e reject_row --save-errors

# load batches
python -m client.client batches load ${working_dir}/batches.csv -e reject_row --save-errors --mapping ${working_dir}/batches_mapping.json

# load assay schema
python -m client.client schema load ${working_dir}/assays_schema.json

# load assay runs schema
python -m client.client schema load ${working_dir}/assay_runs_schema.json

# load assay results schema
python -m client.client schema load ${working_dir}/assay_results_schema.json

# load assays
python -m client.client assays load ${working_dir}/assays.json

# load assay runs
python -m client.client assays runs load ${working_dir}/assay_runs.csv -e reject_row --save-errors -m ${working_dir}/assay_runs_mapping.json

# # load assay results
python -m client.client assays results load ${working_dir}/assay_results.csv -e reject_row --save-errors --mapping ${working_dir}/assay_results_mapping.json

# database stats
python -m client.client database stats