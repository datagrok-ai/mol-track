# MolTrack Client CLI Documentation

The MolTrack Client CLI provides a command-line interface for interacting with the MolTrack API. It supports schema management, compound registration, batch management, properties management, additions management, assays management, and various data operations.

## Installation and Setup

The client is designed to run from the project root directory:

```bash
python mtcli.py [COMMAND] [OPTIONS]
```

## Command Structure

The CLI is organized into several command groups:

- **Schema Commands**: Manage API schema definitions
- **Compound Commands**: Handle compound registration and management
- **Batch Commands**: Manage batch operations and data
- **Additions Commands**: Handle additions management
- **Assays Commands**: Handle assays, assay runs, and assay results management
- **Utility Commands**: General API operations

## Schema Commands

### List Schema

```bash
#List the schema for all entities
python mtcli.py schema list all [OPTIONS]

#List the schema for compounds
python mtcli.py schema list compounds [OPTIONS]

#List the schema for batches
python mtcli.py schema list batches [OPTIONS]

#List the synonyms related to compounds
python mtcli.py schema synonyms compounds [OPTIONS]

#List the synonyms related to batches
python mtcli.py schema synonyms batches [OPTIONS]
```

Lists the current schema for all entities, compounds, or batches.

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

- `--output-format, -o TEXT` Output format: table or json (default: table)
- `--max-rows, -m INTEGER` Maximum number of rows to display in table output (default: None)
- `--output-file, -of TEXT`Path to output file (default: None)

### Load Schema

```bash
python mtcli.py schema load <file_path> [OPTIONS]
```

Adds schema definitions from a JSON file.

**Arguments:**

- `file_path`: Path to the JSON file containing schema data

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Schema File Format:**

```json
{
    "properties": [
        {
            "name": "Molecular Weight",
            "value_type": "double",
            "property_class": "CALCULATED",
            "unit": "g/mol",
            "scope": "COMPOUND",
            "description": "Molecular weight of the compound"
        }
    ],
    "synonym_types": [
        {
            "name": "CAS Number",
            "value_type": "string",
            "property_class": "DECLARED",
            "unit": "",
            "scope": "COMPOUND",
            "pattern": "^\\d{1,7}-\\d{2}-\\d$",
            "description": "CAS Registry Number",
            "semantic_type_id": 1
        }
    ]
}
```

## Compound Commands

### List Compounds

```bash
python mtcli.py compounds list [OPTIONS]
```

Lists compounds using the v1 endpoint.

**Options:**

- `--skip INTEGER`: Number of records to skip (default: 0)
- `--limit INTEGER`: Maximum number of records to return (default: 10)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--output-format, -o TEXT`: Output format: table or json (default: table)
- `--output-file, -of TEXT`: Path to output file (default: None)

### Get Specific Compound

```bash
python mtcli.py compounds get <corporate_compound_id> [OPTIONS]
```
Get a specific compound by corporate_compound_id (friendly name).

**Arguments:**

- `corporate_compound_id`: Corporate Compound ID (friendly name) to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--output-format, -o TEXT`: Output format: table or json (default: table)
- `--output-file, -of TEXT`: Path to output file (default: None)

### Load Compounds from CSV

```bash
python mtcli.py compounds load <csv_file> [OPTIONS]
```

Adds compounds from a CSV file using the `/v1/compounds/` endpoint.

**Arguments:**

- `csv_file`: Path to the CSV file containing compound data

**Options:**

- `--mapping, -m TEXT`: Path to the JSON mapping file (optional)
- `--rows, -r INTEGER`: Number of data rows to process (excludes header row)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--error-handling, -e [reject_all|reject_row]`: Error handling strategy (default: reject_all)
- `--output-format, -o [json|csv]`: Output format (default: json)
- `--dry-run`: Validate data without sending to server
- `--save-errors`: Save error records to a JSON file

**Example Mapping:**

```json
{
    "structure": "smiles",
    "common_name": "compounds_details.common_name",
    "cas": "compounds_details.cas"
}
```

**Usage Examples:**

```bash
# Process all compounds in CSV
python mtcli.py compounds load compounds.csv

# Process only first 5 rows
python mtcli.py compounds load compounds.csv --rows 5

# Use mapping file
python mtcli.py compounds load compounds.csv --mapping mapping.json

# Dry run to validate data
python mtcli.py compounds load compounds.csv --dry-run
```
### Delete specific compound

```bash
python mtcli.py batches delete <corporate_compound_id> [OPTIONS]
```
Delete a specific compound by Corporate Compound ID (friendly name).

**Arguments:**

- `corporate_compound_id`: Corporate Compound ID (friendly name) to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])


## Batch Commands

### List Batches

```bash
python mtcli.py batches list [OPTIONS]
```

Lists batches or gets a specific batch by ID.

**Options:**

- `--skip, -s INTEGER`: Number of records to skip (default: 0)
- `--limit, -l INTEGER`: Maximum number of records to return (default: 10)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--output-format, -o TEXT`: Output format: table or json (default: table)
- `--output-file, -of TEXT`: Path to output file (default: None)

**Usage Examples:**

```bash
# List all batches
python mtcli.py batches list

# List batches with pagination
python mtcli.py batches list --skip 10 --limit 20
```

### Get Specific Batch

```bash
python mtcli.py batches get <corporate_batch_id> [OPTIONS]
```
Get the specific batch by corporate_batch_id (friendly name).

**Arguments:**

- `corporate_batch_id`: Corporate Batch ID (friendly name) to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--output-format, -o TEXT`: Output format: table or json (default: table)
- `--output-file, -of TEXT`: Path to output file (default: None)

### Load Batches from CSV

```bash
python mtcli.py batches load <csv_file> [OPTIONS]
```

Adds batches from a CSV file using the `/v1/batches/` endpoint.

**Arguments:**

- `csv_file`: Path to the CSV file containing batch data

**Options:**

- `--mapping, -m TEXT`: Path to the JSON mapping file (optional)
- `--rows, -r INTEGER`: Number of data rows to process (excludes header row)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--error-handling, -e [reject_all|reject_row]`: Error handling strategy (default: reject_all)
- `--dry-run`: Validate data without sending to server
- `--save-errors`: Save error records to a JSON file

**Example Mapping:**

```json
{
    "structure": "smiles",
    "batch_corporate_id": "batches_details.batch_corporate_id",
    "common_name": "compounds_details.common_name",
    "cas": "compounds_details.cas"
}
```

**Usage Examples:**

```bash
# Process all batches in CSV
python mtcli.py batches load batches.csv

# Process only first 3 rows
python mtcli.py batches load batches.csv --rows 3

# Use mapping file
python mtcli.py batches load batches.csv --mapping mapping.json

# Dry run to validate data
python mtcli.py batches load batches.csv --dry-run
```

### Delete specific batch

```bash
python mtcli.py batches delete <corporate_batch_id> [OPTIONS]
```
Delete a specific batch by Corporate Batch ID (friendly name).

**Arguments:**

- `corporate_batch_id`: Corporate Batch ID (friendly name) to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])


## Additions Commands

### List Additions

```bash
#Load all additions
python mtcli.py additions list all [OPTIONS]

#Load all additions with role of salts
python mtcli.py additions list salts [OPTIONS]

#Load all additions with role of solvates
python mtcli.py additions list solvates [OPTIONS]
```

Lists all additions using the v1 endpoint.

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--output-format, -o TEXT`: Output format: table or json (default: table)
- `--output-file, -of TEXT`: Path to output file (default: None)

### Get specific addition

```bash
python mtcli.py additions get <corporate_batch_id> [OPTIONS]
```
Get the specific addition by addition_id.

**Arguments:**

- `addition_id`:  Addition ID to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--output-format, -o TEXT`: Output format: table or json (default: table)
- `--output-file, -of TEXT`: Path to output file (default: None)

### Load Additions from CSV

```bash
python mtcli.py additions load <csv_file> [OPTIONS]
```

Adds additions from a CSV file using the `/v1/additions/` endpoint.

**Arguments:**

- `csv_file`: Path to the CSV file containing addition data

**Options:**

- `--mapping, -m TEXT`: Path to the JSON mapping file (optional)
- `--rows, -r INTEGER`: Number of data rows to process (excludes header row)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--error-handling, -e [reject_all|reject_row]`: Error handling strategy (default: reject_all)
- `--dry-run`: Validate data without sending to server
- `--save-errors`: Save error records to a JSON file



**Example Mapping:**

```json
{
    "structure": "smiles",
    "common_name": "additions_details.common_name",
    "cas": "additions_details.cas"
}
```

### Update Addition

```bash
python mtcli.py additions update <addition_id> <file_path> [OPTIONS]
```

Updates information for the specified addition.

**Arguments:**

- `addition_id`: Addition ID to update
- `file_path`: Path to the JSON file containing addition update data

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Example Update File:**

```json
{
    "name": "Updated Addition Name",
    "role": "SALT",
    "is_active": true
}
```

### Delete Addition

```bash
python mtcli.py additions delete <addition_id> [OPTIONS]
```

Soft deletes the specified addition (only if no dependent batches exist).

**Arguments:**

- `addition_id`: Addition ID to delete

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Usage Examples:**

```bash
# Process all additions in CSV
python mtcli.py additions load additions.csv

# Process only first 5 rows
python mtcli.py additions load additions.csv --rows 5

# Use mapping file
python mtcli.py additions load additions.csv --mapping mapping.json

# Dry run to validate data
python mtcli.py additions load additions.csv --dry-run

# List all additions
python mtcli.py additions list all

# List all salts
python mtcli.py additions list salts

# List all solvates
python mtcli.py additions list solvates

# Get specific addition
python mtcli.py additions list 123

# Update addition
python mtcli.py additions update 123 update_data.json

# Delete addition
python mtcli.py additions delete 123
```

## Assays Commands

### List Assay data

```bash
#list all assays
python mtcli.py assays list [OPTIONS]

#List all assay runs
python mtcli.py assays runs list [OPTIONS]

#List all assay results
python mtcli.py assays results list [OPTIONS]
```

Lists assays, assays runs or assay results using the v1 endpoint.

**Options:**
- `--skip, -s INTEGER`: Number of records to skip (default: 0)
- `--limit, -l INTEGER`: Maximum number of records to return (default: 10)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--output-format, -o TEXT`: Output format: table or json (default: table)
- `--output-file, -of TEXT`: Path to output file (default: None)

**Usage Examples:**

```bash
# List all assays
python mtcli.py assays list
```

### Load Assay data

```bash
#Load assays
python mtcli.py assays load <file_path> [OPTIONS]

#Load assay runs
python mtcli.py assays load <file_path> [OPTIONS]

#Load assay results
python mtcli.py assays load <file_path> [OPTIONS]
```

Load assay data from a JSON file using the `/v1/assays` endpoint.

**Arguments:**

- `file_path`: Path to the JSON file containing assay data

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Usage Examples:**

```bash
# List all assays
python mtcli.py assays list

# Get specific assay
python mtcli.py assays list 123

# Create assays from JSON
python mtcli.py assays load assays.json

# List all assay runs
python mtcli.py assays runs list

# Get specific assay run
python mtcli.py assays runs list 456

# Create assay runs from CSV
python mtcli.py assays runs load runs.csv --mapping runs_mapping.json

# List all assay results
python mtcli.py assays results list

# Get specific assay result
python mtcli.py assays results list 789

# Create assay results from CSV
python mtcli.py assays results load results.csv --mapping results_mapping.json

# Dry run to validate data
python mtcli.py assays runs load runs.csv --dry-run
```

## Search

```bash

# Search over compounds
python mtcli.py search compounds [OPTIONS]

# Search over batches
python mtcli.py search batches [OPTIONS]

# Search over assays
python mtcli.py search assays [OPTIONS]

# Search over assay runs
python mtcli.py search assay-runs [OPTIONS]

# Search over assay results
python mtcli.py search assay-results [OPTIONS]
```
Advanced search ove compounds, batches, assays, assay-runs, and assay-results.

**Options**
- `--output, -oc TEXT`: Comma-separated list of columns to return
- `--filter, -f TEXT`: Filter as JSON string
- `--aggregations, -a TEXT`: Aggregations as JSON string
- `--output-format, -o TEXT`: Output format: table or json (default: table)
- `--max-rows, -m INTEGER`: Maximum number of rows to display
- `input-file, -if TEXT`: Get serch input from file, if this option is provided options `--output`, `--filter` and `--aggregations` should be ommited
- `--output-file, -of TEXT`: Path to output file (default: None)

**Examples**

```bash
python mtcli.py search assays --filter '{"field": "assays.id","operator": ">","value": 0,"threshold": null}' --output 'assays.id,assays.details.assay format' -o 'csv' -of search_res.csv

python mtcli.py search assays -if search_input.json -of search_res.json
```

**Search input file format**
```json
{
   "output":[
      "compounds.canonical_smiles",
      "compounds.details.chembl",
      "compounds.details.polarSurface"
   ],
   "aggregations":[
      {
         "field": "assay_results.details.ic50",
         "operation": "AVG"
      },
      {
         "field": "assay_results.details.ic50",
         "operation": "COUNT"
      }
   ],
   "filter":{
      "operator":"AND",
      "conditions":[
         {
            "field":"compounds.structure",
            "operator":"IS SIMILAR",
            "value":"Cc1ccc",
            "threshold":0.95
         },
         {
            "operator":"OR",
            "conditions":[
               {
                  "field":"compounds.details.chembl",
                  "operator":"IN",
                  "value":[
                     "CHEMBL123",
                     "CHEMBL123"
                  ]
               },
               {
                  "field":"compounds.details.project",
                  "operator":"=",
                  "value":"My project"
               }
            ]
         },
         {
            "operator":"AND",
            "conditions":[
               {
                  "field":"assay_results.ic50",
                  "operator":">",
                  "value":0.42
               },
               {
                  "field":"assay_results.ec50",
                  "operator":"<",
                  "value":100
               }
            ]
         }
      ]
   }
}
```

## CSV File Format

The client supports CSV files with headers. The first row should contain column names, and subsequent rows should contain data.

**Example CSV Structure:**

```csv
Chemical,CAS #,smiles,corporate_batch_id,corporate_compound_id
"1,3-Benzenedicarboxylic acid",121-91-5,C1=CC(=CC(=C1)C(=O)O)C(=O)O,EPA-001-001,EPA-001
"Celecoxib",169590-42-5,c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N,EPA-120-001,EPA-120
```

## Mapping Files

Mapping files are JSON objects that map CSV column names to API field names.

**Example Mapping:**

```json
{
    "smiles": "smiles",
    "Chemical": "compounds_details.common_name",
    "CAS #": "compounds_details.cas",
    "corporate_batch_id": "batches_details.batch_corporate_id",
    "corporate_compound_id": "compounds_details.corporate_compound_id"
}
```

## Error Handling

The client supports two error handling strategies:

- **reject_all** (default): The entire request fails if any entry is invalid
- **reject_row**: Invalid entries are skipped, while valid ones are processed

## Output Formats

The client supports two output formats:

- **json** (default): Returns structured JSON data
- **csv**: Returns CSV-formatted data

## Common Usage Patterns

### 1. Schema Setup

```bash
# Add schema definitions
python mtcli.py schema load schema.json
```

### 2. Compound Registration

```bash
# Register compounds from CSV
python mtcli.py compounds load compounds.csv --mapping compound_mapping.json
```

### 3. Batch Registration

```bash
# Register batches from CSV
python mtcli.py batches load batches.csv --mapping batch_mapping.json
```

### 4. Additions Registration

```bash
# Register additions from CSV
python mtcli.py additions load additions.csv --mapping additions_mapping.json
```

### 5. Assays Registration

```bash
# Create assays from JSON
python mtcli.py assays load assays.json

# Create assay runs from CSV
python mtcli.py assays runs load runs.csv --mapping runs_mapping.json

# Create assay results from CSV
python mtcli.py assays results load results.csv --mapping results_mapping.json
```

### 6. Data Validation

```bash
# Validate data without sending to server
python mtcli.py compounds load compounds.csv --dry-run
```

### 7. Limited Data Processing

```bash
# Process only first 10 rows for testing
python mtcli.py batches load batches.csv --rows 10
```

## Troubleshooting

### Common Issues

1. **Connection Errors**: Ensure the server is running and accessible at the specified URL
2. **File Not Found**: Verify that CSV and mapping files exist and are accessible
3. **Invalid JSON**: Check that mapping files contain valid JSON
4. **CSV Format**: Ensure CSV files have headers and proper formatting

### Debugging

Use the `--dry-run` option to validate data without sending it to the server:

```bash
python mtcli.py compounds load data.csv --dry-run
```

### Getting Help

Use the `--help` option to get detailed information about any command:

```bash
python mtcli.py --help
python mtcli.py compounds --help
python mtcli.py batches --help
```

## Environment Variables

The client uses the following default configuration:

- **Server URL**: [http://127.0.0.1:8000]
- **Error Handling**: reject_all
- **Output Format**: json

These can be overridden using command-line options.

