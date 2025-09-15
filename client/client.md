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
- **Properties Commands**: Handle properties management
- **Additions Commands**: Handle additions management
- **Assays Commands**: Handle assays, assay runs, and assay results management
- **Utility Commands**: General API operations

## Schema Commands

### List Schema

```bash
python mtcli.py schema list [OPTIONS]
```

Lists the current schema for both compounds and batches.

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

### Create Schema

```bash
python mtcli.py schema create <file_path> [OPTIONS]
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
python mtcli.py compound list [OPTIONS]
```

Lists compounds using the v1 endpoint.

**Options:**

- `--skip INTEGER`: Number of records to skip (default: 0)
- `--limit INTEGER`: Maximum number of records to return (default: 10)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

### Create Compound (Legacy)

```bash
python mtcli.py compound create_old <item_json> [OPTIONS]
```

Creates a compound using the legacy `/compounds` endpoint.

**Arguments:**

- `item_json`: JSON string containing compound data

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

### Create Compounds from CSV

```bash
python mtcli.py compound create <csv_file> [OPTIONS]
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
python mtcli.py compound create compounds.csv

# Process only first 5 rows
python mtcli.py compound create compounds.csv --rows 5

# Use mapping file
python mtcli.py compound create compounds.csv --mapping mapping.json

# Dry run to validate data
python mtcli.py compound create compounds.csv --dry-run
```

## Batch Commands

### List Batches

```bash
python mtcli.py batch list [batch_id] [OPTIONS]
```

Lists batches or gets a specific batch by ID.

**Arguments:**

- `batch_id`: Optional batch ID to retrieve

**Options:**

- `--skip, -s INTEGER`: Number of records to skip (default: 0)
- `--limit, -l INTEGER`: Maximum number of records to return (default: 10)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Usage Examples:**

```bash
# List all batches
python mtcli.py batch list

# List batches with pagination
python mtcli.py batch list --skip 10 --limit 20

# Get specific batch
python mtcli.py batch list 123
```

### Create Batches from CSV

```bash
python mtcli.py batch create <csv_file> [OPTIONS]
```

Adds batches from a CSV file using the `/v1/batches/` endpoint.

**Arguments:**

- `csv_file`: Path to the CSV file containing batch data

**Options:**

- `--mapping, -m TEXT`: Path to the JSON mapping file (optional)
- `--rows, -r INTEGER`: Number of data rows to process (excludes header row)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--error-handling, -e [reject_all|reject_row]`: Error handling strategy (default: reject_all)
- `--output-format, -o [json|csv]`: Output format (default: json)
- `--dry-run`: Validate data without sending to server

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
python mtcli.py batch create batches.csv

# Process only first 3 rows
python mtcli.py batch create batches.csv --rows 3

# Use mapping file
python mtcli.py batch create batches.csv --mapping mapping.json

# Dry run to validate data
python mtcli.py batch create batches.csv --dry-run
```

### Batch List Subcommands

#### Get Batch Properties

```bash
python mtcli.py batch list properties <batch_id> [OPTIONS]
```

Gets all properties for a specific batch.

**Arguments:**

- `batch_id`: Batch ID to get properties for

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

#### Get Batch Synonyms

```bash
python mtcli.py batch list synonyms <batch_id> [OPTIONS]
```

Gets all synonyms for a specific batch.

**Arguments:**

- `batch_id`: Batch ID to get synonyms for

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

#### Get Batch Additions

```bash
python mtcli.py batch list additions <batch_id> [OPTIONS]
```

Gets all additions for a specific batch.

**Arguments:**

- `batch_id`: Batch ID to get additions for

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

## Properties Commands

### List Properties

```bash
python mtcli.py properties list [OPTIONS]
```

Lists all properties.

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

## Additions Commands

### List Additions

```bash
python mtcli.py additions list [OPTIONS]
```

Lists all additions using the v1 endpoint.

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

### Create Additions from CSV

```bash
python mtcli.py additions create <csv_file> [OPTIONS]
```

Adds additions from a CSV file using the `/v1/additions/` endpoint.

**Arguments:**

- `csv_file`: Path to the CSV file containing addition data

**Options:**

- `--mapping, -m TEXT`: Path to the JSON mapping file (optional)
- `--rows, -r INTEGER`: Number of data rows to process (excludes header row)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--error-handling, -e [reject_all|reject_row]`: Error handling strategy (default: reject_all)
- `--output-format, -o [json|csv]`: Output format (default: json)
- `--dry-run`: Validate data without sending to server

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

### Additions List Subcommands

#### List Additions Salts

```bash
python mtcli.py additions list salts [OPTIONS]
```

Lists all additions with role of salts.

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

#### List Additions Solvates

```bash
python mtcli.py additions list solvates [OPTIONS]
```

Lists all additions with role of solvates.

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

#### Get Addition

```bash
python mtcli.py additions list <addition_id> [OPTIONS]
```

Gets all information for a specific addition.

**Arguments:**

- `addition_id`: Addition ID to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Usage Examples:**

```bash
# Process all additions in CSV
python mtcli.py additions create additions.csv

# Process only first 5 rows
python mtcli.py additions create additions.csv --rows 5

# Use mapping file
python mtcli.py additions create additions.csv --mapping mapping.json

# Dry run to validate data
python mtcli.py additions create additions.csv --dry-run

# List all additions
python mtcli.py additions list

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

### List Assays

```bash
python mtcli.py assays list [assay_id] [OPTIONS]
```

Lists assays using the v1 endpoint.

**Arguments:**

- `assay_id`: Optional assay ID to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Usage Examples:**

```bash
# List all assays
python mtcli.py assays list

# Get specific assay
python mtcli.py assays list 123
```

### Create Assays

```bash
python mtcli.py assays create <file_path> [OPTIONS]
```

Creates assays from a JSON file using the `/v1/assays` endpoint.

**Arguments:**

- `file_path`: Path to the JSON file containing assay data

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

### Assay Runs Commands

#### List Assay Runs

```bash
python mtcli.py assays runs list [assay_run_id] [OPTIONS]
```

Lists assay runs using the v1 endpoint.

**Arguments:**

- `assay_run_id`: Optional assay run ID to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Usage Examples:**

```bash
# List all assay runs
python mtcli.py assays runs list

# Get specific assay run
python mtcli.py assays runs list 456
```

#### Create Assay Runs from CSV

```bash
python mtcli.py assays runs create <csv_file> [OPTIONS]
```

Creates assay runs from a CSV file using the `/v1/assay_runs/` endpoint.

**Arguments:**

- `csv_file`: Path to the CSV file containing assay run data

**Options:**

- `--mapping, -m TEXT`: Path to the JSON mapping file (optional)
- `--rows, -r INTEGER`: Number of data rows to process (excludes header row)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--error-handling, -e [reject_all|reject_row]`: Error handling strategy (default: reject_all)
- `--output-format, -o [json|csv]`: Output format (default: json)
- `--dry-run`: Validate data without sending to server

### Assay Results Commands

#### List Assay Results

```bash
python mtcli.py assays results list [assay_result_id] [OPTIONS]
```

Lists assay results using the v1 endpoint.

**Arguments:**

- `assay_result_id`: Optional assay result ID to retrieve

**Options:**

- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])

**Usage Examples:**

```bash
# List all assay results
python mtcli.py assays results list

# Get specific assay result
python mtcli.py assays results list 789
```

#### Create Assay Results from CSV

```bash
python mtcli.py assays results create <csv_file> [OPTIONS]
```

Creates assay results from a CSV file using the `/v1/assay_results/` endpoint.

**Arguments:**

- `csv_file`: Path to the CSV file containing assay result data

**Options:**

- `--mapping, -m TEXT`: Path to the JSON mapping file (optional)
- `--rows, -r INTEGER`: Number of data rows to process (excludes header row)
- `--url TEXT`: Server URL (default: [http://127.0.0.1:8000])
- `--error-handling, -e [reject_all|reject_row]`: Error handling strategy (default: reject_all)
- `--output-format, -o [json|csv]`: Output format (default: json)
- `--dry-run`: Validate data without sending to server

**Usage Examples:**

```bash
# List all assays
python mtcli.py assays list

# Get specific assay
python mtcli.py assays list 123

# Create assays from JSON
python mtcli.py assays create assays.json

# List all assay runs
python mtcli.py assays runs list

# Get specific assay run
python mtcli.py assays runs list 456

# Create assay runs from CSV
python mtcli.py assays runs create runs.csv --mapping runs_mapping.json

# List all assay results
python mtcli.py assays results list

# Get specific assay result
python mtcli.py assays results list 789

# Create assay results from CSV
python mtcli.py assays results create results.csv --mapping results_mapping.json

# Dry run to validate data
python mtcli.py assays runs create runs.csv --dry-run
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
python mtcli.py schema create schema.json
```

### 2. Compound Registration

```bash
# Register compounds from CSV
python mtcli.py compound create compounds.csv --mapping compound_mapping.json
```

### 3. Batch Registration

```bash
# Register batches from CSV
python mtcli.py batch create batches.csv --mapping batch_mapping.json
```

### 4. Additions Registration

```bash
# Register additions from CSV
python mtcli.py additions create additions.csv --mapping additions_mapping.json
```

### 5. Assays Registration

```bash
# Create assays from JSON
python mtcli.py assays create assays.json

# Create assay runs from CSV
python mtcli.py assays runs create runs.csv --mapping runs_mapping.json

# Create assay results from CSV
python mtcli.py assays results create results.csv --mapping results_mapping.json
```

### 6. Data Validation

```bash
# Validate data without sending to server
python mtcli.py compound create compounds.csv --dry-run
```

### 7. Limited Data Processing

```bash
# Process only first 10 rows for testing
python mtcli.py batch create batches.csv --rows 10
```

### 8. Batch Information Retrieval

```bash
# Get batch details
python mtcli.py batch list 123

# Get batch properties
python mtcli.py batch list properties 123

# Get batch synonyms
python mtcli.py batch list synonyms 123

# Get batch additions
python mtcli.py batch list additions 123
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
python mtcli.py compound create data.csv --dry-run
```

### Getting Help

Use the `--help` option to get detailed information about any command:

```bash
python mtcli.py --help
python mtcli.py compound --help
python mtcli.py batch --help
```

## Environment Variables

The client uses the following default configuration:

- **Server URL**: [http://127.0.0.1:8000]
- **Error Handling**: reject_all
- **Output Format**: json

These can be overridden using command-line options.

## Search Output Parameter Formats

The search commands (`search compounds`, `search batches`, `search assay-results`) now support multiple formats for the `--output` parameter.

### Supported Formats

#### 1. Comma-separated string (original format)
```bash
mtcli.py search compounds --output "id,canonical_smiles,common_name"
```

#### 2. JSON file with object containing "output" key
```json
{
  "output": [
    "id",
    "canonical_smiles", 
    "common_name",
    "created_at"
  ]
}
```
```bash
mtcli.py search compounds --output output.json
```

#### 3. JSON file with simple list
```json
[
  "id",
  "batch_regno",
  "compound_id", 
  "notes"
]
```
```bash
mtcli.py search batches --output output_list.json
```

### Examples

#### Compounds Search
```bash
# String format
mtcli.py search compounds --output "id,canonical_smiles" --filter '{"field": "compounds.details.common_name", "operator": "=", "value": "Aspirin"}'

# JSON file format
mtcli.py search compounds --output example_output.json --filter filter.json --output-format table
```

#### Batches Search
```bash
# String format
mtcli.py search batches --output "id,batch_regno,notes" --filter '{"field": "batches.compound_id", "operator": "=", "value": 1}'

# JSON file format
mtcli.py search batches --output example_output_list.json --filter batch_filter.json
```

#### Assay Results Search
```bash
# String format
mtcli.py search assay-results --output "id,value_num,assay_id" --filter '{"field": "assay_results.value_num", "operator": ">", "value": 50}'

# JSON file format
mtcli.py search assay-results --output output.json --filter assay_filter.json --output-format table
```

### File Detection

The system automatically detects whether the `--output` parameter is a file path or a string by checking if the path exists as a file. If the file exists, it's treated as a JSON file; otherwise, it's treated as a comma-separated string.

### Error Handling

- If a JSON file is specified but doesn't exist, an error is shown
- If a JSON file exists but has invalid format, an error is shown
- The JSON file must contain either:
  - A list of column names
  - An object with an "output" key containing a list
  - An object with a "columns" key containing a list
