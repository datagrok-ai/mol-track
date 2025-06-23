# client.py
import sys
import csv
import json
from pathlib import Path
import typer
import requests

# Add parent directory to Python path for imports
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

try:
    from models import SchemaPayload
except ImportError:
    # Fallback for when models is not in the path
    SchemaPayload = None


app = typer.Typer()
schema_app = typer.Typer()
compound_app = typer.Typer()
batch_app = typer.Typer()
batch_list_app = typer.Typer()
properties_app = typer.Typer()
additions_app = typer.Typer()
additions_list_app = typer.Typer()
assays_app = typer.Typer()
assays_runs_app = typer.Typer()
assays_results_app = typer.Typer()
app.add_typer(schema_app, name="schema", help="Schema management commands")
app.add_typer(compound_app, name="compound", help="Compound management commands")
app.add_typer(batch_app, name="batch", help="Batch management commands")
app.add_typer(properties_app, name="properties", help="Properties management commands")
app.add_typer(additions_app, name="additions", help="Additions management commands")
app.add_typer(assays_app, name="assays", help="Assays management commands")
batch_app.add_typer(batch_list_app, name="list", help="List batch information")
additions_app.add_typer(additions_list_app, name="list", help="List additions information")
assays_app.add_typer(assays_runs_app, name="runs", help="Assay runs management")
assays_app.add_typer(assays_results_app, name="results", help="Assay results management")

DEFAULT_SERVER_URL = "http://127.0.0.1:8000"


def validate_file_exists(file_path: str) -> Path:
    """Validate that a file exists and is a file."""
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        typer.echo(f"Error: File '{file_path}' does not exist.", err=True)
        raise typer.Exit(1)

    if not file_path_obj.is_file():
        typer.echo(f"Error: '{file_path}' is not a file.", err=True)
        raise typer.Exit(1)

    return file_path_obj


def load_and_validate_json(file_path: str, model_class=None) -> dict:
    """Load and validate JSON from a file, optionally using a Pydantic model."""
    file_path_obj = validate_file_exists(file_path)

    # Read and parse JSON file
    try:
        with open(file_path_obj, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        typer.echo(f"Error: Invalid JSON format in file '{file_path}': {e}", err=True)
        raise typer.Exit(1)
    except Exception as e:
        typer.echo(f"Error reading file '{file_path}': {e}", err=True)
        raise typer.Exit(1)

    # Validate with Pydantic model if provided
    if model_class:
        if model_class is None:
            typer.echo("⚠️  Warning: Model class is not available for validation. Skipping validation.", err=True)
            return data

        try:
            validated_data = model_class(**data)
            model_name = model_class.__name__
            typer.echo(f"✅ JSON validation passed using {model_name} model!")
            return validated_data.model_dump()
        except Exception as e:
            model_name = model_class.__name__
            typer.echo(f"❌ JSON validation failed using {model_name} model: {e}", err=True)
            raise typer.Exit(1)
    else:
        typer.echo("✅ JSON loaded successfully (no validation model specified)")

    return data


def load_csv_data(file_path: str, max_rows: int | None = None) -> list[dict[str, str]]:
    """Load CSV data from a file and return as list of dictionaries."""
    file_path_obj = validate_file_exists(file_path)

    try:
        with open(file_path_obj, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            data = list(reader)

            if max_rows is not None:
                data = data[:max_rows]
                typer.echo(f"📊 Limited to {len(data)} data rows (requested: {max_rows})")

            return data
    except Exception as e:
        typer.echo(f"Error reading CSV file '{file_path}': {e}", err=True)
        raise typer.Exit(1)


def map_csv_row_to_compound(row: dict[str, str], mapping: dict[str, str]) -> dict:
    """Map a CSV row to compound data using the provided mapping."""
    compound_data = {}

    for csv_field, model_field in mapping.items():
        if csv_field in row:
            value = row[csv_field].strip() if row[csv_field] else None
            if value is not None:
                compound_data[model_field] = value

    return compound_data


def print_response(response):
    if response.status_code == 200:
        print(response.json())
    else:
        print(f"Error: {response.status_code}")


# New utility functions for CSV processing
def validate_and_load_csv_data(
    csv_file: str, entity_type: str, max_rows: int | None = None
) -> tuple[Path, list[dict[str, str]]]:
    """
    Validate CSV file and load data for processing.

    Args:
        csv_file: Path to the CSV file
        entity_type: Type of entity (e.g., "compounds", "batches") for user messages
        max_rows: Maximum number of data rows to load (optional)

    Returns:
        Tuple of (csv_path, csv_data)
    """
    csv_path = validate_file_exists(csv_file)

    typer.echo(f"Loading CSV data from '{csv_file}'...")
    csv_data = load_csv_data(csv_file, max_rows)

    if not csv_data:
        typer.echo("Error: CSV file is empty or has no data rows.", err=True)
        raise typer.Exit(1)

    return csv_path, csv_data


def load_and_validate_mapping(mapping_file: str | None) -> dict[str, any] | None:
    """
    Load and validate mapping file if provided.

    Args:
        mapping_file: Path to the mapping JSON file (optional)

    Returns:
        Mapping data dictionary or None
    """
    mapping_data = None
    if mapping_file:
        typer.echo(f"Loading mapping from '{mapping_file}'...")
        mapping_data = load_and_validate_json(mapping_file)

        if not isinstance(mapping_data, dict):
            typer.echo("Error: Mapping must be a JSON object.", err=True)
            raise typer.Exit(1)

    return mapping_data


def report_csv_information(csv_data: list[dict[str, str]], mapping_data: dict[str, any] | None, entity_type: str):
    """
    Report information about CSV data and mapping.

    Args:
        csv_data: Loaded CSV data
        mapping_data: Optional mapping data
        entity_type: Type of entity for user messages
    """
    csv_headers = set(csv_data[0].keys())
    typer.echo(f"✅ Found {len(csv_data)} {entity_type} in CSV file")
    typer.echo(f"📋 CSV headers: {csv_headers}")

    if mapping_data:
        mapped_headers = set(mapping_data.keys())
        unmapped_headers = csv_headers - mapped_headers
        typer.echo(f"📋 Mapped columns: {len(mapped_headers)}")
        if unmapped_headers:
            typer.echo(f"⚠️  Unmapped columns (will be ignored): {unmapped_headers}")


def send_csv_upload_request(
    csv_path: Path,
    mapping_data: dict[str, any] | None,
    url: str,
    endpoint: str,
    error_handling: str,
    output_format: str,
    entity_type: str,
    csv_data: list[dict[str, str]] | None = None,
) -> None:
    """
    Send CSV upload request to the specified endpoint.

    Args:
        csv_path: Path to the CSV file
        mapping_data: Optional mapping data
        url: Server URL
        endpoint: API endpoint to send request to
        error_handling: Error handling strategy
        output_format: Output format
        entity_type: Type of entity for user messages
        csv_data: Optional CSV data to use instead of reading from file
    """
    try:
        if csv_data is not None:
            # Create temporary CSV file with limited data
            import tempfile

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".csv", delete=False, newline="", encoding="utf-8"
            ) as temp_file:
                if csv_data:
                    writer = csv.DictWriter(temp_file, fieldnames=csv_data[0].keys())
                    writer.writeheader()
                    writer.writerows(csv_data)
                    temp_file_path = temp_file.name
                else:
                    temp_file_path = None
        else:
            temp_file_path = None

        try:
            # Use temporary file if created, otherwise use original
            file_to_send = temp_file_path if temp_file_path else csv_path

            with open(file_to_send, "rb") as f:
                files = {"csv_file": (csv_path.name, f, "text/csv")}

                data = {"error_handling": error_handling, "output_format": output_format}

                if mapping_data:
                    data["mapping"] = json.dumps(mapping_data)

                typer.echo(f"🚀 Sending {csv_path.name} to {url}{endpoint}...")

                response = requests.post(f"{url}{endpoint}", files=files, data=data)

                if response.status_code == 200:
                    result = response.json()
                    typer.echo(f"✅ {entity_type.capitalize()} registered successfully!")

                    # Parse the result based on output format
                    if output_format == "json":
                        if "data" in result:
                            data_list = result["data"]
                            success_count = sum(1 for item in data_list if item.get("registration_status") == "success")
                            error_count = len(data_list) - success_count
                            typer.echo(f"📊 Results: {success_count} successful, {error_count} errors")

                            # Show any errors
                            errors = [item for item in data_list if item.get("registration_status") != "success"]
                            if errors:
                                typer.echo("❌ Errors found:")
                                for error in errors[:5]:  # Show first 5 errors
                                    typer.echo(
                                        f"  - Row {error.get('row', 'unknown')}: {error.get('registration_error_message', 'Unknown error')}"
                                    )
                                if len(errors) > 5:
                                    typer.echo(f"  ... and {len(errors) - 5} more errors")
                        else:
                            typer.echo("✅ Registration completed")
                    else:
                        typer.echo("✅ CSV output received")

                else:
                    typer.echo(f"❌ Error: {response.status_code}")
                    try:
                        error_detail = response.json()
                        typer.echo(f"Details: {json.dumps(error_detail, indent=2)}")
                    except (json.JSONDecodeError, ValueError):
                        typer.echo(f"Response: {response.text}")
        finally:
            # Clean up temporary file
            if temp_file_path and temp_file_path != str(csv_path):
                try:
                    import os

                    os.unlink(temp_file_path)
                except Exception:
                    pass  # Ignore cleanup errors

    except requests.exceptions.ConnectionError:
        typer.echo(f"❌ Error: Could not connect to server at {url}", err=True)
        raise typer.Exit(1)
    except requests.exceptions.RequestException as e:
        typer.echo(f"❌ Error making request: {e}", err=True)
        raise typer.Exit(1)
    except Exception as e:
        typer.echo(f"❌ Error: {e}", err=True)
        raise typer.Exit(1)


# Schema Commands
@schema_app.command("list")
def list_schema_group(url: str = DEFAULT_SERVER_URL):
    """
    List the schema.
    """
    print(f"Listing schema from {url}")
    response = requests.get(f"{url}/v1/schema/compounds")
    typer.echo(f"Compounds schema: {response.json()}")
    response = requests.get(f"{url}/v1/schema/batches")
    typer.echo(f"Batches schema: {response.json()}")


@schema_app.command("create")
def add_schema_from_file(
    file_path: str = typer.Argument(..., help="Path to the JSON file containing schema data"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Add schema from a JSON file.

    The file should contain a JSON object with 'properties' and optionally 'synonym_types' arrays.
    Example format:
    {
        "properties": [
            {
                "name": "Molecular Weight",
                "value_type": "number",
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
    """
    # Load and validate schema using utility function
    schema_data = load_and_validate_json(file_path, SchemaPayload)

    # Send request to server
    try:
        typer.echo(f"Adding schema from file '{file_path}' to {url}...")
        response = requests.post(f"{url}/v1/schema", json=schema_data)

        if response.status_code == 200:
            result = response.json()
            typer.echo("✅ Schema added successfully!")

            if "created" in result:
                created = result["created"]
                if "properties" in created and created["properties"]:
                    typer.echo(f"📋 Created {len(created['properties'])} properties")
                if "synonym_types" in created and created["synonym_types"]:
                    typer.echo(f"🏷️  Created {len(created['synonym_types'])} synonym types")
        else:
            typer.echo(f"❌ Error: {response.status_code}")
            try:
                error_detail = response.json()
                typer.echo(f"Details: {json.dumps(error_detail, indent=2)}")
            except (json.JSONDecodeError, ValueError):
                typer.echo(f"Response: {response.text}")

    except requests.exceptions.ConnectionError:
        typer.echo(f"❌ Error: Could not connect to server at {url}", err=True)
        raise typer.Exit(1)
    except requests.exceptions.RequestException as e:
        typer.echo(f"❌ Error making request: {e}", err=True)
        raise typer.Exit(1)


# Compound Commands
@compound_app.command("list")
def list_compounds_group(skip: int = 0, limit: int = 10, url: str = DEFAULT_SERVER_URL):
    """
    List compounds using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/compounds/?skip={skip}&limit={limit}")
    print_response(response)


@compound_app.command("create_old")
def create_compound_old(item_json: str, url: str = DEFAULT_SERVER_URL):
    """
    Create a compound using the legacy /compounds endpoint.
    """
    try:
        item = json.loads(item_json)
    except json.JSONDecodeError:
        print("Invalid JSON format for item.")
        return

    response = requests.post(f"{url}/compounds/", json=item)
    print_response(response)


@compound_app.command("create")
def add_compounds_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing compound data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json or csv"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
):
    """
    Add compounds from a CSV file using the /v1/compounds/ endpoint.

    The CSV file should contain compound data with headers. If a mapping file is provided,
    it should be a JSON object mapping CSV column names to field names.

    Example mapping:
    {
        "structure": "smiles",
        "common_name": "compounds_details.common_name",
        "cas": "compounds_details.cas"
    }

    If no mapping is provided, the system will attempt to infer the mapping from column names.

    Use --rows to limit the number of data rows processed (header row is excluded).
    """

    # Use utility functions for common steps
    csv_path, csv_data = validate_and_load_csv_data(csv_file, "compounds", rows)
    mapping_data = load_and_validate_mapping(mapping_file)
    report_csv_information(csv_data, mapping_data, "compounds")

    if dry_run:
        typer.echo("✅ Dry run completed successfully!")
        return

    # Send the request using utility function
    send_csv_upload_request(
        csv_path=csv_path,
        mapping_data=mapping_data,
        url=url,
        endpoint="/v1/compounds/",
        error_handling=error_handling,
        output_format=output_format,
        entity_type="compounds",
        csv_data=csv_data if rows is not None else None,
    )


# List Commands
@properties_app.command("list")
def list_properties(url: str = DEFAULT_SERVER_URL):
    """
    List the properties.
    """
    response = requests.get(f"{url}/properties")
    print_response(response)


# Batch Commands
@batch_app.command("list")
def list_batches_group(
    batch_id: int | None = typer.Argument(None, help="Batch ID to retrieve (optional)"),
    skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
    limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    List batches using the v1 endpoint.

    If no batch_id is provided, lists all batches with pagination.
    If batch_id is provided, gets the specific batch.
    """
    if batch_id is not None:
        # Get specific batch
        response = requests.get(f"{url}/v1/batches/{batch_id}")
    else:
        # List all batches
        response = requests.get(f"{url}/v1/batches/?skip={skip}&limit={limit}")
    print_response(response)


@batch_app.command("create")
def create_batches_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing batch data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json or csv"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
):
    """
    Add batches from a CSV file using the /v1/batches/ endpoint.

    The CSV file should contain batch data with headers. If a mapping file is provided,
    it should be a JSON object mapping CSV column names to field names.

    Example mapping:
    {
        "structure": "smiles",
        "batch_corporate_id": "batches_details.batch_corporate_id",
        "common_name": "compounds_details.common_name",
        "cas": "compounds_details.cas"
    }

    If no mapping is provided, the system will attempt to infer the mapping from column names.

    Use --rows to limit the number of data rows processed (header row is excluded).
    """

    # Use utility functions for common steps
    csv_path, csv_data = validate_and_load_csv_data(csv_file, "batches", rows)
    mapping_data = load_and_validate_mapping(mapping_file)
    report_csv_information(csv_data, mapping_data, "batches")

    if dry_run:
        typer.echo("✅ Dry run completed successfully!")
        return

    # Send the request using utility function
    send_csv_upload_request(
        csv_path=csv_path,
        mapping_data=mapping_data,
        url=url,
        endpoint="/v1/batches/",
        error_handling=error_handling,
        output_format=output_format,
        entity_type="batches",
        csv_data=csv_data if rows is not None else None,
    )


# Batch List Commands
@batch_list_app.command("properties")
def get_batch_properties(
    batch_id: int = typer.Argument(..., help="Batch ID to get properties for"), url: str = DEFAULT_SERVER_URL
):
    """
    Get all properties for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/properties")
    print_response(response)


@batch_list_app.command("synonyms")
def get_batch_synonyms(
    batch_id: int = typer.Argument(..., help="Batch ID to get synonyms for"), url: str = DEFAULT_SERVER_URL
):
    """
    Get all synonyms for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/synonyms")
    print_response(response)


@batch_list_app.command("additions")
def get_batch_additions(
    batch_id: int = typer.Argument(..., help="Batch ID to get additions for"), url: str = DEFAULT_SERVER_URL
):
    """
    Get all additions for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/additions")
    print_response(response)


# Additions Commands
@additions_app.command("list")
def list_additions(url: str = DEFAULT_SERVER_URL):
    """
    List the additions.
    """
    response = requests.get(f"{url}/v1/additions")
    print_response(response)


@additions_app.command("update")
def update_addition(
    addition_id: int = typer.Argument(..., help="Addition ID to update"),
    file_path: str = typer.Argument(..., help="Path to the JSON file containing addition update data"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Update information for the specified addition.
    """
    # Load and validate the update data
    update_data = load_and_validate_json(file_path)

    response = requests.put(f"{url}/v1/additions/{addition_id}", json=update_data)
    print_response(response)


@additions_app.command("delete")
def delete_addition(
    addition_id: int = typer.Argument(..., help="Addition ID to delete"), url: str = DEFAULT_SERVER_URL
):
    """
    Soft delete the specified addition (only if no dependent batches exist).
    """
    response = requests.delete(f"{url}/v1/additions/{addition_id}")
    print_response(response)


@additions_app.command("create")
def add_additions_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing addition data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json or csv"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
):
    """
    Add additions from a CSV file using the /v1/additions/ endpoint.

    The CSV file should contain addition data with headers. If a mapping file is provided,
    it should be a JSON object mapping CSV column names to field names.

    Example mapping:
    {
        "structure": "smiles",
        "common_name": "additions_details.common_name",
        "cas": "additions_details.cas"
    }

    If no mapping is provided, the system will attempt to infer the mapping from column names.

    Use --rows to limit the number of data rows processed (header row is excluded).
    """

    # Use utility functions for common steps
    csv_path, csv_data = validate_and_load_csv_data(csv_file, "additions", rows)
    mapping_data = load_and_validate_mapping(mapping_file)
    report_csv_information(csv_data, mapping_data, "additions")

    if dry_run:
        typer.echo("✅ Dry run completed successfully!")
        return

    # Send the request using utility function
    send_csv_upload_request(
        csv_path=csv_path,
        mapping_data=mapping_data,
        url=url,
        endpoint="/v1/additions/",
        error_handling=error_handling,
        output_format=output_format,
        entity_type="additions",
        csv_data=csv_data if rows is not None else None,
    )


# Additions List Commands
@additions_list_app.command("salts")
def list_additions_salts(url: str = DEFAULT_SERVER_URL):
    """
    List all additions with role of salts.
    """
    response = requests.get(f"{url}/v1/additions/salts")
    print_response(response)


@additions_list_app.command("solvates")
def list_additions_solvates(url: str = DEFAULT_SERVER_URL):
    """
    List all additions with role of solvates.
    """
    response = requests.get(f"{url}/v1/additions/solvates")
    print_response(response)


@additions_list_app.command()
def get_addition(addition_id: int = typer.Argument(..., help="Addition ID to retrieve"), url: str = DEFAULT_SERVER_URL):
    """
    Get all information for a specific addition.
    """
    response = requests.get(f"{url}/v1/additions/{addition_id}")
    print_response(response)


# Assays Commands
@assays_app.command("list")
def list_assays(
    assay_id: int | None = typer.Argument(None, help="Assay ID to retrieve (optional)"), url: str = DEFAULT_SERVER_URL
):
    """
    List assays using the v1 endpoint.

    If no assay_id is provided, lists all assays.
    If assay_id is provided, gets the specific assay.
    """
    if assay_id is not None:
        # Get specific assay
        response = requests.get(f"{url}/v1/assays/{assay_id}")
    else:
        # List all assays
        response = requests.get(f"{url}/v1/assays")
    print_response(response)


@assays_app.command("create")
def create_assays(
    file_path: str = typer.Argument(..., help="Path to the JSON file containing assay data"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Create assays from a JSON file using the /v1/assays endpoint.
    """
    # Load and validate the assay data
    assay_data = load_and_validate_json(file_path)

    response = requests.post(f"{url}/v1/assays", json=assay_data)
    print_response(response)


# Assay Runs Commands
@assays_runs_app.command("list")
def list_assay_runs(
    assay_run_id: int | None = typer.Argument(None, help="Assay run ID to retrieve (optional)"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    List assay runs using the v1 endpoint.

    If no assay_run_id is provided, lists all assay runs.
    If assay_run_id is provided, gets the specific assay run.
    """
    if assay_run_id is not None:
        # Get specific assay run
        response = requests.get(f"{url}/v1/assay_runs/{assay_run_id}")
    else:
        # List all assay runs
        response = requests.get(f"{url}/v1/assay_runs")
    print_response(response)


@assays_runs_app.command("create")
def create_assay_runs_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing assay run data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json or csv"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
):
    """
    Create assay runs from a CSV file using the /v1/assay_runs/ endpoint.
    """
    # Use utility functions for common steps
    csv_path, csv_data = validate_and_load_csv_data(csv_file, "assay runs", rows)
    mapping_data = load_and_validate_mapping(mapping_file)
    report_csv_information(csv_data, mapping_data, "assay runs")

    if dry_run:
        typer.echo("✅ Dry run completed successfully!")
        return

    # Send the request using utility function
    send_csv_upload_request(
        csv_path=csv_path,
        mapping_data=mapping_data,
        url=url,
        endpoint="/v1/assay_runs/",
        error_handling=error_handling,
        output_format=output_format,
        entity_type="assay runs",
        csv_data=csv_data if rows is not None else None,
    )


# Assay Results Commands
@assays_results_app.command("list")
def list_assay_results(
    assay_result_id: int | None = typer.Argument(None, help="Assay result ID to retrieve (optional)"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    List assay results using the v1 endpoint.

    If no assay_result_id is provided, lists all assay results.
    If assay_result_id is provided, gets the specific assay result.
    """
    if assay_result_id is not None:
        # Get specific assay result
        response = requests.get(f"{url}/v1/assay_results/{assay_result_id}")
    else:
        # List all assay results
        response = requests.get(f"{url}/v1/assay_results")
    print_response(response)


@assays_results_app.command("create")
def create_assay_results_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing assay result data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json or csv"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
):
    """
    Create assay results from a CSV file using the /v1/assay_results/ endpoint.
    """
    # Use utility functions for common steps
    csv_path, csv_data = validate_and_load_csv_data(csv_file, "assay results", rows)
    mapping_data = load_and_validate_mapping(mapping_file)
    report_csv_information(csv_data, mapping_data, "assay results")

    if dry_run:
        typer.echo("✅ Dry run completed successfully!")
        return

    # Send the request using utility function
    send_csv_upload_request(
        csv_path=csv_path,
        mapping_data=mapping_data,
        url=url,
        endpoint="/v1/assay_results/",
        error_handling=error_handling,
        output_format=output_format,
        entity_type="assay results",
        csv_data=csv_data if rows is not None else None,
    )


if __name__ == "__main__":
    app()
