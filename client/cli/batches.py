import json
import typer
import requests
from client.config.settings import settings
from client.utils.api_helpers import print_response
from client.utils.data_ingest import report_csv_information, send_csv_upload_request
from client.utils.display import display_batches_table
from client.utils.file_utils import load_and_validate_mapping, validate_and_load_csv_data

batch_app = typer.Typer()
batch_list_app = typer.Typer()
batch_app.add_typer(batch_list_app, name="list", help="List batch information")


# Batch Commands
@batch_app.command("list")
def list_batches_group(
    batch_id: int | None = typer.Argument(None, help="Batch ID to retrieve (optional)"),
    skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
    limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
    url: str = settings.API_BASE_URL,
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


@batch_app.command("load")
def create_batches_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing batch data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = settings.API_BASE_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json or csv"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
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
        save_errors=save_errors,
    )


@batch_app.command("get")
def get_batch(
    batch_id: int = typer.Argument(..., help="Batch ID to retrieve"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get a specific batch by ID.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}")

    if response.status_code == 200:
        batch_data = response.json()

        if output_format == "json":
            print(json.dumps(batch_data, indent=2))
        else:
            # Display single batch in table format
            display_batches_table([batch_data])
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


@batch_app.command("delete")
def delete_batch(
    batch_id: int = typer.Argument(..., help="Batch ID to delete"),
    url: str = settings.API_BASE_URL,
):
    """
    Delete a specific batch by ID.
    """
    response = requests.delete(f"{url}/v1/batches/{batch_id}")
    print_response(response)


# Batch List Commands
@batch_list_app.command("properties")
def get_batch_properties(
    batch_id: int = typer.Argument(..., help="Batch ID to get properties for"), url: str = settings.API_BASE_URL
):
    """
    Get all properties for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/properties")
    print_response(response)


@batch_list_app.command("synonyms")
def get_batch_synonyms(
    batch_id: int = typer.Argument(..., help="Batch ID to get synonyms for"), url: str = settings.API_BASE_URL
):
    """
    Get all synonyms for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/synonyms")
    print_response(response)


@batch_list_app.command("additions")
def get_batch_additions(
    batch_id: int = typer.Argument(..., help="Batch ID to get additions for"), url: str = settings.API_BASE_URL
):
    """
    Get all additions for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/additions")
    print_response(response)
