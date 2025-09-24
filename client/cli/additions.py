import json
import typer

from app import models
from client.config import settings
from client.utils.api_helpers import handle_delete_request, handle_get_request, handle_put_request
from client.utils.data_ingest import report_csv_information, send_csv_upload_request
from client.utils.display import display_additions_table
from client.utils.file_utils import (
    load_and_validate_json,
    load_and_validate_mapping,
    validate_and_load_csv_data,
    write_result_to_file,
)


additions_app = typer.Typer()
additions_list_app = typer.Typer()
additions_app.add_typer(additions_list_app, name="list", help="List additions information")


# Additions Commands
@additions_app.command("update")
def update_addition(
    addition_id: int = typer.Argument(..., help="Addition ID to update"),
    file_path: str = typer.Argument(..., help="Path to the JSON file containing addition update data"),
    url: str = settings.API_BASE_URL,
):
    """
    Update information for the specified addition.
    """
    # Load and validate the update data
    update_data = load_and_validate_json(file_path, models.AdditionUpdate)

    endpoint = f"{url}/v1/additions/{addition_id}"
    handle_put_request(endpoint, json_data=update_data)
    typer.echo(f"✅ Addition with id {addition_id} has been successfully updated.")


@additions_app.command("delete")
def delete_addition(
    addition_id: int = typer.Argument(..., help="Addition ID to delete"), url: str = settings.API_BASE_URL
):
    """
    Soft delete the specified addition (only if no dependent batches exist).
    """
    endpoint = f"{url}/v1/additions/{addition_id}"
    handle_delete_request(endpoint)
    typer.echo(f"✅ Addition with id {addition_id} has been successfully deleted.")


@additions_app.command("get")
def get_addition(
    addition_id: int = typer.Argument(..., help="Addition ID to retrieve"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    Get all information for a specific addition.
    """

    endpoint = f"{url}/v1/additions/{addition_id}"
    data = handle_get_request(endpoint)

    if output_format == "json":
        typer.echo(json.dumps(data, indent=2))
    else:
        display_additions_table([data])

    write_result_to_file(data, output_format, output_file)


@additions_app.command("load")
def add_additions_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing addition data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = settings.API_BASE_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
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
        entity_type="additions",
        csv_data=csv_data if rows is not None else None,
        save_errors=save_errors,
    )


# Additions List Commands
@additions_list_app.command("all")
def list_additions(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    List the additions.
    """
    endpoint = f"{url}/v1/additions/"
    data = handle_get_request(endpoint)
    if output_format == "json":
        typer.echo(json.dumps(data, indent=2))
    else:
        display_additions_table(data)
    write_result_to_file(data, output_format, output_file)


@additions_list_app.command("salts")
def list_additions_salts(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    List all additions with role of salts.
    """
    endpoint = f"{url}/v1/additions/salts"
    data = handle_get_request(endpoint)
    if output_format == "json":
        typer.echo(json.dumps(data, indent=2))
    else:
        display_additions_table(data)

    write_result_to_file(data, output_format, output_file)


@additions_list_app.command("solvates")
def list_additions_solvates(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    List all additions with role of solvates.
    """
    endpoint = f"{url}/v1/additions/solvates"
    data = handle_get_request(endpoint)
    if output_format == "json":
        typer.echo(json.dumps(data, indent=2))
    else:
        display_additions_table(data)

    write_result_to_file(data, output_format, output_file)
