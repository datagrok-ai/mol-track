import json
from typing import Optional
import typer

from app import models
from client.config import settings
from client.utils.api_helpers import handle_delete_request, handle_get_request, handle_put_request
from client.utils.data_ingest import report_csv_information, send_csv_upload_request
from client.utils.display import display_additions_table
from client.utils.file_utils import (
    load_and_validate_json,
    validate_and_load_csv_data,
    write_result_to_file,
)

ADDITIONS_ENDPOINTS = {
    "all": ("", "List the additions."),
    "salts": ("/salts", "List all additions with role of salts."),
    "solvates": ("/solvates", "List all additions with role of solvates."),
}


def register_additions_commands(app: typer.Typer):
    for cmd_name, (endpoint_suffix, description) in ADDITIONS_ENDPOINTS.items():

        def command_func(
            url: str = settings.API_BASE_URL,
            output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
            output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
            _endpoint_suffix=endpoint_suffix,
        ):
            endpoint = f"{url}/v1/additions{_endpoint_suffix}"
            data = handle_get_request(endpoint)
            output_data(data, output_format, output_file, display_additions_table)

        app.command(cmd_name, help=description)(command_func)


additions_app = typer.Typer()
additions_list_app = typer.Typer()
additions_app.add_typer(additions_list_app, name="list", help="List additions information")
register_additions_commands(additions_list_app)


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
    Hard delete the specified addition (only if no dependent batches exist).
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
    output_data(data, output_format, output_file, display_additions_table)


@additions_app.command("load")
def add_additions_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing addition data"),
    rows: Optional[int] = typer.Option(
        None, "--rows", "-r", help="Number of data rows to process (excludes header row)"
    ),
    url: str = settings.API_BASE_URL,
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
):
    """
    Add additions from a CSV file.

    Args:
        csv_file: Path to the CSV file.
        rows: Optional number of rows to process (excluding header).
        url: API base URL. Defaults to settings.API_BASE_URL.
        dry_run: Validate data without sending it.
        save_errors: Save errors to JSON file.
    """
    # Use utility functions for common steps
    csv_path, csv_data = validate_and_load_csv_data(csv_file, "additions", rows)
    report_csv_information(csv_data, "additions")

    if dry_run:
        typer.echo("✅ Dry run completed successfully!")
        return

    # Send the request using utility function
    send_csv_upload_request(
        csv_path=csv_path,
        url=url,
        endpoint="/v1/additions/",
        entity_type="additions",
        csv_data=csv_data if rows is not None else None,
        save_errors=save_errors,
    )


def output_data(data, output_format="table", output_file=None, table_display_func=None):
    table_data = data if isinstance(data, list) else [data]

    if output_format == "json":
        typer.echo(json.dumps(data, indent=2))
    elif table_display_func:
        table_display_func(table_data)

    write_result_to_file(data, output_format, output_file)
