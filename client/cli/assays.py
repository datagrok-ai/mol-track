import json
import requests
import typer

from client.config import settings
from client.utils.api_helpers import print_response
from client.utils.data_ingest import report_csv_information, send_csv_upload_request
from client.utils.display import display_assays_table, display_properties_table
from client.utils.file_utils import load_and_validate_json, load_and_validate_mapping, validate_and_load_csv_data


assays_app = typer.Typer()
assays_runs_app = typer.Typer()
assays_results_app = typer.Typer()
assays_app.add_typer(assays_runs_app, name="runs", help="Assay runs management")
assays_app.add_typer(assays_results_app, name="results", help="Assay results management")


# Assays Commands
@assays_app.command("list")
def list_assays(
    skip: int = 0,
    limit: int = 10,
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    List assays using the v1 endpoint.

    If no assay_id is provided, lists all assays.
    If assay_id is provided, gets the specific assay.
    """
    response = requests.get(f"{url}/v1/assays/?skip={skip}&limit={limit}")
    if response.status_code == 200:
        assays_data = response.json()

        if output_format == "json":
            typer.echo(json.dumps(assays_data, indent=2))
        else:
            # Default to table format
            display_assays_table(assays_data)
    else:
        typer.secho(f"Error: {response.status_code}", fg=typer.colors.RED, err=True)
        try:
            error_detail = response.json()
            typer.secho(f"Details: {json.dumps(error_detail, indent=2)}", fg=typer.colors.RED, err=True)
        except (json.JSONDecodeError, ValueError):
            typer.secho(f"Response: {response.text}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


@assays_app.command("load")
def create_assays(
    file_path: str = typer.Argument(..., help="Path to the JSON file containing assay data"),
    url: str = settings.API_BASE_URL,
):
    """
    Load assays from a JSON file using the /v1/assays endpoint.
    """
    # Load and validate the assay data
    assay_data = load_and_validate_json(file_path)

    response = requests.post(f"{url}/v1/assays", json=assay_data)
    print_response(response)


@assays_app.command("get")
def get_assay(
    assay_id: int = typer.Argument(..., help="Assay ID to retrieve"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get a specific assay by ID.
    """
    response = requests.get(f"{url}/v1/assays/{assay_id}")

    if response.status_code == 200:
        assay_data = response.json()

        if output_format == "json":
            print(json.dumps(assay_data, indent=2))
        else:
            # Display single assay in table format
            display_assays_table([assay_data])
            display_properties_table(assay_data["properties"], display_value=True)
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


# Assay Runs Commands
@assays_runs_app.command("list")
def list_assay_runs(
    skip: int = 0,
    limit: int = 10,
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    List assay runs using the v1 endpoint.

    If no assay_run_id is provided, lists all assay runs.
    If assay_run_id is provided, gets the specific assay run.
    """
    response = requests.get(f"{url}/v1/assay_runs/?skip={skip}&limit={limit}")
    if response.status_code == 200:
        assay_runs_data = response.json()

        if output_format == "json":
            typer.echo(json.dumps(assay_runs_data, indent=2))
        else:
            # Default to table format
            display_assays_table(assay_runs_data)
    else:
        typer.secho(f"Error: {response.status_code}", fg=typer.colors.RED, err=True)
        try:
            error_detail = response.json()
            typer.secho(f"Details: {json.dumps(error_detail, indent=2)}", fg=typer.colors.RED, err=True)
        except (json.JSONDecodeError, ValueError):
            typer.secho(f"Response: {response.text}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


@assays_runs_app.command("load")
def create_assay_runs_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing assay run data"),
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
    Load assay runs from a CSV file using the /v1/assay_runs/ endpoint.
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
        save_errors=save_errors,
    )


# Assay Results Commands
@assays_results_app.command("list")
def list_assay_results(
    assay_result_id: int | None = typer.Argument(None, help="Assay result ID to retrieve (optional)"),
    url: str = settings.API_BASE_URL,
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


@assays_results_app.command("load")
def create_assay_results_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing assay result data"),
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
    Load assay results from a CSV file using the /v1/assay_results/ endpoint.
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
        save_errors=save_errors,
    )
