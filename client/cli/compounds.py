import json
import typer
import requests
from client.utils.api_helpers import get_request_and_porocess
from client.utils.data_ingest import report_csv_information, send_csv_upload_request
from client.utils.display import display_compounds_table, display_properties_table
from client.config.settings import settings
from client.utils.file_utils import load_and_validate_mapping, validate_and_load_csv_data

compound_app = typer.Typer()


# Compound Commands
@compound_app.command("list")
def list_compounds_group(
    skip: int = 0,
    limit: int = 10,
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    List compounds using the v1 endpoint.
    """
    url = f"{url}/v1/compounds/?skip={skip}&limit={limit}"
    get_request_and_porocess(url, None, display_compounds_table, output_format)
    # response = requests.get(f"{url}/v1/compounds/?skip={skip}&limit={limit}")

    # if response.status_code == 200:
    #     compounds_data = response.json()

    #     if output_format == "json":
    #         typer.echo(json.dumps(compounds_data, indent=2))
    #     else:
    #         # Default to table format
    #         display_compounds_table(compounds_data)
    # else:
    #     typer.secho(f"Error: {response.status_code}", fg=typer.colors.RED, err=True)
    #     try:
    #         error_detail = response.json()
    #         typer.secho(f"Details: {json.dumps(error_detail, indent=2)}", fg=typer.colors.RED, err=True)
    #     except (json.JSONDecodeError, ValueError):
    #         typer.secho(f"Response: {response.text}", fg=typer.colors.RED, err=True)
    #     raise typer.Exit(code=1)


@compound_app.command("get")
def get_compound(
    corporate_compound_id: str = typer.Argument(..., help="Corporate Compound ID (friendly name) to retrieve"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get a specific compound by corporate_compound_id (friendly name).
    """
    params = {"property_value": corporate_compound_id, "property_name": "corporate_compound_id"}
    response = requests.get(f"{url}/v1/compounds", params=params)

    if response.status_code == 200:
        compound_data = response.json()

        if output_format == "json":
            print(json.dumps(compound_data, indent=2))
        else:
            # Display single compound in table format
            display_compounds_table([compound_data])
            display_properties_table(compound_data["properties"], display_value=True)
    else:
        typer.secho(f"Error: {response.status_code}", fg=typer.colors.RED, err=True)
        try:
            error_detail = response.json()
            typer.secho(f"Details: {json.dumps(error_detail, indent=2)}", fg=typer.colors.RED, err=True)
        except (json.JSONDecodeError, ValueError):
            typer.secho(f"Response: {response.text}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


@compound_app.command("properties")
def get_compound_properties(
    corporate_compound_id: str = typer.Argument(
        ..., help="Corporate Compound ID (friendly name) to get properties for"
    ),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get all properties for a specific compound by corporate_compound_id (friendly name).
    """
    params = {"property_value": corporate_compound_id, "property_name": "corporate_compound_id"}
    response = requests.get(f"{url}/v1/compounds", params=params)

    if response.status_code == 200:
        properties_data = response.json()

        if output_format == "json":
            print(json.dumps(properties_data, indent=2))
        else:
            # Display properties in table format
            display_properties_table(properties_data["properties"], display_value=True)
    else:
        typer.secho(f"Error: {response.status_code}", fg=typer.colors.RED, err=True)
        try:
            error_detail = response.json()
            typer.secho(f"Details: {json.dumps(error_detail, indent=2)}", fg=typer.colors.RED, err=True)
        except (json.JSONDecodeError, ValueError):
            typer.secho(f"Response: {response.text}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


@compound_app.command("synonyms")
def get_compound_synonyms(
    corporate_compound_id: str = typer.Argument(..., help="Corporate Compound ID (friendly name) to get synonyms for"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get all synonyms for a specific compound by corporate_compound_id (friendly name).
    """
    params = {"property_value": corporate_compound_id, "property_name": "corporate_compound_id"}
    response = requests.get(f"{url}/v1/compounds/synonyms", params=params)

    if response.status_code == 200:
        synonyms_data = response.json()

        if output_format == "json":
            print(json.dumps(synonyms_data, indent=2))
        else:
            # Display synonyms in table format
            display_properties_table(synonyms_data, display_value=True)
    else:
        typer.secho(f"Error: {response.status_code}", fg=typer.colors.RED, err=True)
        try:
            error_detail = response.json()
            typer.secho(f"Details: {json.dumps(error_detail, indent=2)}", fg=typer.colors.RED, err=True)
        except (json.JSONDecodeError, ValueError):
            typer.secho(f"Response: {response.text}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


@compound_app.command("delete")
def delete_compound(
    corporate_compound_id: str = typer.Argument(..., help="Corporate Compound ID (friendly name) to delete"),
    url: str = settings.API_BASE_URL,
):
    """
    Delete a specific compound by corporate_compound_id (friendly name).
    """
    response = requests.delete(f"{url}/v1/compounds/{corporate_compound_id}")
    response_dict = response.json()
    if response.status_code == 200:
        typer.echo(f"✅ {response_dict['message']}")
    else:
        typer.secho(
            f"❌ Error updating the rule. Error message:\n {json.dumps(response_dict, indent=2)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)


@compound_app.command("load")
def add_compounds_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing compound data"),
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
        save_errors=save_errors,
    )
