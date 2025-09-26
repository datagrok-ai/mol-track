import requests
import typer
from client.cli.shared import EntityCLI
from client.config.settings import settings
from client.utils.api_helpers import print_response
from client.utils.display import display_batches_table

from app.utils import enums

batch_app = typer.Typer()


class BatchCLI(EntityCLI):
    entity_type = "batches"
    display_fn = staticmethod(display_batches_table)

    def get_endpoint(self) -> str:
        return "v1/batches"


batch_cli = BatchCLI()


@batch_app.command("list")
def list_batches_group(
    skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
    limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    batch_cli.list(skip=skip, limit=limit, url=url, output_format=output_format, output_file=output_file)


@batch_app.command("get")
def get_batch(
    corporate_batch_id: str = typer.Argument(..., help="Corporate Batch ID (friendly name) to retrieve"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    batch_cli.get(corporate_batch_id, url=url, output_format=output_format, output_file=output_file)


@batch_app.command("delete")
def delete_batch(
    corporate_batch_id: str = typer.Argument(..., help="Corporate Batch ID (friendly name) to delete"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
):
    batch_cli.delete(corporate_batch_id, url=url)


@batch_app.command("load")
def create_batches_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing batch data"),
    mapping_file: str | None = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    error_handling: enums.ErrorHandlingOptions = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy"
    ),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
):
    batch_cli.load_from_csv(
        csv_file=csv_file,
        mapping_file=mapping_file,
        rows=rows,
        url=url,
        error_handling=error_handling,
        dry_run=dry_run,
        save_errors=save_errors,
    )


@batch_app.command("additions")
def get_batch_additions(
    corporate_batch_id: str = typer.Argument(..., help="Corporate Batch ID (friendly name) to get additions for"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
):
    params = {"property_value": corporate_batch_id, "property_name": "corporate_batch_id"}
    response = requests.get(f"{url}/v1/batches/additions", params=params)
    # TODO: This will probably need to be formatted
    print_response(response)
