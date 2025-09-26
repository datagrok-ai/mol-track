import requests
import typer

from client.config import settings
from client.utils.api_helpers import print_response
from client.utils.display import display_assays_table
from client.utils.file_utils import load_and_validate_json
from client.cli.shared import EntityCLI


assays_app = typer.Typer()
assays_runs_app = typer.Typer()
assays_results_app = typer.Typer()
assays_app.add_typer(assays_runs_app, name="runs", help="Assay runs management")
assays_app.add_typer(assays_results_app, name="results", help="Assay results management")


# Assays Commands
class AssaysCLI(EntityCLI):
    entity_type = "assays"
    display_fn = staticmethod(display_assays_table)

    def get_endpoint(self) -> str:
        return "v1/assays"


assays_cli = AssaysCLI()


@assays_app.command("list")
def list_assays(
    skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
    limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    assays_cli.list(skip=skip, limit=limit, url=url, output_format=output_format, output_file=output_file)


@assays_app.command("get")
def get_assay(
    assay_id: int = typer.Argument(..., help="Assay ID to retrieve"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    # TODO: Update get to support not only the corporate id but also the id
    assays_cli.get(str(assay_id), url=url, output_format=output_format, output_file=output_file)


@assays_app.command("load")
def create_assays(
    file_path: str = typer.Argument(..., help="Path to the JSON file containing assay data"),
    url: str = settings.API_BASE_URL,
):
    """
    Load assays from a JSON file using the /v1/assays endpoint.
    """
    assay_data = load_and_validate_json(file_path)
    response = requests.post(f"{url}/v1/assays", json=assay_data)
    print_response(response)


# Assay Runs Commands
class AssayRunsCLI(EntityCLI):
    entity_type = "assay_runs"

    @staticmethod
    def display_fn(data):
        display_assays_table(data, assay_entity="run")

    def get_endpoint(self) -> str:
        return "v1/assay_runs"


assay_runs_cli = AssayRunsCLI()


@assays_runs_app.command("list")
def list_assay_runs(
    skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
    limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    assay_runs_cli.list(skip=skip, limit=limit, url=url, output_format=output_format, output_file=output_file)


@assays_runs_app.command("load")
def create_assay_runs_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing assay run data"),
    mapping_file: str | None = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
):
    assay_runs_cli.load_from_csv(
        csv_file=csv_file,
        mapping_file=mapping_file,
        rows=rows,
        url=url,
        error_handling=error_handling,
        dry_run=dry_run,
        save_errors=save_errors,
    )


# Assay Results Commands
class AssayResultsCLI(EntityCLI):
    entity_type = "assay results"

    @staticmethod
    def display_fn(data):
        display_assays_table(data, assay_entity="run")

    def get_endpoint(self) -> str:
        return "v1/assay_results"


assay_results_cli = AssayResultsCLI()


@assays_results_app.command("list")
def list_assay_results(
    skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
    limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    assay_results_cli.list(skip=skip, limit=limit, url=url, output_format=output_format, output_file=output_file)


@assays_results_app.command("load")
def create_assay_results_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing assay result data"),
    mapping_file: str | None = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
):
    assay_results_cli.load_from_csv(
        csv_file=csv_file,
        mapping_file=mapping_file,
        rows=rows,
        url=url,
        error_handling=error_handling,
        dry_run=dry_run,
        save_errors=save_errors,
    )
