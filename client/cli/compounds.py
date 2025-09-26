import typer
from client.config.settings import settings
from client.utils.display import display_compounds_table
from app.utils import enums
from client.cli.shared import EntityCLI

compound_app = typer.Typer()


class CompoundsCLI(EntityCLI):
    entity_type = "compounds"
    display_fn = staticmethod(display_compounds_table)

    def get_endpoint(self) -> str:
        return "v1/compounds"


compound_cli = CompoundsCLI()


@compound_app.command("list")
def list_compounds_group(
    skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
    limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    compound_cli.list(skip=skip, limit=limit, url=url, output_format=output_format, output_file=output_file)


@compound_app.command("get")
def get_compound(
    corporate_compound_id: str = typer.Argument(..., help="Corporate Compound ID (friendly name) to retrieve"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    compound_cli.get(corporate_compound_id, url=url, output_format=output_format, output_file=output_file)


@compound_app.command("delete")
def delete_compound(
    corporate_compound_id: str = typer.Argument(..., help="Corporate Compound ID (friendly name) to delete"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
):
    compound_cli.delete(corporate_compound_id, url=url)


@compound_app.command("load")
def add_compounds_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing compound data"),
    mapping_file: str | None = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
    error_handling: enums.ErrorHandlingOptions = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy"
    ),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
):
    compound_cli.load_from_csv(
        csv_file=csv_file,
        mapping_file=mapping_file,
        rows=rows,
        url=url,
        error_handling=error_handling,
        dry_run=dry_run,
        save_errors=save_errors,
    )
