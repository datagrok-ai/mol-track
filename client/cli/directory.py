from pathlib import Path
import typer

from client.cli.assays import assays_cli, assay_runs_cli, assay_results_cli
from client.cli.batches import batch_cli
from client.cli.compounds import compound_cli
from client.cli.schema import add_schema_from_file
from client.config import settings


directory_app = typer.Typer()


@directory_app.command("load")
def load_directory(
    directory_path: str = typer.Argument(..., help="Path to the directory containing files to load"),
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    url: str = settings.API_BASE_URL,
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    auto_mapping: bool = typer.Option(False, "--auto-mapping", help="Skip local mapping files and use auto-detection"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to JSON files"),
):
    """
    Load contents from a directory in the correct order.
    """
    directory = Path(directory_path)

    if not directory.exists() or not directory.is_dir():
        typer.echo(f"❌ Error: '{directory_path}' is not a valid directory.", err=True)
        raise typer.Exit(1)

    typer.echo(f"📁 Loading contents from directory: {directory_path}")
    typer.echo("=" * 60)

    # Generic loader function
    def load_file(file_name, loader_fn, mapping_name=None, is_schema=False):
        file_path = directory / file_name
        if not file_path.exists():
            typer.echo(f"⏭️  No {file_name} found - skipping")
            return

        typer.echo(f"📥 Loading {file_name}...")
        mapping_path = None
        if mapping_name:
            mapping_file = directory / mapping_name
            if not auto_mapping and mapping_file.exists():
                mapping_path = str(mapping_file)
                typer.echo(f"📝 Using mapping file: {mapping_file.name}")
            elif auto_mapping:
                typer.echo("📝 Skipping mapping file due to --auto-mapping flag")
            else:
                typer.echo("📝 No mapping file found - will use auto-detection")

        try:
            if loader_fn == load_assays_wrapper:
                loader_fn(file_path=str(file_path), url=url)
            elif is_schema:
                loader_fn(str(file_path), url)
            else:
                loader_fn(
                    csv_file=str(file_path),
                    mapping_file=mapping_path,
                    url=url,
                    error_handling=error_handling,
                    dry_run=dry_run,
                    save_errors=save_errors,
                )
            typer.echo(f"✅ {file_name} loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading {file_name}: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)

    def load_assays_wrapper(file_path, url, **kwargs):
        assays_cli.load_assays(file_path=file_path, url=url)

    tasks = [
        ("compounds_schema.json", add_schema_from_file, None, True),
        ("compounds.csv", compound_cli.load_entity, "compounds_mapping.json", False),
        ("batches_schema.json", add_schema_from_file, None, True),
        ("batches.csv", batch_cli.load_entity, "batches_mapping.json", False),
        ("assays_schema.json", add_schema_from_file, None, True),
        ("assay_runs_schema.json", add_schema_from_file, None, True),
        ("assay_results_schema.json", add_schema_from_file, None, True),
        ("assays.json", load_assays_wrapper, None, False),
        ("assay_runs.csv", assay_runs_cli.load_entity, "assay_runs_mapping.json", False),
        ("assay_results.csv", assay_results_cli.load_entity, "assay_results_mapping.json", False),
    ]

    # Execute all tasks
    for file_name, loader_fn, mapping_name, is_schema in tasks:
        load_file(file_name, loader_fn, mapping_name, is_schema)

    typer.echo("\n✅ Directory loading completed!")
