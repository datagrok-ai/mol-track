from pathlib import Path
import typer

from client.cli.assays import create_assay_results_from_csv, create_assay_runs_from_csv, create_assays
from client.cli.batches import create_batches_from_csv
from client.cli.compounds import add_compounds_from_csv
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

    This command will process files in the following order:
    1. compounds_schema.json (if exists) - loads schema definitions
    2. compounds.csv (if exists) - loads compounds with optional mapping
    3. batches_schema.json (if exists) - loads batch schema definitions
    4. batches.csv (if exists) - loads batches with optional mapping
    5. assays_schema.json (if exists) - loads assay schema definitions
    6. assay_runs_schema.json (if exists) - loads assay run schema definitions
    7. assay_results_schema.json (if exists) - loads assay result schema definitions
    8. assays.json (if exists) - loads assays
    9. assay_runs.csv (if exists) - loads assay runs with optional mapping
    10. assay_results.csv (if exists) - loads assay results with optional mapping

    The compounds_mapping.json, batches_mapping.json, assay_runs_mapping.json, and assay_results_mapping.json files will be used as mapping if they exist.
    """
    directory = Path(directory_path)

    if not directory.exists():
        typer.echo(f"❌ Error: Directory '{directory_path}' does not exist.", err=True)
        raise typer.Exit(1)

    if not directory.is_dir():
        typer.echo(f"❌ Error: '{directory_path}' is not a directory.", err=True)
        raise typer.Exit(1)

    typer.echo(f"📁 Loading contents from directory: {directory_path}")
    typer.echo("=" * 60)

    # Step 1: Load schema if compounds_schema.json exists
    schema_file = directory / "compounds_schema.json"
    if schema_file.exists():
        typer.echo(f"📋 Loading schema from {schema_file.name}...")
        try:
            add_schema_from_file(str(schema_file), url)
            typer.echo("✅ Schema loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading schema: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No compounds_schema.json found - skipping schema loading")

    # Step 2: Load compounds if compounds.csv exists
    compounds_file = directory / "compounds.csv"
    if compounds_file.exists():
        typer.echo(f"🧪 Loading compounds from {compounds_file.name}...")

        # Check for mapping file (unless auto_mapping is specified)
        mapping_file = directory / "compounds_mapping.json"
        mapping_path = None
        if not auto_mapping and mapping_file.exists():
            mapping_path = str(mapping_file)
            typer.echo(f"📝 Using mapping file: {mapping_file.name}")
        else:
            if auto_mapping:
                typer.echo("📝 Skipping mapping file due to --auto-mapping flag")
            else:
                typer.echo("📝 No mapping file found - will use auto-detection")

        try:
            add_compounds_from_csv(
                csv_file=str(compounds_file),
                mapping_file=mapping_path,
                url=url,
                error_handling=error_handling,
                output_format="json",
                dry_run=dry_run,
                save_errors=False,
            )
            typer.echo("✅ Compounds loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading compounds: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No compounds.csv found - skipping compounds loading")

    # Step 3: Load batch schema if batches_schema.json exists
    batch_schema_file = directory / "batches_schema.json"
    if batch_schema_file.exists():
        typer.echo(f"📋 Loading batch schema from {batch_schema_file.name}...")
        try:
            add_schema_from_file(str(batch_schema_file), url)
            typer.echo("✅ Batch schema loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading batch schema: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No batches_schema.json found - skipping batch schema loading")

    # Step 4: Load batches if batches.csv exists
    batches_file = directory / "batches.csv"
    if batches_file.exists():
        typer.echo(f"📦 Loading batches from {batches_file.name}...")

        # Check for mapping file (unless auto_mapping is specified)
        batch_mapping_file = directory / "batches_mapping.json"
        batch_mapping_path = None
        if not auto_mapping and batch_mapping_file.exists():
            batch_mapping_path = str(batch_mapping_file)
            typer.echo(f"📝 Using batch mapping file: {batch_mapping_file.name}")
        else:
            if auto_mapping:
                typer.echo("📝 Skipping mapping file due to --auto-mapping flag")
            else:
                typer.echo("📝 No batch mapping file found - will use auto-detection")

        try:
            create_batches_from_csv(
                csv_file=str(batches_file),
                mapping_file=batch_mapping_path,
                url=url,
                error_handling=error_handling,
                output_format="json",
                dry_run=dry_run,
                save_errors=False,
            )
            typer.echo("✅ Batches loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading batches: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No batches.csv found - skipping batches loading")

    # Step 5: Load assay schema if assays_schema.json exists
    assay_schema_file = directory / "assays_schema.json"
    if assay_schema_file.exists():
        typer.echo(f"📋 Loading assay schema from {assay_schema_file.name}...")
        try:
            add_schema_from_file(str(assay_schema_file), url)
            typer.echo("✅ Assay schema loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading assay schema: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No assays_schema.json found - skipping assay schema loading")

    # Step 6: Load assay run schema if assay_runs_schema.json exists
    assay_run_schema_file = directory / "assay_runs_schema.json"
    if assay_run_schema_file.exists():
        typer.echo(f"📋 Loading assay run schema from {assay_run_schema_file.name}...")
        try:
            add_schema_from_file(str(assay_run_schema_file), url)
            typer.echo("✅ Assay run schema loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading assay run schema: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No assay_runs_schema.json found - skipping assay run schema loading")

    # Step 7: Load assay result schema if assay_results_schema.json exists
    assay_result_schema_file = directory / "assay_results_schema.json"
    if assay_result_schema_file.exists():
        typer.echo(f"📋 Loading assay result schema from {assay_result_schema_file.name}...")
        try:
            add_schema_from_file(str(assay_result_schema_file), url)
            typer.echo("✅ Assay result schema loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading assay result schema: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No assay_results_schema.json found - skipping assay result schema loading")

    # Step 8: Load assays if assays.json exists
    assays_file = directory / "assays.json"
    if assays_file.exists():
        typer.echo(f"🧬 Loading assays from {assays_file.name}...")
        try:
            create_assays(file_path=str(assays_file), url=url)
            typer.echo("✅ Assays loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading assays: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No assays.json found - skipping assays loading")

    # Step 9: Load assay runs if assay_runs.csv exists
    assay_runs_file = directory / "assay_runs.csv"
    if assay_runs_file.exists():
        typer.echo(f"🧪 Loading assay runs from {assay_runs_file.name}...")

        # Check for mapping file (unless auto_mapping is specified)
        assay_runs_mapping_file = directory / "assay_runs_mapping.json"
        assay_runs_mapping_path = None
        if not auto_mapping and assay_runs_mapping_file.exists():
            assay_runs_mapping_path = str(assay_runs_mapping_file)
            typer.echo(f"📝 Using assay runs mapping file: {assay_runs_mapping_file.name}")
        else:
            if auto_mapping:
                typer.echo("📝 Skipping mapping file due to --auto-mapping flag")
            else:
                typer.echo("📝 No assay runs mapping file found - will use auto-detection")

        try:
            create_assay_runs_from_csv(
                csv_file=str(assay_runs_file),
                mapping_file=assay_runs_mapping_path,
                url=url,
                error_handling=error_handling,
                output_format="json",
                dry_run=dry_run,
                save_errors=save_errors,
            )
            typer.echo("✅ Assay runs loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading assay runs: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No assay_runs.csv found - skipping assay runs loading")

    # Step 10: Load assay results if assay_results.csv exists
    assay_results_file = directory / "assay_results.csv"
    if assay_results_file.exists():
        typer.echo(f"📊 Loading assay results from {assay_results_file.name}...")

        # Check for mapping file (unless auto_mapping is specified)
        assay_results_mapping_file = directory / "assay_results_mapping.json"
        assay_results_mapping_path = None
        if not auto_mapping and assay_results_mapping_file.exists():
            assay_results_mapping_path = str(assay_results_mapping_file)
            typer.echo(f"📝 Using assay results mapping file: {assay_results_mapping_file.name}")
        else:
            if auto_mapping:
                typer.echo("📝 Skipping mapping file due to --auto-mapping flag")
            else:
                typer.echo("📝 No assay results mapping file found - will use auto-detection")

        try:
            create_assay_results_from_csv(
                csv_file=str(assay_results_file),
                mapping_file=assay_results_mapping_path,
                url=url,
                error_handling=error_handling,
                output_format="json",
                dry_run=dry_run,
                save_errors=save_errors,
            )
            typer.echo("✅ Assay results loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading assay results: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No assay_results.csv found - skipping assay results loading")

    typer.echo("\n✅ Directory loading completed!")
