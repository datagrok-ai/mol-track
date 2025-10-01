import csv
import json
from pathlib import Path

import typer


def validate_file_exists(file_path: str) -> Path:
    """Validate that a file exists and is a file."""
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        typer.echo(f"Error: File '{file_path}' does not exist.", err=True)
        raise typer.Exit(1)

    if not file_path_obj.is_file():
        typer.echo(f"Error: '{file_path}' is not a file.", err=True)
        raise typer.Exit(1)

    return file_path_obj


def load_and_validate_json(file_path: str, model_class=None) -> dict:
    """Load and validate JSON from a file, optionally using a Pydantic model."""
    file_path_obj = validate_file_exists(file_path)

    # Read and parse JSON file
    try:
        with open(file_path_obj, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        typer.echo(f"Error: Invalid JSON format in file '{file_path}': {e}", err=True)
        raise typer.Exit(1)
    except Exception as e:
        typer.echo(f"Error reading file '{file_path}': {e}", err=True)
        raise typer.Exit(1)

    # Validate with Pydantic model if provided
    if model_class:
        if model_class is None:
            typer.echo("⚠️  Warning: Model class is not available for validation. Skipping validation.", err=True)
            return data

        try:
            model_class(**data)
            model_name = model_class.__name__
            typer.echo(f"✅ JSON validation passed using {model_name} model!")
            return data
        except Exception as e:
            model_name = model_class.__name__
            typer.secho(f"❌ JSON validation failed using {model_name} model: {e}", fg=typer.colors.RED, err=True)
            raise typer.Exit(1)
    else:
        typer.echo("✅ JSON loaded successfully (no validation model specified)")

    return data


def load_csv_data(file_path: str, max_rows: int | None = None) -> list[dict[str, str]]:
    """Load CSV data from a file and return as list of dictionaries."""
    file_path_obj = validate_file_exists(file_path)

    try:
        with open(file_path_obj, "r", newline="", encoding="utf-8") as f:
            # Read the first line to get headers and trim them
            header_line = f.readline().strip()
            if not header_line:
                typer.echo("Error: CSV file is empty or has no header row.", err=True)
                raise typer.Exit(1)

            # Parse headers and trim whitespace
            headers = [header.strip() for header in next(csv.reader([header_line]))]

            # Create a reader with the trimmed headers
            reader = csv.DictReader(f, fieldnames=headers)
            data = list(reader)

            if max_rows is not None:
                # Ensure max_rows is an integer
                try:
                    max_rows_int = int(max_rows)
                except (ValueError, TypeError):
                    typer.echo(f"Warning: Invalid max_rows value '{max_rows}', ignoring row limit", err=True)
                    max_rows_int = None

                if max_rows_int is not None:
                    data = data[:max_rows_int]
                    typer.echo(f"📊 Limited to {len(data)} data rows (requested: {max_rows_int})")

            return data
    except Exception as e:
        typer.echo(f"Error reading CSV file '{file_path}': {e}", err=True)
        raise typer.Exit(1)


def validate_and_load_csv_data(
    csv_file: str, entity_type: str, max_rows: int | None = None
) -> tuple[Path, list[dict[str, str]]]:
    """
    Validate CSV file and load data for processing.

    Args:
        csv_file: Path to the CSV file
        entity_type: Type of entity (e.g., "compounds", "batches") for user messages
        max_rows: Maximum number of data rows to load (optional)

    Returns:
        Tuple of (csv_path, csv_data)
    """
    csv_path = validate_file_exists(csv_file)

    typer.echo(f"Loading CSV data from '{csv_file}'...")
    csv_data = load_csv_data(csv_file, max_rows)

    if not csv_data:
        typer.echo("Error: CSV file is empty or has no data rows.", err=True)
        raise typer.Exit(1)

    return csv_path, csv_data


def load_and_validate_mapping(mapping_file: str | None) -> dict[str, any] | None:
    """
    Load and validate mapping file if provided.

    Args:
        mapping_file: Path to the mapping JSON file (optional)

    Returns:
        Mapping data dictionary or None
    """
    mapping_data = None
    if mapping_file:
        typer.echo(f"Loading mapping from '{mapping_file}'...")
        mapping_data = load_and_validate_json(mapping_file)

        if not isinstance(mapping_data, dict):
            typer.echo("Error: Mapping must be a JSON object.", err=True)
            raise typer.Exit(1)

    return mapping_data


def write_result_to_file(data, output_format, output_file, parsed=True):
    """
    Write the response data to a file in the specified format.

    Args:
        response: The response object from the API call
        output_format: The desired output format ("json" or "csv")
        output_file: Path to the output file
    """
    if not output_file:
        return  # No output file specified

    try:
        with open(output_file, "w", encoding="utf-8") as f:
            if output_format == "json" or output_format == "table":
                if not parsed:
                    data = data.json()
                json.dump(data, f, indent=2)
            elif output_format == "csv":
                if not parsed:
                    data = data.text
                f.write(data)
            else:
                typer.secho(
                    f"Warning: Unsupported output format '{output_format}'. No data written.",
                    fg=typer.colors.RED,
                    err=True,
                )
                return
        typer.secho(f"✅ Results written to '{output_file}'")
    except Exception as e:
        typer.secho(f"Error writing results to file '{output_file}': {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)


def load_input_from_file(input_file: str):
    try:
        with open(input_file, "r") as f:
            data = json.load(f)

            output = data.get("output", None)
            filter = data.get("filter", None)
            aggregations = data.get("aggregations", None)
            return output, filter, aggregations

    except Exception as e:
        typer.secho(f"❌ Error parsing file {input_file}: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
