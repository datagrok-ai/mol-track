# client.py
import sys
import csv
import json
from pathlib import Path
import typer
import requests
from rich.console import Console
from rich.table import Table
from datetime import datetime
from sqlalchemy import text
from dateutil.parser import parse as date_parse

# Add parent directory to Python path for imports
# sys.path.append("..")
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

try:
    # from app.models import SchemaPayload
    # from .database import engine, DB_SCHEMA
    from app.models import SchemaPayload
except ImportError:
    # Fallback for when models is not in the path
    # from app.models import SchemaPayload
    print("..app.models import failed")


try:
    from app.setup.database import engine, DB_SCHEMA
except ImportError:
    try:
        from app.database import engine, DB_SCHEMA
    except ImportError:
        # Set to None for when running standalone
        engine = None
        DB_SCHEMA = None

app = typer.Typer()
schema_app = typer.Typer()
compound_app = typer.Typer()
batch_app = typer.Typer()
batch_list_app = typer.Typer()
properties_app = typer.Typer()
additions_app = typer.Typer()
additions_list_app = typer.Typer()
assays_app = typer.Typer()
assays_runs_app = typer.Typer()
assays_results_app = typer.Typer()
database_app = typer.Typer()
directory_app = typer.Typer()
csv_app = typer.Typer()
search_app = typer.Typer()
admin_app = typer.Typer()
app.add_typer(schema_app, name="schema", help="Schema management commands")
app.add_typer(compound_app, name="compounds", help="Compound management commands")
app.add_typer(batch_app, name="batches", help="Batch management commands")
app.add_typer(properties_app, name="properties", help="Property management commands")
app.add_typer(additions_app, name="additions", help="Addition management commands")
app.add_typer(assays_app, name="assays", help="Assays management commands")
app.add_typer(database_app, name="database", help="Database management commands")
app.add_typer(directory_app, name="directory", help="Directory loading commands")
app.add_typer(csv_app, name="csv", help="CSV utility commands")
app.add_typer(search_app, name="search", help="Search functionality")
app.add_typer(admin_app, name="admin", help="Administrative functions")
batch_app.add_typer(batch_list_app, name="list", help="List batch information")
additions_app.add_typer(additions_list_app, name="list", help="List additions information")
assays_app.add_typer(assays_runs_app, name="runs", help="Assay runs management")
assays_app.add_typer(assays_results_app, name="results", help="Assay results management")

DEFAULT_SERVER_URL = "http://127.0.0.1:8000"


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
            validated_data = model_class(**data)
            model_name = model_class.__name__
            typer.echo(f"✅ JSON validation passed using {model_name} model!")
            return validated_data.model_dump()
        except Exception as e:
            model_name = model_class.__name__
            typer.echo(f"❌ JSON validation failed using {model_name} model: {e}", err=True)
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


def map_csv_row_to_compound(row: dict[str, str], mapping: dict[str, str]) -> dict:
    """Map a CSV row to compound data using the provided mapping."""
    compound_data = {}

    for csv_field, model_field in mapping.items():
        if csv_field in row:
            value = row[csv_field].strip() if row[csv_field] else None
            if value is not None:
                compound_data[model_field] = value

    return compound_data


def print_response(response):
    if response.status_code == 200:
        print(response.json())
    else:
        print(f"Error: {response.status_code}")


# New utility functions for CSV processing
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


def report_csv_information(csv_data: list[dict[str, str]], mapping_data: dict[str, any] | None, entity_type: str):
    """
    Report information about CSV data and mapping.

    Args:
        csv_data: Loaded CSV data
        mapping_data: Optional mapping data
        entity_type: Type of entity for user messages
    """
    csv_headers = set(csv_data[0].keys())
    typer.echo(f"✅ Found {len(csv_data)} {entity_type} in CSV file")
    typer.echo(f"📋 CSV headers: {csv_headers}")

    if mapping_data:
        mapped_headers = set(mapping_data.keys())
        unmapped_headers = csv_headers - mapped_headers
        typer.echo(f"📋 Mapped columns: {len(mapped_headers)}")
        if unmapped_headers:
            typer.echo(f"⚠️  Unmapped columns (will be ignored): {unmapped_headers}")


def send_csv_upload_request(
    csv_path: Path,
    mapping_data: dict[str, any] | None,
    url: str,
    endpoint: str,
    error_handling: str,
    output_format: str,
    entity_type: str,
    csv_data: list[dict[str, str]] | None = None,
    save_errors: bool = False,
) -> None:
    """
    Send CSV upload request to the specified endpoint.

    Args:
        csv_path: Path to the CSV file
        mapping_data: Optional mapping data
        url: Server URL
        endpoint: API endpoint to send request to
        error_handling: Error handling strategy
        output_format: Output format
        entity_type: Type of entity for user messages
        csv_data: Optional CSV data to use instead of reading from file
        save_errors: Whether to save error records to a file
    """
    try:
        if csv_data is not None:
            # Create temporary CSV file with limited data
            import tempfile

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".csv", delete=False, newline="", encoding="utf-8"
            ) as temp_file:
                if csv_data:
                    writer = csv.DictWriter(temp_file, fieldnames=csv_data[0].keys())
                    writer.writeheader()
                    writer.writerows(csv_data)
                    temp_file_path = temp_file.name
                else:
                    temp_file_path = None
        else:
            temp_file_path = None

        try:
            # Use temporary file if created, otherwise use original
            file_to_send = temp_file_path if temp_file_path else csv_path

            with open(file_to_send, "rb") as f:
                files = {"csv_file": (csv_path.name, f, "text/csv")}

                data = {"error_handling": error_handling, "output_format": output_format}

                if mapping_data:
                    data["mapping"] = json.dumps(mapping_data)

                typer.echo(f"🚀 Sending {csv_path.name} to {url}{endpoint}...")

                response = requests.post(f"{url}{endpoint}", files=files, data=data)

                if response.status_code == 200:
                    result = response.json()
                    typer.echo(f"✅ {entity_type.capitalize()} registered successfully!")

                    # Parse the result based on output format
                    if output_format == "json":
                        if "data" in result:
                            data_list = result["data"]
                            success_count = sum(1 for item in data_list if item.get("registration_status") == "success")
                            error_count = len(data_list) - success_count
                            typer.echo(f"📊 Results: {success_count} successful, {error_count} errors")

                            # Show any errors
                            errors = [item for item in data_list if item.get("registration_status") != "success"]
                            if errors:
                                typer.echo("❌ Errors found:")
                                for error in errors[:5]:  # Show first 5 errors
                                    typer.echo(
                                        f"  - Row {error.get('row', 'unknown')}: {error.get('registration_error_message', 'Unknown error')}"
                                    )
                                if len(errors) > 5:
                                    typer.echo(f"  ... and {len(errors) - 5} more errors")
                                
                                # Save errors to file if requested
                                if save_errors:
                                    # Generate filename based on endpoint
                                    endpoint_name = endpoint.strip("/").replace("/", "_")
                                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                                    error_filename = f"{endpoint_name}_errors_{timestamp}.json"
                                    
                                    try:
                                        with open(error_filename, "w") as error_file:
                                            json.dump(errors, error_file, indent=2)
                                        typer.echo(f"💾 Error records saved to: {error_filename}")
                                    except Exception as e:
                                        typer.echo(f"⚠️  Warning: Could not save error file: {e}")
                        else:
                            typer.echo("✅ Registration completed")
                    else:
                        typer.echo("✅ CSV output received")

                else:
                    typer.echo(f"❌ Error: {response.status_code}")
                    try:
                        error_detail = response.json()
                        typer.echo(f"Details: {json.dumps(error_detail, indent=2)}")
                    except (json.JSONDecodeError, ValueError):
                        typer.echo(f"Response: {response.text}")
        finally:
            # Clean up temporary file
            if temp_file_path and temp_file_path != str(csv_path):
                try:
                    import os

                    os.unlink(temp_file_path)
                except Exception:
                    pass  # Ignore cleanup errors

    except requests.exceptions.ConnectionError:
        typer.echo(f"❌ Error: Could not connect to server at {url}", err=True)
        raise typer.Exit(1)
    except requests.exceptions.RequestException as e:
        typer.echo(f"❌ Error making request: {e}", err=True)
        raise typer.Exit(1)
    except Exception as e:
        typer.echo(f"❌ Error: {e}", err=True)
        raise typer.Exit(1)


# Table display utility functions
def format_timestamp(timestamp_str: str) -> str:
    """
    Format a timestamp string to a readable format.

    Args:
        timestamp_str: ISO format timestamp string

    Returns:
        Formatted timestamp string or original if parsing fails
    """
    if not timestamp_str:
        return ""

    try:
        # Parse and format the timestamp
        dt = datetime.fromisoformat(timestamp_str.replace("Z", "+00:00"))
        return dt.strftime("%Y-%m-%d %H:%M:%S")
    except (ValueError, AttributeError):
        # If parsing fails, use the original value
        return timestamp_str


def create_rich_table(title: str, columns: list[tuple[str, str, dict]]) -> Table:
    """
    Create a rich table with the specified columns.

    Args:
        title: Table title
        columns: List of (header, style, kwargs) tuples for each column

    Returns:
        Configured Table object
    """
    table = Table(title=title)

    for header, style, kwargs in columns:
        table.add_column(header, style=style, **kwargs)

    return table


def display_data_table(
    data: list[dict],
    title: str,
    columns: list[tuple[str, str, dict]],
    row_extractor: callable,
    show_total: bool = False,
    total_label: str = "TOTAL",
) -> None:
    """
    Display data in a rich table format.

    Args:
        data: List of data dictionaries
        title: Table title
        columns: List of (header, style, kwargs) tuples for each column
        row_extractor: Function that extracts row values from a data item
        show_total: Whether to show a total row
        total_label: Label for the total row
    """
    console = Console()
    table = create_rich_table(title, columns)

    total_value = 0

    for item in data:
        row_values = row_extractor(item)
        table.add_row(*row_values)

        # Update total if needed and if the last value is numeric
        if show_total and len(row_values) > 0:
            try:
                total_value += float(str(row_values[-1]).replace(",", ""))
            except (ValueError, TypeError):
                pass

    # Add total row if requested
    if show_total:
        table.add_row("", "")
        table.add_row(total_label, f"{total_value:,.0f}", style="bold")

    console.print(table)


def get_table_row_counts(specific_tables: list[str] | None = None) -> dict[str, int]:
    """
    Get row counts for database tables.

    Args:
        specific_tables: List of specific table names to count. If None, gets all tables.

    Returns:
        Dictionary mapping table names to row counts
    """
    # Use the already-imported engine and DB_SCHEMA
    if engine is None or DB_SCHEMA is None:
        raise ImportError("Database connection not available - engine or DB_SCHEMA is None")

    row_counts = {}

    with engine.connect() as connection:
        # Get all tables in the schema
        tables_query = text("""
            SELECT tablename 
            FROM pg_tables 
            WHERE schemaname = :schema
            ORDER BY tablename
        """)

        tables_result = connection.execute(tables_query, {"schema": DB_SCHEMA})
        all_tables = [row[0] for row in tables_result]

        # Filter tables if specific_tables is provided
        if specific_tables is not None:
            all_tables = [table for table in all_tables if table in specific_tables]

        # TODO build a union query so that theree is one traversal of the database
        # Get actual row counts for each table
        for table in all_tables:
            count_query = text(f"SELECT COUNT(*) FROM {table}")
            count_result = connection.execute(count_query)
            row_counts[table] = count_result.scalar() or 0

    return row_counts


# Schema Commands
@schema_app.command("list")
def list_schema_group(url: str = DEFAULT_SERVER_URL):
    """
    List the schema.
    """
    print(f"Listing schema from {url}")
    response = requests.get(f"{url}/v1/schema/compounds")
    typer.echo(f"Compounds schema:\n{json.dumps(response.json(), indent=2)}")
    response = requests.get(f"{url}/v1/schema/batches")
    typer.echo(f"Batches schema:\n{json.dumps(response.json(), indent=2)}")


@schema_app.command("load")
def add_schema_from_file(
    file_path: str = typer.Argument(..., help="Path to the JSON file containing schema data"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Add schema from a JSON file.

    The file should contain a JSON object with 'properties' and optionally 'synonym_types' arrays.
    Example format:
    {
        "properties": [
            {
                "name": "Molecular Weight",
                "value_type": "number",
                "property_class": "CALCULATED",
                "unit": "g/mol",
                "scope": "COMPOUND",
                "description": "Molecular weight of the compound"
            }
        ],
        "synonym_types": [
            {
                "name": "CAS Number",
                "value_type": "string",
                "property_class": "DECLARED",
                "unit": "",
                "scope": "COMPOUND",
                "pattern": "^\\d{1,7}-\\d{2}-\\d$",
                "description": "CAS Registry Number",
                "semantic_type_id": 1
            }
        ]
    }
    """
    # Load and validate schema using utility function
    schema_data = load_and_validate_json(file_path, SchemaPayload)

    # Send request to server
    try:
        typer.echo(f"Adding schema from file '{file_path}' to {url}...")
        response = requests.post(f"{url}/v1/schema", json=schema_data)

        if response.status_code == 200:
            result = response.json()
            typer.echo("✅ Schema added successfully!")

            # Report detailed statistics
            if "created" in result:
                created = result["created"]

                # Properties reporting
                if "properties" in created:
                    properties_created = len(created["properties"]) if created["properties"] else 0
                    typer.echo(f"📋 Properties created: {properties_created}")

                # Synonym types reporting
                if "synonym_types" in created:
                    synonyms_created = len(created["synonym_types"]) if created["synonym_types"] else 0
                    typer.echo(f"🏷️  Synonym types created: {synonyms_created}")

            # Report skipped items if available
            if "skipped" in result:
                skipped = result["skipped"]

                if "properties" in skipped:
                    properties_skipped = len(skipped["properties"]) if skipped["properties"] else 0
                    if properties_skipped > 0:
                        typer.echo(f"⏭️  Properties skipped: {properties_skipped}")

                if "synonym_types" in skipped:
                    synonyms_skipped = len(skipped["synonym_types"]) if skipped["synonym_types"] else 0
                    if synonyms_skipped > 0:
                        typer.echo(f"⏭️  Synonym types skipped: {synonyms_skipped}")

            # Report total counts from input file
            input_properties = len(schema_data.get("properties", []))
            input_synonyms = len(schema_data.get("synonym_types", []))

            if input_properties > 0 or input_synonyms > 0:
                typer.echo(f"📊 Summary: {input_properties} properties and {input_synonyms} synonym types processed")
        else:
            typer.echo(f"❌ Error: {response.status_code}")
            try:
                error_detail = response.json()
                typer.echo(f"Details: {json.dumps(error_detail, indent=2)}")
            except (json.JSONDecodeError, ValueError):
                typer.echo(f"Response: {response.text}")

    except requests.exceptions.ConnectionError:
        typer.echo(f"❌ Error: Could not connect to server at {url}", err=True)
        raise typer.Exit(1)
    except requests.exceptions.RequestException as e:
        typer.echo(f"❌ Error making request: {e}", err=True)
        raise typer.Exit(1)


# Compound Commands
@compound_app.command("list")
def list_compounds_group(
    skip: int = 0,
    limit: int = 10,
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    List compounds using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/compounds/?skip={skip}&limit={limit}")

    if response.status_code == 200:
        compounds_data = response.json()

        if output_format == "json":
            print(json.dumps(compounds_data, indent=2))
        else:
            # Default to table format
            display_compounds_table(compounds_data)
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


def display_compounds_table(compounds_data):
    """Display compounds data in a rich table format."""
    columns = [
        ("ID", "cyan", {"no_wrap": True}),
        ("MolRegNo", "blue", {"no_wrap": True}),
        ("Common Name", "green", {"no_wrap": True}),
        ("CAS", "yellow", {"no_wrap": True}),
        ("SMILES", "magenta", {"no_wrap": True}),
        ("Synonyms", "white", {"no_wrap": True}),
        ("Created At", "red", {}),
    ]

    def extract_compound_row(compound):
        # Extract common name, CAS, and synonyms from properties
        common_name = ""
        cas_number = ""
        synonyms = []

        for prop in compound.get("properties", []):
            prop_name = prop.get("name", "").lower()
            prop_value = prop.get("value_string", "")
            semantic_type_id = prop.get("semantic_type_id")

            if prop_name == "common name":
                common_name = prop_value
            elif prop_name == "cas":
                cas_number = prop_value
            elif semantic_type_id == 1 and prop_value:  # Synonyms have semantic_type_id=1
                synonyms.append(f"{prop.get('name', '')}: {prop_value}")

        # Join synonyms with semicolon
        synonyms_str = "; ".join(synonyms[:3])  # Limit to first 3 synonyms
        if len(synonyms) > 3:
            synonyms_str += f" (+{len(synonyms) - 3} more)"

        # Truncate long fields for better display
        smiles = compound.get("canonical_smiles", "")
        if len(smiles) > 35:
            smiles = smiles[:32] + "..."

        if len(common_name) > 25:
            common_name = common_name[:22] + "..."

        if len(synonyms_str) > 30:
            synonyms_str = synonyms_str[:27] + "..."

        # Get molregno if available
        molregno = compound.get("molregno", "")

        return [
            str(compound.get("id", "")),
            str(molregno) if molregno else "",
            common_name,
            cas_number,
            smiles,
            synonyms_str,
            format_timestamp(compound.get("created_at", "")),
        ]

    display_data_table(data=compounds_data, title="Compounds", columns=columns, row_extractor=extract_compound_row)


def display_batches_table(batches_data):
    """Display batches data in a rich table format."""
    columns = [
        ("ID", "cyan", {"no_wrap": True}),
        ("Batch Regno", "blue", {"no_wrap": True}),
        ("Compound ID", "green", {"no_wrap": True}),
        ("Notes", "yellow", {"no_wrap": True}),
        ("Created At", "red", {}),
    ]

    def extract_batch_row(batch):
        return [
            str(batch.get("id", "")),
            str(batch.get("batch_regno", "")),
            str(batch.get("compound_id", "")),
            batch.get("notes", "")[:30] + "..." if batch.get("notes") and len(batch.get("notes", "")) > 30 else batch.get("notes", ""),
            format_timestamp(batch.get("created_at", "")),
        ]

    display_data_table(data=batches_data, title="Batches", columns=columns, row_extractor=extract_batch_row)


def display_assays_table(assays_data):
    """Display assays data in a rich table format."""
    columns = [
        ("ID", "cyan", {"no_wrap": True}),
        ("Name", "blue", {"no_wrap": True}),
        ("Description", "green", {"no_wrap": True}),
        ("Created At", "red", {}),
    ]

    def extract_assay_row(assay):
        return [
            str(assay.get("id", "")),
            assay.get("name", ""),
            assay.get("description", "")[:30] + "..." if assay.get("description") and len(assay.get("description", "")) > 30 else assay.get("description", ""),
            format_timestamp(assay.get("created_at", "")),
        ]

    display_data_table(data=assays_data, title="Assays", columns=columns, row_extractor=extract_assay_row)


@compound_app.command("get")
def get_compound(
    compound_id: int = typer.Argument(..., help="Compound ID to retrieve"),
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get a specific compound by ID.
    """
    response = requests.get(f"{url}/v1/compounds/{compound_id}")

    if response.status_code == 200:
        compound_data = response.json()

        if output_format == "json":
            print(json.dumps(compound_data, indent=2))
        else:
            # Display single compound in table format
            display_compounds_table([compound_data])
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


@compound_app.command("properties")
def get_compound_properties(
    compound_id: int = typer.Argument(..., help="Compound ID to get properties for"),
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get all properties for a specific compound.
    """
    response = requests.get(f"{url}/v1/compounds/{compound_id}/properties")

    if response.status_code == 200:
        properties_data = response.json()

        if output_format == "json":
            print(json.dumps(properties_data, indent=2))
        else:
            # Display properties in table format
            display_properties_table(properties_data)
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


@compound_app.command("synonyms")
def get_compound_synonyms(
    compound_id: int = typer.Argument(..., help="Compound ID to get synonyms for"),
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get all synonyms for a specific compound.
    """
    response = requests.get(f"{url}/v1/compounds/{compound_id}/synonyms")

    if response.status_code == 200:
        synonyms_data = response.json()

        if output_format == "json":
            print(json.dumps(synonyms_data, indent=2))
        else:
            # Display synonyms in table format
            display_properties_table(synonyms_data)
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


@compound_app.command("delete")
def delete_compound(
    compound_id: int = typer.Argument(..., help="Compound ID to delete"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Delete a specific compound by ID.
    """
    response = requests.delete(f"{url}/v1/compounds/{compound_id}")
    print_response(response)


@compound_app.command("load")
def add_compounds_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing compound data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
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


# List Commands
@properties_app.command("list")
def list_properties(
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    List the properties.
    """
    response = requests.get(f"{url}/properties")

    if response.status_code == 200:
        properties_data = response.json()

        if output_format == "json":
            print(json.dumps(properties_data, indent=2))
        else:
            # Default to table format
            display_properties_table(properties_data)
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


def display_properties_table(properties_data):
    """Display properties data in a rich table format."""
    columns = [
        ("Name", "cyan", {"no_wrap": True}),
        ("Scope", "magenta", {}),
        ("Value Type", "green", {}),
        ("Semantic Type ID", "yellow", {}),
        ("Created At", "blue", {}),
    ]

    def extract_property_row(prop):
        return [
            prop.get("name", ""),
            prop.get("scope", ""),
            prop.get("value_type", ""),
            str(prop.get("semantic_type_id", "")),
            format_timestamp(prop.get("created_at", "")),
        ]

    display_data_table(data=properties_data, title="Properties", columns=columns, row_extractor=extract_property_row)


# Batch Commands
@batch_app.command("list")
def list_batches_group(
    batch_id: int | None = typer.Argument(None, help="Batch ID to retrieve (optional)"),
    skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
    limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    List batches using the v1 endpoint.

    If no batch_id is provided, lists all batches with pagination.
    If batch_id is provided, gets the specific batch.
    """
    if batch_id is not None:
        # Get specific batch
        response = requests.get(f"{url}/v1/batches/{batch_id}")
    else:
        # List all batches
        response = requests.get(f"{url}/v1/batches/?skip={skip}&limit={limit}")
    print_response(response)


@batch_app.command("load")
def create_batches_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing batch data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json or csv"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
    save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
):
    """
    Add batches from a CSV file using the /v1/batches/ endpoint.

    The CSV file should contain batch data with headers. If a mapping file is provided,
    it should be a JSON object mapping CSV column names to field names.

    Example mapping:
    {
        "structure": "smiles",
        "batch_corporate_id": "batches_details.batch_corporate_id",
        "common_name": "compounds_details.common_name",
        "cas": "compounds_details.cas"
    }

    If no mapping is provided, the system will attempt to infer the mapping from column names.

    Use --rows to limit the number of data rows processed (header row is excluded).
    """

    # Use utility functions for common steps
    csv_path, csv_data = validate_and_load_csv_data(csv_file, "batches", rows)
    mapping_data = load_and_validate_mapping(mapping_file)
    report_csv_information(csv_data, mapping_data, "batches")

    if dry_run:
        typer.echo("✅ Dry run completed successfully!")
        return

    # Send the request using utility function
    send_csv_upload_request(
        csv_path=csv_path,
        mapping_data=mapping_data,
        url=url,
        endpoint="/v1/batches/",
        error_handling=error_handling,
        output_format=output_format,
        entity_type="batches",
        csv_data=csv_data if rows is not None else None,
        save_errors=save_errors,
    )


@batch_app.command("get")
def get_batch(
    batch_id: int = typer.Argument(..., help="Batch ID to retrieve"),
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Get a specific batch by ID.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}")

    if response.status_code == 200:
        batch_data = response.json()

        if output_format == "json":
            print(json.dumps(batch_data, indent=2))
        else:
            # Display single batch in table format
            display_batches_table([batch_data])
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


@batch_app.command("delete")
def delete_batch(
    batch_id: int = typer.Argument(..., help="Batch ID to delete"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Delete a specific batch by ID.
    """
    response = requests.delete(f"{url}/v1/batches/{batch_id}")
    print_response(response)


# Batch List Commands
@batch_list_app.command("properties")
def get_batch_properties(
    batch_id: int = typer.Argument(..., help="Batch ID to get properties for"), url: str = DEFAULT_SERVER_URL
):
    """
    Get all properties for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/properties")
    print_response(response)


@batch_list_app.command("synonyms")
def get_batch_synonyms(
    batch_id: int = typer.Argument(..., help="Batch ID to get synonyms for"), url: str = DEFAULT_SERVER_URL
):
    """
    Get all synonyms for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/synonyms")
    print_response(response)


@batch_list_app.command("additions")
def get_batch_additions(
    batch_id: int = typer.Argument(..., help="Batch ID to get additions for"), url: str = DEFAULT_SERVER_URL
):
    """
    Get all additions for a specific batch using the v1 endpoint.
    """
    response = requests.get(f"{url}/v1/batches/{batch_id}/additions")
    print_response(response)


# Additions Commands
@additions_app.command("list")
def list_additions(url: str = DEFAULT_SERVER_URL):
    """
    List the additions.
    """
    response = requests.get(f"{url}/v1/additions")
    print_response(response)


@additions_app.command("update")
def update_addition(
    addition_id: int = typer.Argument(..., help="Addition ID to update"),
    file_path: str = typer.Argument(..., help="Path to the JSON file containing addition update data"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Update information for the specified addition.
    """
    # Load and validate the update data
    update_data = load_and_validate_json(file_path)

    response = requests.put(f"{url}/v1/additions/{addition_id}", json=update_data)
    print_response(response)


@additions_app.command("delete")
def delete_addition(
    addition_id: int = typer.Argument(..., help="Addition ID to delete"), url: str = DEFAULT_SERVER_URL
):
    """
    Soft delete the specified addition (only if no dependent batches exist).
    """
    response = requests.delete(f"{url}/v1/additions/{addition_id}")
    print_response(response)


@additions_app.command("load")
def add_additions_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing addition data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json or csv"),
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
        output_format=output_format,
        entity_type="additions",
        csv_data=csv_data if rows is not None else None,
        save_errors=save_errors,
    )


# Additions List Commands
@additions_list_app.command("salts")
def list_additions_salts(url: str = DEFAULT_SERVER_URL):
    """
    List all additions with role of salts.
    """
    response = requests.get(f"{url}/v1/additions/salts")
    print_response(response)


@additions_list_app.command("solvates")
def list_additions_solvates(url: str = DEFAULT_SERVER_URL):
    """
    List all additions with role of solvates.
    """
    response = requests.get(f"{url}/v1/additions/solvates")
    print_response(response)


@additions_list_app.command()
def get_addition(addition_id: int = typer.Argument(..., help="Addition ID to retrieve"), url: str = DEFAULT_SERVER_URL):
    """
    Get all information for a specific addition.
    """
    response = requests.get(f"{url}/v1/additions/{addition_id}")
    print_response(response)


# Assays Commands
@assays_app.command("list")
def list_assays(
    assay_id: int | None = typer.Argument(None, help="Assay ID to retrieve (optional)"), url: str = DEFAULT_SERVER_URL
):
    """
    List assays using the v1 endpoint.

    If no assay_id is provided, lists all assays.
    If assay_id is provided, gets the specific assay.
    """
    if assay_id is not None:
        # Get specific assay
        response = requests.get(f"{url}/v1/assays/{assay_id}")
    else:
        # List all assays
        response = requests.get(f"{url}/v1/assays")
    print_response(response)


@assays_app.command("load")
def create_assays(
    file_path: str = typer.Argument(..., help="Path to the JSON file containing assay data"),
    url: str = DEFAULT_SERVER_URL,
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
    url: str = DEFAULT_SERVER_URL,
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
    assay_run_id: int | None = typer.Argument(None, help="Assay run ID to retrieve (optional)"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    List assay runs using the v1 endpoint.

    If no assay_run_id is provided, lists all assay runs.
    If assay_run_id is provided, gets the specific assay run.
    """
    if assay_run_id is not None:
        # Get specific assay run
        response = requests.get(f"{url}/v1/assay_runs/{assay_run_id}")
    else:
        # List all assay runs
        response = requests.get(f"{url}/v1/assay_runs")
    print_response(response)


@assays_runs_app.command("load")
def create_assay_runs_from_csv(
    csv_file: str = typer.Argument(..., help="Path to the CSV file containing assay run data"),
    mapping_file: str = typer.Option(None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"),
    rows: int | None = typer.Option(None, "--rows", "-r", help="Number of data rows to process (excludes header row)"),
    url: str = DEFAULT_SERVER_URL,
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
    url: str = DEFAULT_SERVER_URL,
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
    url: str = DEFAULT_SERVER_URL,
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


# Database Commands
@database_app.command("stats")
def database_stats():
    """
    Show database table statistics with row counts.
    """
    try:
        # Get row counts for all tables using shared utility
        row_counts = get_table_row_counts()

        # Convert to the format expected by display function
        table_stats = [(DB_SCHEMA, table_name, count) for table_name, count in row_counts.items()]

        # Display results in a table
        display_database_stats_table(table_stats, DB_SCHEMA)

    except Exception as e:
        typer.echo(f"❌ Error connecting to database: {e}", err=True)
        raise typer.Exit(1)


@database_app.command("clean")
def database_clean(
    force: bool = typer.Option(False, "--force", "-f", help="Skip confirmation prompt"), url: str = DEFAULT_SERVER_URL
):
    """
    Clean the database by deleting all data rows in dependency order.

    This will delete rows from the following tables in order:
    - assay_results
    - assay_run_details
    - assay_runs
    - assay_properties
    - assay_details
    - assays
    - batch_additions
    - batch_details
    - batches
    - additions
    - compound_details
    - compounds
    - properties (except corporate_compound_id and corporate_batch_id properties)

    The following tables are preserved:
    - semantic_types
    - settings
    - users
    """
    try:
        # Tables to clean in dependency order (child tables first)
        tables_to_clean = [
            "assay_results",
            "assay_run_details",
            "assay_runs",
            "assay_properties",
            "assay_details",
            "assays",
            "batch_additions",
            "batch_details",
            "batches",
            "additions",
            "compound_details",
            "compounds",
            "properties",
        ]

        # Tables to preserve
        tables_to_preserve = ["semantic_types", "settings", "users"]

        # Show what will be cleaned
        typer.echo("🗑️  Database Clean Operation")
        typer.echo("=" * 50)
        typer.echo("Tables to be cleaned (in dependency order):")
        for i, table in enumerate(tables_to_clean, 1):
            typer.echo(f"  {i:2d}. {table}")

        typer.echo("\nTables to be preserved:")
        for table in tables_to_preserve:
            typer.echo(f"  ✓ {table}")

        # typer.echo(f"\nDatabase: {DB_NAME} (schema: {DB_SCHEMA})")
        # typer.echo(f"Host: {DB_HOST}:{DB_PORT}")

        # Get row counts before deletion
        typer.echo("\n📊 Getting current row counts...")
        row_counts = get_table_row_counts(tables_to_clean)

        # Display current row counts in a rich table
        typer.echo("\nCurrent row counts:")
        table_data = [
            {"table": table, "count": row_counts[table]} for table in tables_to_clean if row_counts[table] > 0
        ]

        if table_data:
            # TODO: look for optimizations
            # Create rich table for row counts
            columns = [("Table Name", "cyan", {"no_wrap": True}), ("Row Count", "red", {"justify": "right"})]

            def extract_count_row(item):
                return [item["table"], f"{item['count']:,}"]

            display_data_table(
                data=table_data,
                title="Tables to be cleaned",
                columns=columns,
                row_extractor=extract_count_row,
                show_total=True,
                total_label="TOTAL ROWS TO DELETE",
            )
        else:
            typer.echo("  No tables have data to delete.")

        total_rows = sum(row_counts.values())
        if total_rows == 0:
            typer.echo("\n✅ Database is already clean - no rows to delete.")
            return

        # typer.echo(f"\nTotal rows to delete: {total_rows:,}")

        # Confirm unless --force is used
        if not force:
            typer.echo("\n⚠️  WARNING: This operation will permanently delete all data from the specified tables!")
            confirm = typer.confirm("Are you sure you want to proceed?")
            if not confirm:
                typer.echo("❌ Operation cancelled.")
                raise typer.Exit(0)

        # Perform deletion
        typer.echo("\n🧹 Starting database cleanup...")
        deleted_counts = {}

        with engine.connect() as connection:
            for table in tables_to_clean:
                count = row_counts[table]
                if count > 0:
                    try:
                        typer.echo(f"  Deleting from {table}...", nl=False)
                        
                        # Special handling for properties table to preserve corporate IDs
                        if table == "properties":
                            # Delete all properties except corporate_compound_id and corporate_batch_id
                            delete_query = text("""
                                DELETE FROM properties 
                                WHERE name NOT IN ('corporate_compound_id', 'corporate_batch_id')
                            """)
                            result = connection.execute(delete_query)
                            deleted_count = result.rowcount
                        else:
                            # Use TRUNCATE for better performance on large tables
                            # TRUNCATE is faster than DELETE as it doesn't generate individual row events
                            delete_query = text(f"TRUNCATE TABLE {table} RESTART IDENTITY CASCADE")
                            connection.execute(delete_query)
                            # TRUNCATE doesn't return rowcount, so use the pre-operation count
                            deleted_count = count
                        
                        deleted_counts[table] = deleted_count
                        typer.echo(f" ✅ {deleted_count:,} rows deleted")
                    except Exception as e:
                        typer.echo(f" ❌ Error: {e}")
                        deleted_counts[table] = 0
                else:
                    typer.echo(f"  Skipping {table} (already empty)")
                    deleted_counts[table] = 0

            # Commit the transaction after all deletions are complete
            connection.commit()

        # Summary
        typer.echo("\n✅ Database cleanup completed!")

        # Create summary table for deleted rows
        summary_data = [
            {"table": table, "deleted": deleted_counts[table]} for table in tables_to_clean if deleted_counts[table] > 0
        ]

        if summary_data:
            # TODO: look for optimizations
            typer.echo("\nSummary of deleted rows:")
            columns = [("Table Name", "cyan", {"no_wrap": True}), ("Rows Deleted", "green", {"justify": "right"})]

            def extract_summary_row(item):
                return [item["table"], f"{item['deleted']:,}"]

            display_data_table(
                data=summary_data,
                title="Deletion Summary",
                columns=columns,
                row_extractor=extract_summary_row,
                show_total=True,
                total_label="TOTAL ROWS DELETED",
            )
        else:
            typer.echo("\nNo rows were deleted.")

    except Exception as e:
        typer.echo(f"❌ Error during database cleanup: {e}", err=True)
        raise typer.Exit(1)


def display_database_stats_table(table_stats, schema_name):
    """Display database statistics in a rich table format."""
    columns = [("Table Name", "cyan", {"no_wrap": True}), ("Row Count", "green", {"justify": "right"})]

    def extract_stats_row(stat):
        table_name = stat[1]
        row_count = stat[2] or 0
        return [table_name, f"{row_count:,}"]

    display_data_table(
        data=table_stats,
        title=f"Database Statistics - Schema: {schema_name}",
        columns=columns,
        row_extractor=extract_stats_row,
        show_total=True,
    )


@directory_app.command("load")
def load_directory(
    directory_path: str = typer.Argument(..., help="Path to the directory containing files to load"),
    error_handling: str = typer.Option(
        "reject_all", "--error-handling", "-e", help="Error handling strategy: reject_all or reject_row"
    ),
    url: str = DEFAULT_SERVER_URL,
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
            create_assays(
                file_path=str(assays_file),
                url=url
            )
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
                save_errors=save_errors
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
                save_errors=save_errors
            )
            typer.echo("✅ Assay results loaded successfully!")
        except Exception as e:
            typer.echo(f"❌ Error loading assay results: {e}", err=True)
            if error_handling == "reject_all":
                raise typer.Exit(1)
    else:
        typer.echo("⏭️  No assay_results.csv found - skipping assay results loading")
    
    typer.echo("\n✅ Directory loading completed!")


@csv_app.command("analyze")
def analyze_csv_columns(
    input_csv: str = typer.Argument(..., help="Path to the input CSV file"),
    output_csv: str = typer.Argument(..., help="Path to the output CSV file")
):
    """
    Analyze a CSV file's columns and output a CSV with:
    - column name
    - detected type (double, int, datetime, or string)
    - number of null values
    - number of unique values
    """
    def detect_type(values):
        non_nulls = [v for v in values if v.strip() != ""]
        if not non_nulls:
            return "string"
        # Try int
        try:
            for v in non_nulls:
                int(v)
            return "int"
        except Exception:
            pass
        # Try float
        try:
            for v in non_nulls:
                float(v)
            return "double"
        except Exception:
            pass
        # Try datetime, but only if most values parse and values are not all short strings
        parsed_count = 0
        for v in non_nulls:
            try:
                date_parse(v)
                parsed_count += 1
            except Exception:
                pass
        # Only call it datetime if >80% parse and the average value length is at least 8
        if parsed_count / len(non_nulls) > 0.8 and (sum(len(v) for v in non_nulls) / len(non_nulls)) > 8:
            return "datetime"
        return "string"

    with open(input_csv, newline='', encoding='utf-8') as f:
        reader = list(csv.DictReader(f))
        if not reader:
            typer.echo("No data rows found.")
            raise typer.Exit(1)
        columns = reader[0].keys()
        results = []
        for col in columns:
            values = [row[col] for row in reader]
            col_type = detect_type(values)
            null_count = sum(1 for v in values if v.strip() == "")
            unique_count = len(set(v for v in values if v.strip() != ""))
            results.append([col, col_type, null_count, unique_count])
    with open(output_csv, "w", newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["column", "type", "null_count", "unique_count"])
        writer.writerows(results)
    typer.echo(f"Analysis complete. Results written to {output_csv}")


# Search Commands
@search_app.command("compounds-exact")
def search_compounds_exact(
    query_smiles: str = typer.Argument(..., help="SMILES string for exact search"),
    standardization_steps: str = typer.Option(None, "--standardization", help="Standardization steps (comma-separated)"),
    hash_mol: str = typer.Option(None, "--hash-mol", help="Molecular hash for search"),
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Search for compounds using exact SMILES matching.
    """
    payload = {
        "query_smiles": query_smiles
    }
    
    if standardization_steps:
        payload["standardization_steps"] = [step.strip() for step in standardization_steps.split(",")]
    
    if hash_mol:
        payload["hash_mol"] = hash_mol

    response = requests.post(f"{url}/v1/search/compounds/exact", json=payload)

    if response.status_code == 200:
        compounds_data = response.json()

        if output_format == "json":
            print(json.dumps(compounds_data, indent=2))
        else:
            # Display compounds in table format
            display_compounds_table(compounds_data)
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


@search_app.command("compounds-structure")
def search_compounds_structure(
    search_type: str = typer.Argument(..., help="Search type: substructure, tautomer, stereo, similarity, connectivity"),
    query_smiles: str = typer.Argument(..., help="SMILES string for structure search"),
    search_parameters: str = typer.Option(None, "--parameters", help="Additional search parameters as JSON string"),
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Search for compounds using structure-based search.
    """
    payload = {
        "search_type": search_type,
        "query_smiles": query_smiles
    }
    
    if search_parameters:
        try:
            payload["search_parameters"] = json.loads(search_parameters)
        except json.JSONDecodeError:
            typer.echo("Error: Invalid JSON format for search parameters", err=True)
            raise typer.Exit(1)

    response = requests.post(f"{url}/v1/search/compounds/structure", json=payload)

    if response.status_code == 200:
        compounds_data = response.json()

        if output_format == "json":
            print(json.dumps(compounds_data, indent=2))
        else:
            # Display compounds in table format
            display_compounds_table(compounds_data)
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


@search_app.command("complex")
def search_complex(
    file_path: str = typer.Argument(..., help="Path to JSON file containing complex query"),
    url: str = DEFAULT_SERVER_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: json only"),
):
    """
    Perform a complex query across multiple tables.
    """
    # Load and validate the query data
    query_data = load_and_validate_json(file_path)

    response = requests.post(f"{url}/v1/search/complex", json=query_data)

    if response.status_code == 200:
        results_data = response.json()

        if output_format == "json":
            print(json.dumps(results_data, indent=2))
        else:
            print(json.dumps(results_data, indent=2))  # Default to JSON for complex queries
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")


# Admin Commands
@admin_app.command("compound-matching-rule")
def update_compound_matching_rule(
    rule: str = typer.Argument(..., help="Compound matching rule: ALL_LAYERS, STEREO_INSENSITIVE_LAYERS, TAUTOMER_INSENSITIVE_LAYERS"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Update the compound matching rule.
    """
    payload = {"rule": rule}
    response = requests.patch(f"{url}/v1/admin/compound-matching-rule", data=payload)
    print_response(response)


@admin_app.command("institution-id-pattern")
def update_institution_id_pattern(
    scope: str = typer.Argument(..., help="Scope: BATCH or COMPOUND"),
    pattern: str = typer.Argument(..., help="Pattern for generating IDs (e.g., 'DG-{:05d}')"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Update the pattern for generating corporate IDs.
    """
    payload = {"scope": scope, "pattern": pattern}
    response = requests.patch(f"{url}/v1/admin/institution-id-pattern", data=payload)
    print_response(response)


@admin_app.command("molregno-sequence-start")
def set_molregno_sequence_start(
    start_value: int = typer.Argument(..., help="Starting value for molregno sequence"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Set the starting value for the molregno sequence.
    """
    payload = {"start_value": start_value}
    response = requests.patch(f"{url}/v1/admin/molregno-sequence-start", data=payload)
    print_response(response)


@admin_app.command("batchregno-sequence-start")
def set_batchregno_sequence_start(
    start_value: int = typer.Argument(..., help="Starting value for batchregno sequence"),
    url: str = DEFAULT_SERVER_URL,
):
    """
    Set the starting value for the batchregno sequence.
    """
    payload = {"start_value": start_value}
    response = requests.patch(f"{url}/v1/admin/batchregno-sequence-start", data=payload)
    print_response(response)


if __name__ == "__main__":
    app()
