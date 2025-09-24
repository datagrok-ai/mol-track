import json
import typer
import requests
from client.utils.api_helpers import handle_get_request
from client.utils.display import display_properties_table
from client.utils.file_utils import load_and_validate_json, write_result_to_file
from app.models import SchemaPayload
from client.config.settings import settings


schema_app = typer.Typer()
schema_list_app = typer.Typer()
schema_app.add_typer(schema_list_app, name="list", help="List schema information")
schema_synonyms_app = typer.Typer()
schema_app.add_typer(schema_synonyms_app, name="synonyms", help="Get schema synonym")


# Schema Commands
@schema_list_app.command("compounds")
def list_compounds_schema(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    List the schema information for compounds.
    """
    print(f"Listing schema from {url}")
    endpoint = f"{url}/v1/schema/compounds"
    get_schema_data(endpoint, output_format, max_rows, output_file)


@schema_list_app.command("batches")
def list_batches_schema(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    List the schema information for batches.
    """
    print(f"Listing schema from {url}")
    endpoint = f"{url}/v1/schema/batches"
    get_schema_data(endpoint, output_format, max_rows, output_file)


@schema_list_app.command("all")
def list_all_schema(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    List all the schema information.
    """
    print(f"Listing schema from {url}")
    endpoint = f"{url}/v1/schema"
    get_schema_data(endpoint, output_format, max_rows, output_file)


@schema_synonyms_app.command("compounds")
def get_compounds_synonyms(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    List synonyms for compounds.
    """
    print(f"Listing schema from {url}")
    endpoint = f"{url}/v1/schema/compounds/synonyms"
    get_schema_data(endpoint, output_format, max_rows, output_file)


@schema_synonyms_app.command("batches")
def get_batches_synonyms(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    output_file: str = typer.Option(None, "--output-file", "-of", help="Path to output file"),
):
    """
    List synonyms for batches.
    """
    print(f"Listing schema from {url}")
    endpoint = f"{url}/v1/schema/batches/synonyms"
    get_schema_data(endpoint, output_format, max_rows, output_file, synonym=True)


@schema_app.command("load")
def add_schema_from_file(
    file_path: str = typer.Argument(..., help="Path to the JSON file containing schema data"),
    url: str = settings.API_BASE_URL,
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
                "entity_type": "COMPOUND",
                "description": "Molecular weight of the compound"
            }
        ],
        "synonym_types": [
            {
                "name": "CAS Number",
                "value_type": "string",
                "property_class": "DECLARED",
                "unit": "",
                "entity_type": "COMPOUND",
                "pattern": "^\\d{1,7}-\\d{2}-\\d$",
                "description": "CAS Registry Number",
                "semantic_type_id": 1
            }
        ]
    }
    """
    # Load and validate schema using utility function
    schema_data = load_and_validate_json(file_path, SchemaPayload)
    typer.echo(schema_data)

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
        typer.echo(f"❌ Error making request: {str(e)}", err=True)
        raise typer.Exit(1)
    except Exception as e:
        typer.echo(f"❌ Error making request: {str(e)}", err=True)
        raise typer.Exit(1)


def get_schema_data(endpoint, output_format, max_rows, output_file, synonym=False):
    schema = handle_get_request(endpoint)

    if output_format == "json":
        typer.echo(f"Schema:\n{json.dumps(schema, indent=2)}")
    else:
        if isinstance(schema, list):
            schema_props = schema
        elif "properties" in schema:
            schema_props = schema["properties"]
        if synonym and "synonym_types" in schema:
            schema_props = schema.get("synonym_types", [])
        typer.echo("Schema:")
        display_properties_table(schema_props, max_rows=max_rows)
    write_result_to_file(schema, output_format, output_file)
