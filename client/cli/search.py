import json
import requests
import typer

from client.config import settings
from client.utils.api_helpers import run_advanced_search
from client.utils.display import display_compounds_table
from client.utils.file_utils import load_and_validate_json


search_app = typer.Typer()


# Search Commands
@search_app.command("compounds-exact")
def search_compounds_exact(
    query_smiles: str = typer.Argument(..., help="SMILES string for exact search"),
    standardization_steps: str = typer.Option(
        None, "--standardization", help="Standardization steps (comma-separated)"
    ),
    hash_mol: str = typer.Option(None, "--hash-mol", help="Molecular hash for search"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Search for compounds using exact SMILES matching.
    """
    payload = {"query_smiles": query_smiles}

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
    search_type: str = typer.Argument(
        ..., help="Search type: substructure, tautomer, stereo, similarity, connectivity"
    ),
    query_smiles: str = typer.Argument(..., help="SMILES string for structure search"),
    search_parameters: str = typer.Option(None, "--parameters", help="Additional search parameters as JSON string"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
):
    """
    Search for compounds using structure-based search.
    """
    payload = {"search_type": search_type, "query_smiles": query_smiles}

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
    url: str = settings.API_BASE_URL,
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


@search_app.command("compounds")
def search_compounds(
    output: str = typer.Option(
        ...,
        help="Comma-separated list of columns to return or path to JSON file (e.g. 'id,canonical_smiles' or output.json)",
    ),
    filter: str = typer.Option(None, help="Filter as JSON string or path to JSON file (see docs for format)"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: table, json, or csv"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
):
    """
    Advanced search for compounds using /v1/search/compounds endpoint.

    The --output parameter can be:
    - A comma-separated string: "id,canonical_smiles,common_name"
    - A JSON file path containing a list: ["id", "canonical_smiles", "common_name"]
    - A JSON file path containing an object: {"output": ["id", "canonical_smiles"]}

    Example usage:
      mtcli.py search compounds --output "id,canonical_smiles" --filter '{"field": "compounds.details.common_name", "operator": "=", "value": "Aspirin"}'
      mtcli.py search compounds --output output.json --filter filter.json --output-format table
      mtcli.py search compounds --output "id,canonical_smiles" --output-format csv
    """
    run_advanced_search("compounds", "/v1/search/compounds", output, filter, url, output_format, max_rows=max_rows)


@search_app.command("batches")
def search_batches(
    output: str = typer.Option(
        ...,
        help="Comma-separated list of columns to return or path to JSON file (e.g. 'id,batch_regno' or output.json)",
    ),
    filter: str = typer.Option(None, help="Filter as JSON string or path to JSON file (see docs for format)"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: table, json, or csv"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
):
    """
    Advanced search for batches using /v1/search/batches endpoint.

    The --output parameter can be:
    - A comma-separated string: "id,batch_regno,notes"
    - A JSON file path containing a list: ["id", "batch_regno", "notes"]
    - A JSON file path containing an object: {"output": ["id", "batch_regno"]}

    Example usage:
      mtcli.py search batches --output "id,batch_regno" --filter '{"field": "batches.compound_id", "operator": "=", "value": 1}'
      mtcli.py search batches --output output.json --filter batch_filter.json --output-format table
      mtcli.py search batches --output "id,batch_regno" --output-format csv
    """
    run_advanced_search("batches", "/v1/search/batches", output, filter, url, output_format, max_rows=max_rows)


@search_app.command("assay-results")
def search_assay_results(
    output: str = typer.Option(
        ..., help="Comma-separated list of columns to return or path to JSON file (e.g. 'id,value_num' or output.json)"
    ),
    filter: str = typer.Option(None, help="Filter as JSON string or path to JSON file (see docs for format)"),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: table, json, or csv"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
):
    """
    Advanced search for assay results using /v1/search/assay-results endpoint.

    The --output parameter can be:
    - A comma-separated string: "id,value_num,assay_id"
    - A JSON file path containing a list: ["id", "value_num", "assay_id"]
    - A JSON file path containing an object: {"output": ["id", "value_num"]}

    Example usage:
      mtcli.py search assay-results --output "id,value_num" --filter '{"field": "assay_results.value_num", "operator": ">", "value": 50}'
      mtcli.py search assay-results --output output.json --filter assay_filter.json --output-format table
      mtcli.py search assay-results --output "id,value_num" --output-format csv
    """
    run_advanced_search(
        "assay_results", "/v1/search/assay-results", output, filter, url, output_format, max_rows=max_rows
    )
