import typer

from client.config import settings
from client.utils.api_helpers import run_advanced_search
import click


search_app = typer.Typer()


@search_app.command("compounds")
def search_compounds(
    output: str = typer.Option(
        None,
        "--output",
        "-oc",
        help="Comma-separated list of columns to return or path to JSON file (e.g. 'id,canonical_smiles' or output.json)",
    ),
    filter: str = typer.Option(
        None, "--filter", "-f", help="Filter as JSON string or path to JSON file (see docs for format)"
    ),
    aggregations: str = typer.Option(
        None, "--aggregations", "-a", help="Aggregations as JSON string or path to JSON file (see docs for format)"
    ),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: table, json, or csv"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    input_file: str = typer.Option(None, "--input-file", "-if", help="Get serch input from file"),
    output_file: str = typer.Option(
        None, "--output-file", "-of", help="File to write output to (json, csv, or parquet)"
    ),
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
    validate_mutually_exclusive(click.get_current_context())
    run_advanced_search(
        "compounds",
        "/v1/search/compounds",
        output,
        aggregations,
        filter,
        input_file,
        url,
        output_file,
        output_format,
        max_rows=max_rows,
    )


@search_app.command("batches")
def search_batches(
    output: str = typer.Option(
        None,
        "--output",
        "-oc",
        help="Comma-separated list of columns to return or path to JSON file (e.g. 'id,canonical_smiles' or output.json)",
    ),
    filter: str = typer.Option(
        None, "--filter", "-f", help="Filter as JSON string or path to JSON file (see docs for format)"
    ),
    aggregations: str = typer.Option(
        None, "--aggregations", "-a", help="Aggregations as JSON string or path to JSON file (see docs for format)"
    ),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: table, json, or csv"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    input_file: str = typer.Option(None, "--input-file", "-if", help="Get serch input from file"),
    output_file: str = typer.Option(
        None, "--output-file", "-of", help="File to write output to (json, csv, or parquet)"
    ),
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

    validate_mutually_exclusive(click.get_current_context())
    run_advanced_search(
        "batches",
        "/v1/search/batches",
        output,
        aggregations,
        filter,
        input_file,
        url,
        output_file,
        output_format,
        max_rows=max_rows,
    )


@search_app.command("assay-results")
def search_assay_results(
    output: str = typer.Option(
        None,
        "--output",
        "-oc",
        help="Comma-separated list of columns to return or path to JSON file (e.g. 'id,canonical_smiles' or output.json)",
    ),
    filter: str = typer.Option(
        None, "--filter", "-f", help="Filter as JSON string or path to JSON file (see docs for format)"
    ),
    aggregations: str = typer.Option(
        None, "--aggregations", "-a", help="Aggregations as JSON string or path to JSON file (see docs for format)"
    ),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: table, json, or csv"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    input_file: str = typer.Option(None, "--input-file", "-if", help="Get serch input from file"),
    output_file: str = typer.Option(
        None, "--output-file", "-of", help="File to write output to (json, csv, or parquet)"
    ),
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
    validate_mutually_exclusive(click.get_current_context())
    run_advanced_search(
        "assay_results",
        "/v1/search/assay-results",
        output,
        aggregations,
        filter,
        input_file,
        url,
        output_file,
        output_format,
        max_rows=max_rows,
    )


@search_app.command("assays")
def search_assays(
    output: str = typer.Option(
        None,
        "--output",
        "-oc",
        help="Comma-separated list of columns to return or path to JSON file (e.g. 'id,canonical_smiles' or output.json)",
    ),
    filter: str = typer.Option(
        None, "--filter", "-f", help="Filter as JSON string or path to JSON file (see docs for format)"
    ),
    aggregations: str = typer.Option(
        None, "--aggregations", "-a", help="Aggregations as JSON string or path to JSON file (see docs for format)"
    ),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: table, json, or csv"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    input_file: str = typer.Option(None, "--input-file", "-if", help="Get serch input from file"),
    output_file: str = typer.Option(
        None, "--output-file", "-of", help="File to write output to (json, csv, or parquet)"
    ),
):
    """
    Advanced search for assays using /v1/search/assays endpoint.

    The --output parameter can be:
    - A comma-separated string: "id,value_num,assay_id"
    - A JSON file path containing a list: ["id", "value_num", "assay_id"]
    - A JSON file path containing an object: {"output": ["id", "value_num"]}

    Example usage:
      mtcli.py search assays --output "id,value_num" --filter '{"field": "assay_results.value_num", "operator": ">", "value": 50}'
      mtcli.py search assays --output output.json --filter assay_filter.json --output-format table
      mtcli.py search assays --output "id,value_num" --output-format csv
    """
    validate_mutually_exclusive(click.get_current_context())
    run_advanced_search(
        "assays",
        "/v1/search/assays",
        output,
        aggregations,
        filter,
        input_file,
        url,
        output_file,
        output_format,
        max_rows=max_rows,
    )


@search_app.command("assay-runs")
def search_assay_runs(
    output: str = typer.Option(
        None,
        "--output",
        "-oc",
        help="Comma-separated list of columns to return or path to JSON file (e.g. 'id,canonical_smiles' or output.json)",
    ),
    filter: str = typer.Option(
        None, "--filter", "-f", help="Filter as JSON string or path to JSON file (see docs for format)"
    ),
    aggregations: str = typer.Option(
        None, "--aggregations", "-a", help="Aggregations as JSON string or path to JSON file (see docs for format)"
    ),
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("json", "--output-format", "-o", help="Output format: table, json, or csv"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
    input_file: str = typer.Option(None, "--input-file", "-if", help="Get serch input from file"),
    output_file: str = typer.Option(
        None, "--output-file", "-of", help="File to write output to (json, csv, or parquet)"
    ),
):
    """
    Advanced search for assays using /v1/search/assays endpoint.

    The --output parameter can be:
    - A comma-separated string: "id,value_num,assay_id"
    - A JSON file path containing a list: ["id", "value_num", "assay_id"]
    - A JSON file path containing an object: {"output": ["id", "value_num"]}

    Example usage:
      mtcli.py search assays --output "id,value_num" --filter '{"field": "assay_results.value_num", "operator": ">", "value": 50}'
      mtcli.py search assays --output output.json --filter assay_filter.json --output-format table
      mtcli.py search assays --output "id,value_num" --output-format csv
    """
    validate_mutually_exclusive(click.get_current_context())
    run_advanced_search(
        "assay_runs",
        "/v1/search/assay-runs",
        output,
        aggregations,
        filter,
        input_file,
        url,
        output_file,
        output_format,
        max_rows=max_rows,
    )


def validate_mutually_exclusive(ctx: typer.Context):
    input_file = ctx.params.get("input_file")
    filter = ctx.params.get("filter")
    output = ctx.params.get("output")
    aggregations = ctx.params.get("aggregations")
    if input_file and (filter or output or aggregations):
        raise typer.BadParameter("You cannot use --input-file and --filter, --output, and --aggregations together")
