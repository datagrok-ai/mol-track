import typer

from client.config import settings
from client.utils.api_helpers import run_advanced_search
import click

search_app = typer.Typer()

SEARCH_ENTITIES = {
    "compounds": {
        "endpoint": "/v1/search/compounds",
        "doc": """Advanced search for compounds using /v1/search/compounds endpoint.

The --output parameter can be:
- A comma-separated string: "id,canonical_smiles,common_name"
- A JSON file path containing a list: ["id", "canonical_smiles", "common_name"]
- A JSON file path containing an object: {"output": ["id", "canonical_smiles"]}

Example usage:
    mtcli.py search compounds --output "compounds.id,compounds.canonical_smiles" --filter '{"field": "compounds.details.common_name", "operator": "=", "value": "Aspirin"}'
    mtcli.py search compounds --output output.json --filter filter.json --output-format table
    mtcli.py search compounds --output "compounds.id,compounds.canonical_smiles" --output-format csv
""",
    },
    "batches": {
        "endpoint": "/v1/search/batches",
        "doc": """Advanced search for batches using /v1/search/batches endpoint.

The --output parameter can be:
- A comma-separated string: "id,batch_regno,notes"
- A JSON file path containing a list: ["id", "batch_regno", "notes"]
- A JSON file path containing an object: {"output": ["id", "batch_regno"]}

Example usage:
    mtcli.py search batches --output "batches.id,batches.batch_regno" --filter '{"field": "batches.compound_id", "operator": "=", "value": 1}'
    mtcli.py search batches --output output.json --filter batch_filter.json --output-format table
    mtcli.py search batches --output "batches.id,batches.batch_regno" --output-format csv
""",
    },
    "assay_results": {
        "endpoint": "/v1/search/assay-results",
        "doc": """Advanced search for assay results using /v1/search/assay-results endpoint.

The --output parameter can be:
- A comma-separated string: "id,value_num,assay_id"
- A JSON file path containing a list: ["id", "value_num", "assay_id"]
- A JSON file path containing an object: {"output": ["id", "value_num"]}

Example usage:
    mtcli.py search assay_results --output "assay_results.id,assay_results.details.Reported CLint" --filter '{"field": "assay_results.details.Reported CLint", "operator": ">", "value": 50}'
    mtcli.py search assay_results --output output.json --filter assay_result_filter.json --output-format table
    mtcli.py search assay_results --output "assay_results.id,assay_results.details.Reported CLint" --output-format csv
""",
    },
    "assays": {
        "endpoint": "/v1/search/assays",
        "doc": """Advanced search for assays using /v1/search/assays endpoint.

The --output parameter can be:
- A comma-separated string: "id,value_num,assay_id"
- A JSON file path containing a list: ["id", "value_num", "assay_id"]
- A JSON file path containing an object: {"output": ["id", "value_num"]}

Example usage:
    mtcli.py search assays --output "assays.id,assays.details.assay format" --filter '{"field": "assays.details.assay format", "operator": "=", "value": "in cellulo"}'
    mtcli.py search assays --output output.json --filter assay_filter.json --output-format table
    mtcli.py search assays --output "assays.id,assays.details.assay format" --output-format csv
""",
    },
    "assay_runs": {
        "endpoint": "/v1/search/assay-runs",
        "doc": """Advanced search for assay runs using /v1/search/assay-runs endpoint.

The --output parameter can be:
- A comma-separated string: "id,value_num,assay_id"
- A JSON file path containing a list: ["id", "value_num", "assay_id"]
- A JSON file path containing an object: {"output": ["id", "value_num"]}

Example usage:
    mtcli.py search assay_runs --output "assay_runs.id,assay_runs.details.Cell Species" --filter '{"field": "assay_runs.details.Cell Species", "operator": "=", "value": "Human"}'
    mtcli.py search assay_runs --output output.json --filter assay_run_filter.json --output-format table
    mtcli.py search assay_runs --output "assay_runs.id,assay_runs.details.Cell Species" --output-format csv
""",
    },
}


def create_search_command(app: typer.Typer, search_entities: dict):
    """
    Dynamically create and register search commands for all entities in search_entities.
    """

    for entity_name, info in search_entities.items():
        endpoint = info["endpoint"]
        doc = info["doc"]

        def command_func(
            output: str = typer.Option(
                None,
                "--output",
                "-oc",
                help="Comma-separated list of columns to return or path to JSON file",
            ),
            filter: str = typer.Option(None, "--filter", "-f", help="Filter as JSON string or path to JSON file"),
            aggregations: str = typer.Option(
                None, "--aggregations", "-a", help="Aggregations as JSON string or path to JSON file"
            ),
            url: str = settings.API_BASE_URL,
            output_format: str = typer.Option(
                "json", "--output-format", "-o", help="Output format: table, json, or csv"
            ),
            max_rows: int = typer.Option(
                None, "--max-rows", "-m", help="Maximum number of rows to display in table output"
            ),
            input_file: str = typer.Option(None, "--input-file", "-if", help="Get search input from file"),
            output_file: str = typer.Option(
                None, "--output-file", "-of", help="File to write output to (json, csv, or parquet)"
            ),
            _endpoint=endpoint,
            _entity_name=entity_name,
        ):
            validate_mutually_exclusive(click.get_current_context())
            run_advanced_search(
                _entity_name,
                _endpoint,
                output,
                aggregations,
                filter,
                input_file,
                url,
                output_file,
                output_format,
                max_rows=max_rows,
            )

        command_func.__doc__ = doc
        app.command(entity_name, help=doc)(command_func)


create_search_command(search_app, SEARCH_ENTITIES)


def validate_mutually_exclusive(ctx: typer.Context):
    input_file = ctx.params.get("input_file")
    filter = ctx.params.get("filter")
    output = ctx.params.get("output")
    aggregations = ctx.params.get("aggregations")
    if input_file and (filter or output or aggregations):
        raise typer.BadParameter("You cannot use --input-file and --filter, --output, and --aggregations together")
