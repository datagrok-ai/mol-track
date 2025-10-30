import json
from typing import Any, Callable, Dict, Optional
import requests
from sqlalchemy import text
import typer
from requests.exceptions import RequestException, Timeout


from client.config import settings
from client.utils.data_ingest import parse_arg
from client.utils.display import display_search_csv, display_search_table
from client.utils.file_utils import write_result_to_file, load_input_from_file

try:
    from app.models import SearchRequest
    from app.setup.database import engine, DB_SCHEMA
except ImportError:
    SearchRequest = None
    engine = None
    DB_SCHEMA = None


def print_response(response):
    if response.status_code == 200:
        print(response.json())
    else:
        print(f"Error: {response.status_code}: {response.json()}")


def validate_search_request(level, output, filter_obj, aggregations, output_format, limit):
    """
    Validate the search request using the SearchRequest model if available.
    """
    if SearchRequest is not None:
        try:
            req = SearchRequest(
                level=level,
                output=output,
                filter=filter_obj,
                aggregations=aggregations,
                output_format=output_format,
                limit=limit,
            )
            return req.model_dump()
        except Exception as e:
            typer.secho(f"❌ SearchRequest validation failed: {e}", fg=typer.colors.RED, err=True)
            raise typer.Exit(1)
    else:
        return {
            "level": level,
            "output": output,
            "filter": filter_obj,
            "aggregations": aggregations,
            "output_format": output_format,
            "limit": limit,
        }


def run_advanced_search(
    level, endpoint, output, aggregations, filter, input_file, url, output_file, cli_output_format, max_rows=None
):
    """
    Shared logic for advanced search commands.
    """
    if input_file:
        output, filter, aggregations = load_input_from_file(input_file)
    else:
        output = parse_arg(output, arg_type="output", default_value=[], allow_comma_separated=True)
        filter = filter.replace("'", '"')
        filter = parse_arg(filter, arg_type="json", default_value=None, allow_comma_separated=False)
        aggregations = parse_arg(aggregations, arg_type="json", default_value=None, allow_comma_separated=False)
    search_output_format = cli_output_format if cli_output_format != "table" else "json"
    payload = validate_search_request(level, output, filter, aggregations, search_output_format, max_rows)
    if output_file:
        file_format = output_file.split(".")[-1]
        if file_format != search_output_format:
            typer.secho("❌ Output file extention must match --output-format", fg=typer.colors.RED, err=True)
            raise typer.Exit(1)

    response = requests.post(f"{url}{endpoint}", json=payload)
    if response.status_code == 200:
        write_result_to_file(response, cli_output_format, output_file, parsed=False)
        if cli_output_format == "json":
            resp = response.json()
            print(json.dumps(response.json(), indent=2))
        elif cli_output_format == "csv":
            resp = response.text
            display_search_csv(resp, max_rows=max_rows)
        else:
            resp = response.json()
            display_search_table(resp, max_rows=max_rows)
    else:
        typer.secho(f"❌ Error: {response.status_code}", fg=typer.colors.RED, err=True)
        try:
            error_detail = response.json()
            typer.secho(f"Details: {json.dumps(error_detail, indent=2)}", fg=typer.colors.RED, err=True)
        except Exception:
            typer.secho(f"Response: {response.text}", fg=typer.colors.RED, err=True)


def get_table_row_counts(specific_tables: Optional[list[str]] = None) -> dict[str, int]:
    """
    Get row counts for database tables.

    Args:
        specific_tables: List of specific table names to count. If None, counts all tables.

    Returns:
        Dictionary mapping table names to row counts.
    """
    if engine is None or DB_SCHEMA is None:
        raise ImportError("Database connection not available - engine or DB_SCHEMA is None")

    with engine.connect() as connection:
        # Get all tables in the schema
        tables_query = text("""
            SELECT tablename
            FROM pg_tables
            WHERE schemaname = :schema
            ORDER BY tablename
        """)
        tables_result = connection.execute(tables_query, {"schema": DB_SCHEMA})
        all_tables = [row[0] for row in tables_result] + ["rdk.mols"]

        # Filter tables if specific_tables is provided
        if specific_tables:
            all_tables = [table for table in all_tables if table in specific_tables]

        if not all_tables:
            return {}

        union_queries = " UNION ALL ".join(
            f"SELECT '{table}' AS table_name, COUNT(*) AS row_count FROM {table}" for table in all_tables
        )
        count_query = text(union_queries)
        result = connection.execute(count_query)
        return {row["table_name"]: row["row_count"] for row in result.mappings()}


def handle_request(method: Callable, endpoint: str, **kwargs) -> Dict[str, Any]:
    try:
        response = method(endpoint, timeout=settings.REQUEST_TIMEOUT, **kwargs)
        response.raise_for_status()
    except (RequestException, Timeout) as e:
        try:
            err_detail = response.json().get("detail")
            typer.secho(f"❌ Error: {err_detail}", fg=typer.colors.RED, err=True)
        except Exception:
            typer.secho(f"❌ Request failed: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    try:
        return response.json()
    except (json.JSONDecodeError, ValueError):
        typer.secho(f"❌ Failed to parse JSON. Response: {response.text}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


def handle_get_request(endpoint: str, params: Optional[Dict[str, Any]] = None):
    return handle_request(requests.get, endpoint, params=params)


def handle_delete_request(endpoint: str):
    return handle_request(requests.delete, endpoint)


def handle_put_request(endpoint: str, json_data: Dict[str, Any]):
    return handle_request(requests.put, endpoint, json=json_data)
