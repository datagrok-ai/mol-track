import json
from typing import Any, Dict
import requests
from sqlalchemy import text
import typer
from requests.exceptions import RequestException, Timeout


from client.config import settings
from client.utils.data_ingest import parse_arg
from client.utils.display import display_search_csv, display_search_table

try:
    from app.models import SearchRequest, SearchResponse
    from app.setup.database import engine, DB_SCHEMA
except ImportError:
    print("app.models import failed")
    SearchRequest = None
    SearchResponse = None
    engine = None
    DB_SCHEMA = None


def print_response(response):
    if response.status_code == 200:
        print(response.json())
    else:
        print(f"Error: {response.status_code}: {response.json()}")


def validate_search_request(level, output, filter_obj):
    """
    Validate the search request using the SearchRequest model if available.
    """
    if SearchRequest is not None:
        try:
            req = SearchRequest(level=level, output=output, filter=filter_obj)
            return req.model_dump()
        except Exception as e:
            typer.echo(f"❌ SearchRequest validation failed: {e}", err=True)
            raise typer.Exit(1)
    else:
        return {"level": level, "output": output, "filter": filter_obj}


def validate_search_response(data):
    """
    Validate the search response using the SearchResponse model if available.
    """
    if SearchResponse is not None:
        try:
            resp = SearchResponse(**data)
            return resp
        except Exception as e:
            typer.echo(f"❌ SearchResponse validation failed: {e}", err=True)
            raise typer.Exit(1)
    else:
        return data


def run_advanced_search(level, endpoint, output, filter, url, output_format, max_rows=None):
    """
    Shared logic for advanced search commands.
    """
    output_list = parse_arg(output, arg_type="output", default_value=[], allow_comma_separated=True)
    filter_obj = parse_arg(filter, arg_type="json", default_value=None, allow_comma_separated=False)
    payload = validate_search_request(level, output_list, filter_obj)
    response = requests.post(f"{url}{endpoint}", json=payload)
    if response.status_code == 200:
        resp = validate_search_response(response.json())
        if output_format == "json":
            print(json.dumps(response.json(), indent=2))
        elif output_format == "csv":
            display_search_csv(resp, output_list, max_rows=max_rows)
        else:
            display_search_table(resp, output_list, max_rows=max_rows)
    else:
        typer.echo(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            typer.echo(f"Details: {json.dumps(error_detail, indent=2)}")
        except Exception:
            typer.echo(f"Response: {response.text}")


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
        all_tables = [row[0] for row in tables_result] + ["rdk.mols"]

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


def handle_get_request(endpoint: str, params: Dict[str, Any] = None):
    try:
        response = requests.get(endpoint, params=params, timeout=settings.REQUEST_TIMEOUT)
        response.raise_for_status()
    except (RequestException, Timeout) as e:
        typer.secho(f"❌ Error: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    try:
        return response.json()
    except (json.JSONDecodeError, ValueError):
        typer.secho(f"❌ Failed to parse JSON. Response: {response.text}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


def handle_delete_request(endpoint: str):
    try:
        response = requests.delete(endpoint)
        resp_json = response.json()
        response.raise_for_status()
        return resp_json
    except (RequestException, Timeout):
        typer.secho(f"❌ Error: {resp_json['detail']}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except (json.JSONDecodeError, ValueError):
        typer.secho(f"❌ Failed to parse JSON. Response: {response.text}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
