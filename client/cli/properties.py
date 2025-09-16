import json
import requests
import typer

from client.config import settings
from client.utils.display import display_properties_table


properties_app = typer.Typer()


# List Commands
@properties_app.command("list")
def list_properties(
    url: str = settings.API_BASE_URL,
    output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
    max_rows: int = typer.Option(None, "--max-rows", "-m", help="Maximum number of rows to display in table output"),
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
            display_properties_table(properties_data, max_rows=max_rows)
    else:
        print(f"Error: {response.status_code}")
        try:
            error_detail = response.json()
            print(f"Details: {json.dumps(error_detail, indent=2)}")
        except (json.JSONDecodeError, ValueError):
            print(f"Response: {response.text}")
