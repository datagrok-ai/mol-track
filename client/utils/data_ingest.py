import csv
from datetime import datetime
import json
import os
from pathlib import Path
import tempfile
from typing import Optional
import requests
import typer


def report_csv_information(
    csv_data: list[dict[str, str]],
    entity_type: str,
    mapping_data: Optional[dict[str, any]] = None,
):
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
    url: str,
    endpoint: str,
    entity_type: str,
    mapping_data: Optional[dict[str, any]] = None,
    error_handling: Optional[str] = None,
    csv_data: Optional[list[dict[str, str]]] = None,
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
    output_format = "json"
    try:
        if csv_data is not None:
            # Create temporary CSV file with limited data
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
                file_field = "csv_file" if entity_type == "additions" else "file"
                files = {file_field: (csv_path.name, f, "text/csv")}

                data = {"error_handling": error_handling, "output_format": output_format}

                if mapping_data:
                    data["mapping"] = json.dumps(mapping_data)

                typer.echo(f"🚀 Sending {csv_path.name} to {url}{endpoint}...")

                response = requests.post(f"{url}{endpoint}", files=files, data=data)

                if response.status_code == 200:
                    typer.echo(f"✅ {entity_type.capitalize()} registered successfully!")

                    # Parse the result based on output format
                    data_list = response.json()
                    if entity_type == "additions":
                        data_list = data_list.get("additions", [])
                    success_count = sum(
                        1
                        for item in data_list
                        if item.get("registration_status", "") == "success" or item.get("status", "") == "success"
                    )
                    error_count = len(data_list) - success_count
                    typer.echo(f"📊 Results: {success_count} successful, {error_count} errors")

                    # Show any errors
                    errors = {}
                    for i, item in enumerate(data_list):
                        if item.get("registration_status") != "success":
                            errors[i] = item["registration_error_message"]
                    if errors:
                        typer.echo("❌ Errors found:")
                        for key, value in errors.items():
                            if key == 5:
                                break
                            typer.echo(f"  - Row {key}: {value}")
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


def parse_arg(arg, arg_type="json", default_value=None, allow_comma_separated=False):
    """
    Parse an argument, which can be a JSON string, a path to a JSON file, or a comma-separated string.

    Args:
        arg: The argument to parse
        arg_type: Type of argument ("json" or "output")
        default_value: Value to return if arg is None
        allow_comma_separated: Whether to allow comma-separated string parsing as fallback

    Returns:
        Parsed data (dict for JSON, list for output)
    """
    if arg is None:
        return default_value

    # Try to parse as JSON string
    try:
        return json.loads(arg)
    except Exception:
        # If JSON parsing fails and comma-separated is allowed, try that
        if allow_comma_separated:
            return [col.strip() for col in arg.split(",") if col.strip()]
        else:
            typer.echo(f"Error parsing {arg_type}: Invalid JSON format", err=True)
            raise typer.Exit(1)
