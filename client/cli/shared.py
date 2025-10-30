from abc import ABC, abstractmethod
import typer
import json
from pathlib import Path
from client.config.settings import settings
from client.utils.api_helpers import handle_get_request, handle_delete_request
from client.utils.file_utils import (
    validate_and_load_csv_data,
    write_result_to_file,
    load_and_validate_mapping,
)
from client.utils.data_ingest import send_csv_upload_request, report_csv_information
from app.utils import enums


class EntityCLI(ABC):
    """Base interface for entity CLI commands with built-in list/get/delete/load."""

    entity_type: str
    display_fn: callable

    @abstractmethod
    def get_endpoint(self) -> str:
        """Return the API endpoint for this entity."""
        return f"v1/{self.entity_type}"

    @abstractmethod
    def fetch_entity(self, entity_id: str, url: str):
        """
        Fetch an entity from the server. Must be implemented by subclasses.
        """
        pass

    def register_commands(self, app: typer.Typer):
        @app.command("list")
        def _list(
            skip: int = typer.Option(0, "--skip", "-s", help="Number of records to skip"),
            limit: int = typer.Option(10, "--limit", "-l", help="Maximum number of records to return"),
            url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
            output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
            output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
        ):
            endpoint = f"{url}/{self.get_endpoint()}/?skip={skip}&limit={limit}"
            data = handle_get_request(endpoint)

            if output_format == "json":
                typer.echo(json.dumps(data, indent=2))
            elif self.display_fn:
                self.display_fn(data)
            write_result_to_file(data, output_format, output_file)

        @app.command("get")
        def _get(
            entity_id: str = typer.Argument(
                ..., help=f"{self.entity_type[:-1].capitalize()} ID (friendly name) to retrieve"
            ),
            url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
            output_format: str = typer.Option("table", "--output-format", "-o", help="Output format: table or json"),
            output_file: str | None = typer.Option(None, "--output-file", "-of", help="Path to output file"),
        ):
            data = self.fetch_entity(entity_id, url)

            if output_format == "json":
                typer.echo(json.dumps(data, indent=2))
            elif self.display_fn:
                self.display_fn([data])
                if "properties" in data:
                    from client.utils.display import display_properties_table

                    display_properties_table(data["properties"], display_value=True)

            write_result_to_file(data, output_format, output_file)

        @app.command("delete")
        def _delete(
            entity_id: str = typer.Argument(
                ..., help=f"{self.entity_type[:-1].capitalize()} ID (friendly name) to delete"
            ),
            url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
        ):
            endpoint = f"{url}/{self.get_endpoint()}/{entity_id}"
            handle_delete_request(endpoint)
            typer.echo(f"✅ {self.entity_type[:-1].capitalize()} with ID {entity_id} deleted successfully.")

        @app.command("load")
        def _load(
            csv_file: Path = typer.Argument(..., help="Path to the CSV file containing entity data"),
            mapping_file: Path | None = typer.Option(
                None, "--mapping", "-m", help="Path to the JSON mapping file (optional)"
            ),
            rows: int | None = typer.Option(
                None, "--rows", "-r", help="Number of data rows to process (excludes header row)"
            ),
            url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
            error_handling: enums.ErrorHandlingOptions = typer.Option(
                enums.ErrorHandlingOptions.reject_all, "--error-handling", "-e", help="Error handling strategy"
            ),
            dry_run: bool = typer.Option(False, "--dry-run", help="Validate data without sending to server"),
            save_errors: bool = typer.Option(False, "--save-errors", help="Save error records to a JSON file"),
        ):
            self.load_entity(
                csv_file=csv_file,
                mapping_file=mapping_file,
                rows=rows,
                url=url,
                error_handling=error_handling,
                dry_run=dry_run,
                save_errors=save_errors,
            )

    def load_entity(
        self,
        csv_file: Path,
        url: str,
        mapping_file: Path | None = None,
        rows: int | None = None,
        error_handling: enums.ErrorHandlingOptions = enums.ErrorHandlingOptions.reject_all,
        dry_run: bool = False,
        save_errors: bool = False,
    ):
        error_handling = (
            error_handling.value if isinstance(error_handling, enums.ErrorHandlingOptions) else error_handling
        )
        csv_path, csv_data = validate_and_load_csv_data(csv_file, self.entity_type, rows)
        mapping_data = load_and_validate_mapping(mapping_file) if mapping_file else None
        report_csv_information(csv_data, self.entity_type, mapping_data)

        if dry_run:
            typer.echo("✅ Dry run completed successfully!")
            return

        print("error handling", error_handling)
        print(type(error_handling))

        send_csv_upload_request(
            csv_path=csv_path,
            mapping_data=mapping_data,
            url=url,
            endpoint=f"/{self.get_endpoint()}/",
            error_handling=error_handling,
            entity_type=self.entity_type,
            csv_data=csv_data if rows else None,
            save_errors=save_errors,
        )
