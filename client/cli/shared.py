from abc import ABC, abstractmethod
import typer
import json
from client.utils.data_ingest import send_csv_upload_request, report_csv_information
from client.utils.api_helpers import handle_get_request, handle_delete_request
from client.utils.file_utils import (
    validate_and_load_csv_data,
    write_result_to_file,
)


class EntityCLI(ABC):
    """Base interface for entity CLI commands."""

    entity_type: str
    display_fn: callable

    @abstractmethod
    def get_endpoint(self) -> str:
        """Return the API endpoint for this entity."""
        pass

    def list(
        self,
        skip: int = 0,
        limit: int = 10,
        url: str = None,
        output_format: str = "table",
        output_file: str | None = None,
    ):
        endpoint = f"{url}/{self.get_endpoint()}/?skip={skip}&limit={limit}"
        data = handle_get_request(endpoint)
        if output_format == "json":
            typer.echo(json.dumps(data, indent=2))
        else:
            self.display_fn(data)
        write_result_to_file(data, output_format, output_file)

    def get(self, entity_id: str, url: str = None, output_format: str = "table", output_file: str | None = None):
        param_name = f"corporate_{self.entity_type[:-1]}_id"
        params = {"property_value": entity_id, "property_name": param_name}
        endpoint = f"{url}/{self.get_endpoint()}"
        data = handle_get_request(endpoint, params=params)
        if output_format == "json":
            typer.echo(json.dumps(data, indent=2))
        else:
            self.display_fn([data])
            if "properties" in data:
                from client.utils.display import display_properties_table

                display_properties_table(data["properties"], display_value=True)
        write_result_to_file(data, output_format, output_file)

    def delete(self, entity_id: str, url: str = None):
        endpoint = f"{url}/{self.get_endpoint()}/{entity_id}"
        handle_delete_request(endpoint)
        typer.echo(f"✅ {self.entity_type[:-1].capitalize()} with ID {entity_id} deleted successfully.")

    def load_from_csv(
        self,
        csv_file: str,
        mapping_file: str | None = None,
        rows: int | None = None,
        url: str = None,
        error_handling: str = "reject_all",
        dry_run: bool = False,
        save_errors: bool = False,
    ):
        csv_path, csv_data = validate_and_load_csv_data(csv_file, self.entity_type, rows)
        mapping_data = None
        if mapping_file:
            from client.utils.file_utils import load_and_validate_mapping

            mapping_data = load_and_validate_mapping(mapping_file)
        report_csv_information(csv_data, mapping_data, self.entity_type)

        if dry_run:
            typer.echo("✅ Dry run completed successfully!")
            return

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
