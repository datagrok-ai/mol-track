import requests
import typer

from client.config import settings
from client.utils.api_helpers import print_response
from client.utils.display import display_assays_table
from client.utils.file_utils import load_and_validate_json
from client.cli.shared import EntityCLI


assays_app = typer.Typer()
assays_runs_app = typer.Typer()
assays_results_app = typer.Typer()
assays_app.add_typer(assays_runs_app, name="runs", help="Assay runs management")
assays_app.add_typer(assays_results_app, name="results", help="Assay results management")


# Assays Commands
class AssaysCLI(EntityCLI):
    entity_type = "assays"

    @staticmethod
    def display_fn(data):
        display_assays_table(data)

    def get_endpoint(self) -> str:
        return "v1/assays"

    def fetch_entity(self, entity_id: str, url: str):
        from client.utils.api_helpers import handle_get_request

        return handle_get_request(f"{url}/{self.get_endpoint()}/{entity_id}")

    def load_assays(self, file_path: str, url: str):
        """
        Load assays from a JSON file using the /v1/assays endpoint.
        """
        assay_data = load_and_validate_json(file_path)
        response = requests.post(f"{url}/{self.get_endpoint()}", json=assay_data)
        print_response(response)

    def register_commands(self, app: typer.Typer):
        super().register_commands(app)

        @app.command("load")
        def load_assay_entities(
            file_path: str = typer.Argument(..., help="Path to the JSON file containing assay data"),
            url: str = typer.Option(settings.API_BASE_URL, help="API base URL"),
        ):
            self.load_assays(file_path, url)


assays_cli = AssaysCLI()
assays_cli.register_commands(assays_app)


# Assay Runs Commands
class AssayRunsCLI(EntityCLI):
    entity_type = "assay_runs"

    @staticmethod
    def display_fn(data):
        display_assays_table(data, assay_entity="run")

    def get_endpoint(self) -> str:
        return "v1/assay_runs"

    def fetch_entity(self, entity_id: str, url: str):
        from client.utils.api_helpers import handle_get_request

        return handle_get_request(f"{url}/{self.get_endpoint()}/{entity_id}")


assay_runs_cli = AssayRunsCLI()
assay_runs_cli.register_commands(assays_runs_app)


# Assay Results Commands
class AssayResultsCLI(EntityCLI):
    entity_type = "assay_results"

    @staticmethod
    def display_fn(data):
        display_assays_table(data, assay_entity="result")

    def get_endpoint(self) -> str:
        return "v1/assay_results"

    def fetch_entity(self, entity_id: str, url: str):
        from client.utils.api_helpers import handle_get_request

        return handle_get_request(f"{url}/{self.get_endpoint()}/{entity_id}")


assay_results_cli = AssayResultsCLI()
assay_results_cli.register_commands(assays_results_app)
