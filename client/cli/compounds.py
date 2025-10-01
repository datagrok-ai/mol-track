import typer
from client.utils.api_helpers import handle_get_request
from client.utils.display import display_compounds_table
from client.cli.shared import EntityCLI

compound_app = typer.Typer()


class CompoundsCLI(EntityCLI):
    entity_type = "compounds"
    display_fn = staticmethod(display_compounds_table)

    def get_endpoint(self) -> str:
        return "v1/compounds"

    def fetch_entity(self, entity_id: str, url: str):
        params = {"property_value": entity_id, "property_name": "corporate_compound_id"}
        endpoint = f"{url}/{self.get_endpoint()}"
        return handle_get_request(endpoint, params=params)


compound_cli = CompoundsCLI()
compound_cli.register_commands(compound_app)
