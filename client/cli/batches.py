import typer
from client.cli.shared import EntityCLI
from client.utils.api_helpers import handle_get_request
from client.utils.display import display_batches_table


batch_app = typer.Typer()


class BatchCLI(EntityCLI):
    entity_type = "batches"
    display_fn = staticmethod(display_batches_table)

    def get_endpoint(self) -> str:
        return "v1/batches"

    def fetch_entity(self, entity_id: str, url: str):
        params = {"property_value": entity_id, "property_name": "corporate_batch_id"}
        endpoint = f"{url}/{self.get_endpoint()}"
        return handle_get_request(endpoint, params=params)


batch_cli = BatchCLI()
batch_cli.register_commands(batch_app)
