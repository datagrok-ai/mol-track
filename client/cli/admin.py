import requests
import typer

from client.config import settings
from client.utils.api_helpers import print_response


admin_app = typer.Typer()


# Admin Commands
@admin_app.command("compound-matching-rule")
def update_compound_matching_rule(
    rule: str = typer.Argument(
        ..., help="Compound matching rule: ALL_LAYERS, STEREO_INSENSITIVE_LAYERS, TAUTOMER_INSENSITIVE_LAYERS"
    ),
    url: str = settings.API_BASE_URL,
):
    """
    Update the compound matching rule.
    """
    payload = {"rule": rule}
    response = requests.patch(f"{url}/v1/admin/compound-matching-rule", data=payload)
    print_response(response)


@admin_app.command("institution-id-pattern")
def update_institution_id_pattern(
    entity_type: str = typer.Argument(..., help="entity_type: BATCH or COMPOUND"),
    pattern: str = typer.Argument(..., help="Pattern for generating IDs (e.g., 'DG-{:05d}')"),
    url: str = settings.API_BASE_URL,
):
    """
    Update the pattern for generating corporate IDs.
    """
    payload = {"entity_type": entity_type, "pattern": pattern}
    response = requests.patch(f"{url}/v1/admin/institution-id-pattern", data=payload)
    print_response(response)


@admin_app.command("molregno-sequence-start")
def set_molregno_sequence_start(
    start_value: int = typer.Argument(..., help="Starting value for molregno sequence"),
    url: str = settings.API_BASE_URL,
):
    """
    Set the starting value for the molregno sequence.
    """
    payload = {"start_value": start_value}
    response = requests.patch(f"{url}/v1/admin/molregno-sequence-start", data=payload)
    print_response(response)


@admin_app.command("batchregno-sequence-start")
def set_batchregno_sequence_start(
    start_value: int = typer.Argument(..., help="Starting value for batchregno sequence"),
    url: str = settings.API_BASE_URL,
):
    """
    Set the starting value for the batchregno sequence.
    """
    payload = {"start_value": start_value}
    response = requests.patch(f"{url}/v1/admin/batchregno-sequence-start", data=payload)
    print_response(response)
