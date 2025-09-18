import json
import requests
import typer

from client.config import settings


admin_app = typer.Typer()


# Admin Commands
@admin_app.command("set-compound-matching-rule")
def update_compound_matching_rule(
    rule: str = typer.Argument(
        ..., help="Compound matching rule: ALL_LAYERS, STEREO_INSENSITIVE_LAYERS, TAUTOMER_INSENSITIVE_LAYERS"
    ),
    url: str = settings.API_BASE_URL,
):
    """
    Update the compound matching rule.
    """
    payload = {"name": "COMPOUND_MATCHING_RULE", "value": rule}
    response = requests.patch(f"{url}/v1/admin/settings", data=payload)
    response_dict = response.json()
    if response.status_code == 200:
        typer.echo(f"✅ {response_dict['message']}")
    else:
        typer.secho(
            f"❌ Error updating the rule. Error message:\n {json.dumps(response_dict, indent=2)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)


@admin_app.command("set-institution-id-pattern")
def update_institution_id_pattern(
    entity_type: str = typer.Argument(..., help="entity_type: BATCH or COMPOUND"),
    pattern: str = typer.Argument(..., help="Pattern for generating IDs (e.g., 'DG-{:05d}')"),
    url: str = settings.API_BASE_URL,
):
    """
    Update the pattern for generating corporate IDs.
    """
    if entity_type not in ["COMPOUND", "BATCH"]:
        typer.secho(
            "❌ Error updating the pattern. entity_type must be BATCH or COMPOUND", fg=typer.colors.RED, err=True
        )
        raise typer.Exit(code=1)
    settings_map = {"COMPOUND": "CORPORATE_COMPOUND_ID_PATTERN", "BATCH": "CORPORATE_BATCH_ID_PATTERN"}
    payload = {"name": settings_map[entity_type], "value": pattern}
    response = requests.patch(f"{url}/v1/admin/settings", data=payload)
    response_dict = response.json()
    if response.status_code == 200:
        typer.echo(f"✅ {response_dict['message']}")
    else:
        typer.secho(
            f"❌ Error updating the rule. Error message:\n {json.dumps(response_dict, indent=2)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)


@admin_app.command("set-compound-sequence-start")
def set_molregno_sequence_start(
    start_value: int = typer.Argument(..., help="Starting value for molregno sequence"),
    url: str = settings.API_BASE_URL,
):
    """
    Set the starting value for the molregno sequence.
    """
    payload = {"name": "COMPOUND_SEQUENCE_START", "value": start_value}
    response = requests.patch(f"{url}/v1/admin/settings", data=payload)
    response_dict = response.json()
    if response.status_code == 200:
        typer.echo(f"✅ {response_dict['message']}")
    else:
        typer.secho(
            f"❌ Error updating the rule. Error message:\n {json.dumps(response_dict, indent=2)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)


@admin_app.command("set-batch-sequence-start")
def set_batchregno_sequence_start(
    start_value: int = typer.Argument(..., help="Starting value for batchregno sequence"),
    url: str = settings.API_BASE_URL,
):
    """
    Set the starting value for the batchregno sequence.
    """
    payload = {"name": "BATCH_SEQUENCE_START", "value": start_value}
    response = requests.patch(f"{url}/v1/admin/settings", data=payload)
    response_dict = response.json()
    if response.status_code == 200:
        typer.echo(f"✅ {response_dict['message']}")
    else:
        typer.secho(
            f"❌ Error updating the rule. Error message:\n {json.dumps(response_dict, indent=2)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)
