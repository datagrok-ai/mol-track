import json
import requests
import typer

from client.config import settings
from app.utils import enums

admin_app = typer.Typer()


def _update_admin_setting(name: enums.SettingName | str, value: str | int, url: str = settings.API_BASE_URL):
    """Helper for updating admin settings."""
    payload = {"name": name.value if hasattr(name, "value") else name, "value": value}
    response = requests.patch(f"{url}/v1/admin/settings", data=payload)
    response_dict = response.json()

    if response.ok:
        typer.echo(f"✅ {response_dict.get('message', 'Setting updated successfully.')}")
    else:
        typer.secho(
            f"❌ Failed to update setting `{payload['name']}`.\n"
            f"Server response:\n{json.dumps(response_dict, indent=2)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)


# Admin Commands
@admin_app.command("set-compound-matching-rule")
def update_compound_matching_rule(
    rule: enums.CompoundMatchingRule = typer.Argument(..., help="Compound matching rule"),
    url: str = settings.API_BASE_URL,
):
    """Update the compound matching rule."""
    _update_admin_setting(enums.SettingName.COMPOUND_MATCHING_RULE, rule, url)


@admin_app.command("set-institution-id-pattern")
def update_institution_id_pattern(
    entity_type: enums.EntityTypeReduced = typer.Argument(..., help="Entity type"),
    pattern: str = typer.Argument(..., help="Pattern for generating IDs (e.g., 'DG-{:05d}')"),
    url: str = settings.API_BASE_URL,
):
    """Update the pattern for generating corporate IDs."""
    if entity_type not in enums.EntityTypeReduced._value2member_map_:
        typer.secho(
            f"❌ Invalid entity_type `{entity_type}`. Must be one of: "
            f"{', '.join(e.value for e in enums.EntityTypeReduced)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)

    settings_map = {
        enums.EntityTypeReduced.COMPOUND: enums.SettingName.CORPORATE_COMPOUND_ID_PATTERN,
        enums.EntityTypeReduced.BATCH: enums.SettingName.CORPORATE_BATCH_ID_PATTERN,
    }

    _update_admin_setting(settings_map[enums.EntityTypeReduced(entity_type)], pattern, url)


@admin_app.command("set-compound-sequence-start")
def set_molregno_sequence_start(
    start_value: int = typer.Argument(..., help="Starting value for molregno sequence"),
    url: str = settings.API_BASE_URL,
):
    """Set the starting value for the molregno sequence."""
    _update_admin_setting(enums.SettingName.COMPOUND_SEQUENCE_START, start_value, url)


@admin_app.command("set-batch-sequence-start")
def set_batchregno_sequence_start(
    start_value: int = typer.Argument(..., help="Starting value for batchregno sequence"),
    url: str = settings.API_BASE_URL,
):
    """Set the starting value for the batchregno sequence."""
    _update_admin_setting(enums.SettingName.BATCH_SEQUENCE_START, start_value, url)


@admin_app.command("set-compound-friendly-name")
def set_compound_friendly_name(
    friendly_name: str = typer.Argument(..., help="Friendly name for corporate_compound_id"),
    url: str = settings.API_BASE_URL,
):
    """Set the friendly name for corporate_compound_id."""
    _update_admin_setting(enums.SettingName.CORPORATE_COMPOUND_ID_FRIENDLY_NAME, friendly_name, url)


@admin_app.command("set-batch-friendly-name")
def set_batch_friendly_name(
    friendly_name: str = typer.Argument(..., help="Friendly name for corporate_batch_id"),
    url: str = settings.API_BASE_URL,
):
    """Set the friendly name for corporate_batch_id."""
    _update_admin_setting(enums.SettingName.CORPORATE_BATCH_ID_FRIENDLY_NAME, friendly_name, url)
