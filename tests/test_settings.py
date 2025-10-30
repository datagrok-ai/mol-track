import pytest
from app.utils import enums

COMPOUND_RULES = enums.CompoundMatchingRule
ENTITY = enums.EntityTypeReduced

test_cases = [
    (
        enums.SettingName.COMPOUND_MATCHING_RULE.value,
        COMPOUND_RULES.ALL_LAYERS.value,
        200,
        f"Compound matching rule is already set to {COMPOUND_RULES.ALL_LAYERS.value}",
    ),
    (
        enums.SettingName.COMPOUND_MATCHING_RULE.value,
        COMPOUND_RULES.STEREO_INSENSITIVE_LAYERS.value,
        200,
        f"Compound matching rule updated from {COMPOUND_RULES.ALL_LAYERS.value} to {COMPOUND_RULES.STEREO_INSENSITIVE_LAYERS.value}",
    ),
    (
        enums.SettingName.CORPORATE_COMPOUND_ID_PATTERN.value,
        "DG-{:09d}",
        200,
        f"Corporate ID pattern for {ENTITY.COMPOUND.value} updated",
    ),
    (enums.SettingName.CORPORATE_BATCH_ID_PATTERN.value, "INVALID-PATTERN", 400, "Invalid pattern format"),
    (
        enums.SettingName.CORPORATE_COMPOUND_ID_FRIENDLY_NAME.value,
        "MyCompoundID",
        200,
        f"Friendly name for {ENTITY.COMPOUND.value} updated",
    ),
    (
        enums.SettingName.CORPORATE_BATCH_ID_FRIENDLY_NAME.value,
        "MyBatchID",
        200,
        f"Friendly name for {ENTITY.BATCH.value} updated",
    ),
    (enums.SettingName.COMPOUND_SEQUENCE_START.value, "1000", 200, "moltrack.molregno_seq set to 1000"),
    (enums.SettingName.BATCH_SEQUENCE_START.value, "500", 200, "moltrack.batch_regno_seq set to 500"),
]


@pytest.mark.parametrize("setting_name, value, expected_status, expected_message_substr", test_cases)
def test_update_settings(client, setting_name, value, expected_status, expected_message_substr):
    response = client.patch("/v1/admin/settings", data={"name": setting_name, "value": value})
    message = response.json().get("message", "") or response.json().get("detail", "")
    assert response.status_code == expected_status
    assert expected_message_substr in message
