import pytest
from app.utils.enums import ScopeClassReduced, CompoundMatchingRule
from tests.conftest import client


def test_update_institution_id_pattern_valid(client):
    # Test with valid pattern
    response = client.patch(
        "/v1/admin/institution-id-pattern",
        data={
            "scope": "COMPOUND",
            "pattern": "DG-{:09d}"
        }
    )

    assert response.status_code == 200
    assert "Corporate ID pattern for ScopeClassReduced.COMPOUND updated" in response.json()["message"]
    assert "DG-000000001" in response.json()["message"]  # Check that formatting works


def test_update_institution_id_pattern_invalid_format(client):
    # Test invalid pattern (doesn't match expected regex)
    response = client.patch(
        "/v1/admin/institution-id-pattern",
        data={
            "scope": "BATCH",
            "pattern": "INVALID-PATTERN"
        }
    )

    assert response.status_code == 400
    assert "Invalid pattern format" in response.json()["detail"]


def test_update_compound_matching_rule_already_set(client):
    # Assuming the default value in the database is 'ALL_LAYERS'
    response = client.patch(
        "/v1/admin/compound-matching-rule",
        data={"rule": CompoundMatchingRule.ALL_LAYERS.value} 
    )

    assert response.status_code == 200
    assert f"Compound matching rule is already set to {CompoundMatchingRule.ALL_LAYERS.value}" in response.json()["message"]


def test_update_compound_matching_rule_success(client):
    # Update to a different rule succesfully
    response = client.patch(
        "/v1/admin/compound-matching-rule",
        data={"rule": CompoundMatchingRule.STEREO_INSENSITIVE_LAYERS.value} 
    )

    assert response.status_code == 200
    assert f"Compound matching rule updated from {CompoundMatchingRule.ALL_LAYERS.value} to {CompoundMatchingRule.STEREO_INSENSITIVE_LAYERS.value}" in response.json()["message"]
