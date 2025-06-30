import pytest
from app.utils.enums import ScopeClassReduced
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
    print(response.json())
    assert response.status_code == 400
    assert "Invalid pattern format" in response.json()["detail"]