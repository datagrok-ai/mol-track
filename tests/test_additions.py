import pytest
from tests.conftest import DATA_DIR
import pandas as pd
from app.utils import enums


expected_additions = pd.read_csv(DATA_DIR / "additions.csv").to_dict(orient="records")


def _get_all_additions(client):
    response = client.get("/v1/additions/")
    assert response.status_code == 200
    additions = response.json()
    assert isinstance(additions, list)
    return additions


@pytest.mark.parametrize("expected", expected_additions)
def test_get_all_additions(client, preload_additions, expected):
    additions = _get_all_additions(client)
    api_names = [a["name"] for a in additions]
    assert expected["name"] in api_names


def test_get_salts(client, preload_additions):
    response = client.get("/v1/additions/salts")
    assert response.status_code == 200
    salts = response.json()
    assert isinstance(salts, list)
    for item in salts:
        item = client.get(f"/v1/additions/{item['id']}").json()
        assert item["role"] == enums.AdditionsRole.SALT.value
    assert len(salts) == 20


def test_get_solvates(client, preload_additions):
    response = client.get("/v1/additions/solvates")
    assert response.status_code == 200
    solvates = response.json()
    assert isinstance(solvates, list)
    for item in solvates:
        item = client.get(f"/v1/additions/{item['id']}").json()
        assert item["role"] == enums.AdditionsRole.SOLVATE.value
    assert len(solvates) == 10


@pytest.mark.parametrize("expected", expected_additions)
def test_get_addition_by_id(client, preload_additions, expected):
    additions = _get_all_additions(client)

    api_addition = next(a for a in additions if a["name"] == expected["name"])

    response = client.get(f"/v1/additions/{api_addition['id']}")
    assert response.status_code == 200
    addition = response.json()

    assert addition["id"] == api_addition["id"]
    assert addition["name"] == expected["name"]


@pytest.fixture
def first_addition_id(client, preload_additions):
    additions = _get_all_additions(client)
    return additions[0]["id"]


@pytest.mark.parametrize(
    "update_payload",
    [
        ({"name": "updated addition name"}),
        ({"description": "updated addition description"}),
    ],
)
def test_put_addition_by_id(client, first_addition_id, update_payload):
    response = client.put(f"/v1/additions/{first_addition_id}", json=update_payload)
    assert response.status_code == 200
    data = response.json()

    for key, value in update_payload.items():
        assert data[key] == value, f"Expected {key} to be {value}, but got {data[key]}"


def test_delete_addition_by_id(client, first_addition_id):
    response_delete = client.delete(f"/v1/additions/{first_addition_id}")
    assert response_delete.status_code == 200

    response_get = client.get(f"/v1/additions/{first_addition_id}")
    assert response_get.status_code == 404

    response_get_data = response_get.json()
    assert "detail" in response_get_data
    assert response_get_data["detail"] == "Addition not found"
