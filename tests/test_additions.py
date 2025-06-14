import pandas as pd
import pytest
from tests.conftest import load_csv_file
import enums


@pytest.fixture
def additions_csv():
    file_path = "demo-data/additions.csv"
    csv_data = load_csv_file(file_path)
    return file_path, csv_data


def test_create_additions(client, uploaded_additions):
    actual_additions, expected_names = uploaded_additions
    actual_names = [addition["name"] for addition in actual_additions]
    assert sorted(actual_names) == sorted(expected_names)


def test_get_all_additions(client, additions_csv):
    file_path, _ = additions_csv

    response = client.get("/v1/additions/")
    assert response.status_code == 200
    additions = response.json()
    assert isinstance(additions, list)

    expected_names = pd.read_csv(file_path)["name"].tolist()
    actual_names = [item["name"] for item in additions]
    assert sorted(actual_names) == sorted(expected_names)


def test_get_salts(client):
    response = client.get("/v1/additions/salts")
    assert response.status_code == 200
    salts = response.json()
    assert isinstance(salts, list)
    for item in salts:
        item = client.get(f"/v1/additions/{item['id']}").json()
        assert item["role"] == enums.AdditionsRole.SALT.value
    assert len(salts) == 20


def test_get_solvates(client):
    response = client.get("/v1/additions/solvates")
    assert response.status_code == 200
    solvates = response.json()
    assert isinstance(solvates, list)
    for item in solvates:
        item = client.get(f"/v1/additions/{item['id']}").json()
        assert item["role"] == enums.AdditionsRole.SOLVATE.value
    assert len(solvates) == 10
