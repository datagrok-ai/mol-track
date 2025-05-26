import pandas as pd
import pytest
from tests.conftest import load_csv_file
import enums


@pytest.fixture
def additions_csv():
    file_path = "demo-data/additions_test.csv"
    csv_data = load_csv_file(file_path)
    return file_path, csv_data


def test_create_additions_csv(client, additions_csv):
    file_path, csv_data = additions_csv
    files = {"file": (file_path, csv_data, "text/csv")}
    response = client.post("/v1/additions/", files=files)

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"
    assert "created" in data

    expected_additions = pd.read_csv(file_path)["name"].tolist()
    actual_additions = [addition["name"] for addition in data["created"]["additions"]]
    assert sorted(actual_additions) == sorted(expected_additions)


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
