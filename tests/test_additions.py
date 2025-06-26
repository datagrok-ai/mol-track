from tests.conftest import DATA_DIR
import pandas as pd
from app.utils import enums


def test_get_all_additions(client, preload_additions):
    response = client.get("/v1/additions/")
    assert response.status_code == 200
    additions = response.json()
    assert isinstance(additions, list)

    expected_names = pd.read_csv(DATA_DIR / "additions.csv")["name"].tolist()
    actual_names = [item["name"] for item in additions]
    assert sorted(actual_names) == sorted(expected_names)


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
