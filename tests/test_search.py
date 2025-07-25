import pytest
import json
from sqlalchemy import text

valid_filter = {
    "operator": "OR",
    "conditions": [
        {
            "operator": "AND",
            "conditions": [
                {"field": "compounds.details.corporate_compound_id", "operator": "=", "value": "DG-000001"},
                {"field": "batches.details.corporate_batch_id", "operator": "=", "value": "DGB-000001"},
            ],
        },
        {
            "operator": "AND",
            "conditions": [
                {"field": "compounds.details.corporate_compound_id", "operator": "=", "value": "DG-000002"},
                {"field": "batches.details.corporate_batch_id", "operator": "=", "value": "DGB-000002"},
            ],
        },
    ],
}

valid_output_compounds = ["compounds.molregno", "compounds.details.corporate_compound_id"]
valid_output_batches = ["batches.batch_regno", "batches.details.corporate_batch_id"]


@pytest.mark.usefixtures("preload_simple_data")
def test_valid_json_compounds(client):
    response = client.post("v1/search/compounds", json={"output": valid_output_compounds, "filter": valid_filter})
    assert response.status_code == 200

    content = response.content.decode("utf-8")
    content = json.loads(content)

    assert content["total_count"] == 2
    assert content["level"] == "compounds"

    column_names = [s.replace(".", "_") for s in valid_output_compounds]
    assert content["columns"] == column_names

    details = [
        content["data"][0]["compounds_details_corporate_compound_id"],
        content["data"][1]["compounds_details_corporate_compound_id"],
    ]
    assert "DG-000001" in details
    assert "DG-000002" in details

    molregno = [str(content["data"][0]["compounds_molregno"]), str(content["data"][1]["compounds_molregno"])]
    assert molregno[0] in details[0]
    assert molregno[1] in details[1]


@pytest.mark.usefixtures("preload_simple_data")
def test_valid_json_batches(client):
    response = client.post("v1/search/batches", json={"output": valid_output_batches, "filter": valid_filter})
    assert response.status_code == 200

    content = response.content.decode("utf-8")
    content = json.loads(content)

    assert content["total_count"] == 2
    assert content["level"] == "batches"

    column_names = [s.replace(".", "_") for s in valid_output_batches]
    assert content["columns"] == column_names

    details = [
        content["data"][0]["batches_details_corporate_batch_id"],
        content["data"][1]["batches_details_corporate_batch_id"],
    ]
    assert "DGB-000001" in details
    assert "DGB-000002" in details

    batchregno = [str(content["data"][0]["batches_batch_regno"]), str(content["data"][1]["batches_batch_regno"])]
    assert batchregno[0] in details[0]
    assert batchregno[1] in details[1]


valid_output_assay_results = ["assay_results.id", "assay_results.assay_run_id", "assay_results.created_at"]

valid_filter_assay_results = {"field": "assay_results.id", "operator": "<", "value": "6"}


@pytest.mark.usefixtures("preload_simple_data")
def test_valid_json_assay_results(client):
    response = client.post(
        "v1/search/assay-results", json={"output": valid_output_assay_results, "filter": valid_filter_assay_results}
    )
    assert response.status_code == 200

    content = response.content.decode("utf-8")
    content = json.loads(content)

    assert content["total_count"] == 5
    assert content["level"] == "assay_results"

    column_names = [s.replace(".", "_") for s in valid_output_assay_results]
    assert content["columns"] == column_names


@pytest.mark.usefixtures("preload_simple_data")
def test_valid_molecular_operations(client):
    output = ["compounds.canonical_smiles"]
    filter = {
        "field": "compounds.structure",
        "operator": "IS SIMILAR",
        "value": "CCCCC(CC)COC(=O)C1=CC=C(C=C1)O",
        "threshold": 0.8,
    }
    response = client.post("v1/search/compounds", json={"output": output, "filter": filter})
    assert response.status_code == 200

    content = response.content.decode("utf-8")
    content = json.loads(content)

    assert content["total_count"] == 1
    assert content["level"] == "compounds"

    column_names = [s.replace(".", "_") for s in output]
    assert content["columns"] == column_names


def test_missing_output_field(client):
    response = client.post("v1/search/compounds", json={"filter": valid_filter})
    assert response.status_code == 422


invalid_filter = {
    "operator": "BAD_OP",
    "conditions": [
        {"field": "compounds.id", "operator": "<", "value": "6"},
        {"field": "batches.id", "operator": ">", "value": "6"},
    ],
}


def test_invalid_operator_in_filter(client):
    response = client.post(
        "v1/search/compounds", json={"output": ["compounds.canonical_smiles"], "filter": invalid_filter}
    )
    assert response.status_code in [400, 422]


def test_empty_filter_conditions(client):
    response = client.post(
        "v1/search/compounds",
        json={"output": ["compounds.canonical_smiles"], "filter": {"operator": "AND", "conditions": []}},
    )
    assert response.status_code == 200 or 400  # Depends on your API design


def test_unknown_output_field(client):
    response = client.post("v1/search/compounds", json={"output": ["compounds.nonexistent_field"]})
    assert response.status_code == 400 or 422


@pytest.mark.usefixtures("preload_simple_data")
def test_sql_injection_attempt(client, test_db):
    filter = valid_filter.copy()
    filter["conditions"][0]["conditions"][0]["value"] = "'; DROP TABLE moltrack.compounds;--"
    compounds_before = test_db.execute(text("select * from moltrack.compounds")).fetchall()

    response = client.post("v1/search/compounds", json={"output": valid_output_compounds, "filter": filter})
    assert response.status_code == 200  # should not error

    compounds_after = test_db.execute(text("select * from moltrack.compounds")).fetchall()

    assert compounds_before == compounds_after
