# Import common data from conftest.py (fixtures are automatically available)
from typing import List
import pytest
import pandas as pd
import io
import json
from app.utils import enums


def upload_compounds(client, smiles_list: List[str]) -> None:
    """Uploads a list of SMILES using the /v1/compounds/ endpoint."""
    response = client.get("/v1/compounds/")
    assert response.status_code == 200
    existing = response.json()
    max_id = existing[0]["id"] if existing else 0

    df = pd.DataFrame({"smiles": smiles_list})
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    csv_file = io.BytesIO(csv_bytes)

    mapping = {"smiles": "smiles"}
    files = {"csv_file": ("compounds.csv", csv_file, "text/csv")}
    data = {
        "error_handling": enums.ErrorHandlingOptions.reject_row.value,
        "mapping": json.dumps(mapping),
    }

    response = client.post("/v1/compounds/", files=files, data=data)
    assert response.status_code == 200, f"Failed to upload compounds: {response.text}"

    first_new_id = max_id + 1
    return list(range(first_new_id, first_new_id + len(smiles_list)))


@pytest.mark.skip(reason="This needs to be corrected to work with new search. Old search is no longer available")
def test_create_compound_and_get_hash(client):
    """Test creating a compound and retrieving its hash using the exact search endpoint"""
    smiles = ["C[C@@](F)(Cl)c1cc2ccc[nH]c-2n1"]

    # Step 1: Create the compound
    upload_compounds(client, smiles)

    # Step 2: Use the /v1/search/compounds/exact endpoint to retrieve the hash_mol
    search_payload = {"query_smiles": smiles[0]}
    search_response = client.post("/v1/search/compounds/exact", json=search_payload)
    assert search_response.status_code == 200, f"Unexpected status code: {search_response.status_code}"
    search_results = search_response.json()

    # Verify the hash_mol matches the expected value
    assert len(search_results) == 1, "Expected exactly one result"
    mol_hash = search_results[0]["hash_mol"]
    expected_hash = "9871427f8720e3b4b3219964869a6301377f6908"
    assert mol_hash == expected_hash, f"Expected hash {expected_hash}, got {mol_hash}"


@pytest.fixture()
def predefined_compounds(client):
    """Fixture to insert predefined compounds into the database."""
    smiles_list = [
        "C[C@@](F)(Cl)c1cc2ccc[nH]c-2n1",  # Tautomer1 R
        "C[C@@](F)(Cl)c1cc2cccnc2[nH]1",  # Tautomer1 S
        "C[C@](F)(Cl)c1cc2ccc[nH]c-2n1",  # Tautomer2 R
        "C[C@](F)(Cl)c1cc2cccnc2[nH]1",  # Tautomer2 S
    ]

    compound_ids = upload_compounds(client, smiles_list)
    yield compound_ids


@pytest.mark.skip(reason="This needs to be corrected to work with new search. Old search is no longer available")
def test_search_compound_structure_tautomer(client, predefined_compounds):
    """Test tautomer search using the /v1/search/compounds/structure endpoint"""
    smiles_list = [
        "C[C@@](F)(Cl)c1cc2ccc[nH]c-2n1",  # Tautomer1 R
        "C[C@@](F)(Cl)c1cc2cccnc2[nH]1",  # Tautomer1 S
        "C[C@](F)(Cl)c1cc2ccc[nH]c-2n1",  # Tautomer2 R
        "C[C@](F)(Cl)c1cc2cccnc2[nH]1",  # Tautomer2 S
    ]

    search_payload = {
        "search_type": "tautomer",
        "query_smiles": smiles_list[0],
        "search_parameters": {},
    }
    search_response = client.post("/v1/search/compounds/structure", json=search_payload)
    assert search_response.status_code == 200, f"Unexpected status code: {search_response.status_code}"
    search_results = search_response.json()
    result_ids = [compound["id"] for compound in search_results]

    assert predefined_compounds[0] in result_ids, "Tautomer1 R should match the query"
    assert predefined_compounds[1] in result_ids, "Tautomer1 S should match the query"
    assert predefined_compounds[2] not in result_ids, "Tautomer2 R should NOT match the query"
    assert predefined_compounds[3] not in result_ids, "Tautomer2 S should NOT match the query"


@pytest.mark.skip(reason="This needs to be corrected to work with new search. Old search is no longer available")
def test_search_compound_structure_stereo(client, predefined_compounds):
    """Test stereo search using the /v1/search/compounds/structure endpoint"""
    smiles_list = [
        "C[C@@](F)(Cl)c1cc2ccc[nH]c-2n1",  # Tautomer1 R
        "C[C@@](F)(Cl)c1cc2cccnc2[nH]1",  # Tautomer1 S
        "C[C@](F)(Cl)c1cc2ccc[nH]c-2n1",  # Tautomer2 R
        "C[C@](F)(Cl)c1cc2cccnc2[nH]1",  # Tautomer2 S
    ]

    search_payload = {
        "search_type": "stereo",
        "query_smiles": smiles_list[0],
        "search_parameters": {},
    }
    search_response = client.post("/v1/search/compounds/structure", json=search_payload)
    assert search_response.status_code == 200, f"Unexpected status code: {search_response.status_code}"
    search_results = search_response.json()
    result_ids = [compound["id"] for compound in search_results]

    assert predefined_compounds[0] in result_ids, "Tautomer1 R should match the query"
    assert predefined_compounds[2] in result_ids, "Tautomer2 R should match the query"
    assert predefined_compounds[1] not in result_ids, "Tautomer1 S should NOT match the query"
    assert predefined_compounds[3] not in result_ids, "Tautomer2 S should NOT match the query"


@pytest.mark.skip(reason="This needs to be corrected to work with new search. Old search is no longer available")
def test_search_compound_structure_connectivity(client, predefined_compounds):
    """Test connectivity search using the /v1/search/compounds/structure endpoint"""
    smiles_list = [
        "C[C@@](F)(Cl)c1cc2ccc[nH]c-2n1",  # Tautomer1 R
        "C[C@@](F)(Cl)c1cc2cccnc2[nH]1",  # Tautomer1 S
        "C[C@](F)(Cl)c1cc2ccc[nH]c-2n1",  # Tautomer2 R
        "C[C@](F)(Cl)c1cc2cccnc2[nH]1",  # Tautomer2 S
    ]

    search_payload = {
        "search_type": "connectivity",
        "query_smiles": smiles_list[0],
        "search_parameters": {},
    }
    search_response = client.post("/v1/search/compounds/structure", json=search_payload)
    assert search_response.status_code == 200, f"Unexpected status code: {search_response.status_code}"
    search_results = search_response.json()
    result_ids = [compound["id"] for compound in search_results]

    assert predefined_compounds[0] in result_ids, "Tautomer1 R should match the query"
    assert predefined_compounds[1] in result_ids, "Tautomer1 S should match the query"
    assert predefined_compounds[2] in result_ids, "Tautomer2 R should match the query"
    assert predefined_compounds[3] in result_ids, "Tautomer2 S should match the query"
