# Import common data from conftest.py (fixtures are automatically available)
import pytest

from tests.conftest import aspirin_smiles, aspirin_smiles_noncanonical, caffeine_smiles

# Test data
test_compound_data = {
    "smiles": aspirin_smiles,  # Aspirin SMILES
    "is_archived": False,
}


def test_create_compound(client):
    """Test creating a new compound via the API endpoint"""
    response = client.post("/compounds/", json=test_compound_data)

    # Check response status code and content
    assert response.status_code == 200
    data = response.json()

    # Verify the returned data contains the expected values
    assert "canonical_smiles" in data
    assert "inchi" in data
    assert "inchikey" in data
    assert data["is_archived"] == test_compound_data["is_archived"]

    # Verify the compound has an ID (was saved to the database)
    assert "id" in data
    assert isinstance(data["id"], int)


def test_create_compound_duplicate_smiles(client):
    """Test that creating a compound with a duplicate SMILES string fails appropriately"""
    # First create a compound
    client.post("/compounds/", json=test_compound_data)

    # Try to create another compound with the same SMILES
    response = client.post("/compounds/", json=test_compound_data)

    # Check if the API handles duplicates appropriately
    assert response.status_code == 400
    assert "already exists" in response.json()["detail"].lower()


def test_create_compound_invalid_data(client):
    """Test creating a compound with invalid data"""
    invalid_data = {"smiles": "invalid data"}
    response = client.post("/compounds/", json=invalid_data)

    # Check if the API properly validates the request
    assert response.status_code == 400  # Unprocessable Entity
    assert "Invalid SMILES string" in response.json()["detail"]


def test_get_compound_after_creation(client):
    """Test that a created compound can be retrieved"""
    # Create a compound
    create_response = client.post("/compounds/", json=test_compound_data)
    compound_id = create_response.json()["id"]

    # Then try to get it
    get_response = client.get(f"/compounds/{compound_id}")

    # Check response
    assert get_response.status_code == 200
    data = get_response.json()

    # Verify the returned data matches what we created
    assert data["id"] == compound_id
    assert "canonical_smiles" in data
    assert "inchi" in data
    assert "inchikey" in data


def test_list_compounds(client):
    """Test listing all compounds"""

    client.post("/compounds/", json={"smiles": aspirin_smiles})
    client.post("/compounds/", json={"smiles": caffeine_smiles})

    # Get all compounds
    response = client.get("/compounds/")

    # Check response
    assert response.status_code == 200
    data = response.json()

    # Should have at least 2 compounds
    assert len(data) >= 2

    # Check if compounds are in the list (we can't check exact SMILES since they're canonicalized)
    assert len(data) >= 2


def test_create_compound_with_equivalent_smiles(client):
    """Test that creating a compound with equivalent but different SMILES representation fails"""
    # First create a compound with canonical SMILES
    client.post("/compounds/", json={"smiles": aspirin_smiles})

    # Try to create another compound with non-canonical but chemically equivalent SMILES
    response = client.post("/compounds/", json={"smiles": aspirin_smiles_noncanonical})

    # This should fail since they represent the same molecule
    assert response.status_code == 400
    assert "already exists" in response.json()["detail"].lower()


def test_substructure_search(client):
    """Test substructure search using the API endpoint"""
    # Define test molecules
    substructure = "Oc1c(C(O)=O)cccc1"
    mol1 = "CC(Oc1c(C(O)=O)cccc1)=O"  # Should match
    mol2 = "C(Oc1c(C(O)=O)c(C)ccc1)=O"  # Should match
    mol3 = "CC(Oc1ccccc1)=O"  # Should NOT match

    # Create compounds in the database
    response1 = client.post("/compounds/", json={"smiles": mol1})
    response2 = client.post("/compounds/", json={"smiles": mol2})
    response3 = client.post("/compounds/", json={"smiles": mol3})

    # Get compound IDs
    compound1_id = response1.json()["id"]
    compound2_id = response2.json()["id"]
    compound3_id = response3.json()["id"]

    # Test the API endpoint with substructure parameter
    response = client.get(f"/compounds/?substructure={substructure}")
    assert response.status_code == 200
    data = response.json()

    # Extract IDs from API response
    result_ids = [compound["id"] for compound in data]

    # Verify correct results
    assert compound1_id in result_ids, "Molecule 1 should match the substructure"
    assert compound2_id in result_ids, "Molecule 2 should match the substructure"
    assert compound3_id not in result_ids, "Molecule 3 should NOT match the substructure"


def test_create_compound_and_get_hash(client):
    """Test creating a compound and retrieving its hash using the exact search endpoint"""
    compound_data = {
        "smiles": "C[C@@](F)(Cl)c1cc2ccc[nH]c-2n1",
        "is_archived": False,
    }

    # Step 1: Create the compound
    response = client.post("/compounds/", json=compound_data)
    assert response.status_code == 200, f"Unexpected status code: {response.status_code}"
    # compound_id = response.json()["id"]

    # Step 2: Use the /search/compounds/exact endpoint to retrieve the hash_mol
    search_payload = {"query_smiles": compound_data["smiles"]}
    search_response = client.post("/search/compounds/exact", json=search_payload)
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

    compound_ids = []
    for smiles in smiles_list:
        response = client.post("/compounds/", json={"smiles": smiles, "is_archived": False})
        assert response.status_code == 200, f"Failed to create compound for SMILES: {smiles}"
        compound_ids.append(response.json()["id"])

    yield compound_ids

    # Cleanup after tests
    for compound_id in compound_ids:
        client.delete(f"/compounds/{compound_id}")


def test_search_compound_structure_tautomer(client, predefined_compounds):
    """Test tautomer search using the /search/compounds/structure endpoint"""
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
    search_response = client.post("/search/compounds/structure", json=search_payload)
    assert search_response.status_code == 200, f"Unexpected status code: {search_response.status_code}"
    search_results = search_response.json()
    result_ids = [compound["id"] for compound in search_results]

    assert predefined_compounds[0] in result_ids, "Tautomer1 R should match the query"
    assert predefined_compounds[1] in result_ids, "Tautomer1 S should match the query"
    assert predefined_compounds[2] not in result_ids, "Tautomer2 R should NOT match the query"
    assert predefined_compounds[3] not in result_ids, "Tautomer2 S should NOT match the query"


def test_search_compound_structure_stero(client, predefined_compounds):
    """Test stereo search using the /search/compounds/structure endpoint"""
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
    search_response = client.post("/search/compounds/structure", json=search_payload)
    assert search_response.status_code == 200, f"Unexpected status code: {search_response.status_code}"
    search_results = search_response.json()
    result_ids = [compound["id"] for compound in search_results]

    assert predefined_compounds[0] in result_ids, "Tautomer1 R should match the query"
    assert predefined_compounds[2] in result_ids, "Tautomer2 R should match the query"
    assert predefined_compounds[1] not in result_ids, "Tautomer1 S should NOT match the query"
    assert predefined_compounds[3] not in result_ids, "Tautomer2 S should NOT match the query"


def test_search_compound_structure_connectivity(client, predefined_compounds):
    """Test connectivity search using the /search/compounds/structure endpoint"""
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
    search_response = client.post("/search/compounds/structure", json=search_payload)
    assert search_response.status_code == 200, f"Unexpected status code: {search_response.status_code}"
    search_results = search_response.json()
    result_ids = [compound["id"] for compound in search_results]

    assert predefined_compounds[0] in result_ids, "Tautomer1 R should match the query"
    assert predefined_compounds[1] in result_ids, "Tautomer1 S should match the query"
    assert predefined_compounds[2] in result_ids, "Tautomer2 R should match the query"
    assert predefined_compounds[3] in result_ids, "Tautomer2 S should match the query"
