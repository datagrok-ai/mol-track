"""
Tests for advanced search functionality
"""

import pytest
from fastapi.testclient import TestClient
from models import AtomicCondition, CompareOp

def test_search_compounds_empty_filter(client, preload_black_data):
    search_request = {
        "output": ["compounds.canonical_smiles", "compounds.inchikey", "compounds.formula"],
        "filter": None  # No filter - should return all compounds
    }
    response = client.post("/search/compounds", json=search_request)
    assert response.status_code == 200
    data = response.json()
    assert len(data["data"]) > 0  # Should return some compounds

def test_search_batches_empty_filter(client, preload_black_data):
    """Test basic batch search without filter to see if the endpoint works"""
    search_request = {
        "output": ["batches.details.batch_corporate_id"],
        "filter": None  # No filter - should return all batches
    }
    response = client.post("/search/batches", json=search_request)
    assert response.status_code == 200
    print(response.text)
    data = response.json()
    assert len(data["data"]) > 0  # Should return some batches

def test_search_batches_by_batch_corporate_id(client, preload_black_data):
    search_request = {
        "output": ["batches.details.batch_corporate_id"],
        "filter": {
            "field": "batches.details.batch_corporate_id",
            "operator": "=",
            "value": "EPA-001-001"
        }
    }
    response = client.post("/search/batches", json=search_request)
    print(f"Response status: {response.status_code}")
    print(f"Response body: {response.text}")
    assert response.status_code == 200
    data = response.json()
    assert len(data["data"]) == 1
    # The field resolver converts dots to underscores in column names
    assert data["data"][0]["batches_details_batch_corporate_id"] == "EPA-001-001"


def test_search_assays_no_filter(client, preload_black_data):
    search_request = {
        "output": ["assay_results.details.Mean HTC recovery"],
        "filter": None
    }
    response = client.post("/search/assay-results", json=search_request)
    print(f"Response status: {response.status_code}")
    print(f"Response body: {response.text}")
    assert response.status_code == 200
    data = response.json()
    assert len(data["data"]) == 238   # Hmmm, returns only 165 but there are 238 assay results in the file


def test_search_assays_by_species(client, preload_black_data):
    # First, let's see what properties are available by searching with no filter
    debug_request = {
        "output": ["assay_results.details.Mean HTC recovery"],
        "filter": None
    }
    debug_response = client.post("/search/explain", json=debug_request)
    print(f"Debug response: {debug_response.text}")
    
    # Now create the filter for Trout species and concentration < 0.2
    search_request = {
        "output": ["assay_results.details.Mean HTC recovery"],
        "filter": {
          "field": "assay_run_details.Cell Species",
          "operator": "=",
          "value": "Trout"
        }
    }
    response = client.post("/search/assay-results", json=search_request)
    print(f"Response status: {response.status_code}")
    print(f"Response body: {response.text}")
    assert response.status_code == 200
    data = response.json()
    print(f"Found {len(data['data'])} results")
    assert len(data["data"]) == 80


def test_search_assays_by_species_and_concentration(client, preload_black_data):
    # First, let's see what properties are available by searching with no filter
    debug_request = {
        "output": ["assay_results.details.Mean HTC recovery"],
        "filter": None
    }
    debug_response = client.post("/search/explain", json=debug_request)
    print(f"Debug response: {debug_response.text}")
    
    # Now create the filter for Trout species and concentration < 0.2
    search_request = {
        "output": ["assay_results.details.Mean HTC recovery"],
        "filter": {
            "operator": "AND",
            "conditions": [
                {
                    "field": "assay_results.details.Species",
                    "operator": "=",
                    "value": "Trout"
                },
                {
                    "field": "assay_results.details.Dosed Concentration (µM)",
                    "operator": "<",
                    "value": 0.2
                }
            ]
        }
    }
    response = client.post("/search/assay-results", json=search_request)
    print(f"Response status: {response.status_code}")
    print(f"Response body: {response.text}")
    assert response.status_code == 200
    data = response.json()
    print(f"Found {len(data['data'])} results")
    assert len(data["data"]) == 40  # Should find some results


def test_search_compounds_simple(client, preload_black_data):
    """
    Test the /search/compounds endpoint with a simple compound equals query.
    
    This test:
    1. Uses the preload_black_data fixture to load test data
    2. Searches for a specific compound by canonical SMILES
    3. Verifies the search returns the expected result
    """
    # Test data: search for "1,3-Benzenedicarboxylic acid" which has SMILES "C1=CC(=CC(=C1)C(=O)O)C(=O)O"
    search_request = {
        "output": ["compounds.canonical_smiles", "compounds.inchikey", "compounds.formula"],
        "filter": {
            "field": "compounds.canonical_smiles",
            "operator": "=",
            "value": "C1=CC(=CC(=C1)C(=O)O)C(=O)O"
        }
    }
    
    # Make the search request
    response = client.post("/search/compounds", json=search_request)
    
    assert response.status_code == 200
    data = response.json()
    
    # Verify response structure
    assert "data" in data
    assert "total_count" in data
    assert "level" in data
    assert "columns" in data
    assert data["level"] == "compounds"
    assert data["total_count"] == 1
    
    # Verify columns
    expected_columns = ["compounds_canonical_smiles", "compounds_inchikey", "compounds_formula"]
    assert data["columns"] == expected_columns
    
    # Verify data
    assert len(data["data"]) == 1
    result = data["data"][0]
    assert result["compounds_canonical_smiles"] == "C1=CC(=CC(=C1)C(=O)O)C(=O)O"
    assert "compounds_inchikey" in result
    assert "compounds_formula" in result
