from sqlalchemy.orm import Session
from fastapi import HTTPException
from rdkit import Chem
from typing import List
from app import models

from typing import Dict, Any
from rdkit.Chem.RegistrationHash import HashLayer
from app.utils.chemistry_utils import standardize_mol, generate_hash_layers, generate_uuid_from_string


def search_compounds_exact(query_smiles: str, search_parameters: models.ExactSearchParameters, db: Session):
    """
    Perform an exact search for compounds.

    - **query_smiles**: The SMILES string to search against.
    - **search_parameters**: Parameters for the exact search.
    """
    if not search_parameters:
        raise HTTPException(status_code=400, detail="Search parameters are required for exact search")
    exact_params = models.ExactSearchParameters(**search_parameters)
    return db.query(db=db, query_smiles=query_smiles, fields=exact_params.field)


def search_compounds_substructure(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Placeholder for substructure search.
    """
    raise NotImplementedError("Substructure search is not implemented yet.")


def search_compounds_similarity(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Placeholder for similarity search.
    """
    raise NotImplementedError("Similarity search is not implemented yet.")


def get_standardized_mol_and_layers(query_smiles: str, http_errors: bool = False) -> dict:
    missing_msg = "Query SMILES is required"
    invalid_msg = "Invalid SMILES string"

    if not query_smiles:
        if http_errors:
            raise HTTPException(status_code=400, detail=missing_msg)
        else:
            raise ValueError(missing_msg)

    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        if http_errors:
            raise HTTPException(status_code=400, detail=invalid_msg)
        else:
            raise ValueError(invalid_msg)

    standardized_mol = standardize_mol(mol)
    mol_layers = generate_hash_layers(standardized_mol)
    return mol_layers


def search_compounds_by_hash(
    db: Session, query_smiles: str, hash_layer: Any, hash_attr_name: str
) -> List[models.Compound]:
    mol_layers = get_standardized_mol_and_layers(query_smiles, http_errors=True)
    hash_value = generate_uuid_from_string(mol_layers[hash_layer])
    return db.query(models.Compound).filter(getattr(models.Compound, hash_attr_name) == hash_value).all()


def search_compounds_tautomer(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Perform a tautomer search for compounds. For a given molecule, find all compounds that have the same tautomer hash regardless of their stereochemistry.
    """
    return search_compounds_by_hash(db, query_smiles, HashLayer.TAUTOMER_HASH, "hash_tautomer")


def search_compounds_stereo(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Perform a stereo search for compounds regardless tautomeric state."""
    return search_compounds_by_hash(db, query_smiles, HashLayer.NO_STEREO_SMILES, "hash_no_stereo_smiles")


def search_compounds_connectivity(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Perform a connectivity search for compounds.
    This will find compounds that have the same connectivity regardless of their stereochemistry or tautomeric state.
    """
    return search_compounds_by_hash(db, query_smiles, HashLayer.NO_STEREO_TAUTOMER_HASH, "hash_no_stereo_tautomer")
