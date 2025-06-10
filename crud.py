import random
import uuid
from sqlalchemy.orm import Session, joinedload
from fastapi import HTTPException
from rdkit import Chem
from typing import List
from sqlalchemy import insert, text
from datetime import datetime, timezone
import models as models
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import main
import enums


from typing import Type, Dict, Any
from rdkit.Chem.RegistrationHash import HashLayer, GetMolHash
from chemistry_utils import standardize_mol, generate_hash_layers, generate_uuid_from_string

# Handle both package imports and direct execution
try:
    # When imported as a package (for tests)
    from . import models
except ImportError:
    # When run directly
    import models


# Compound CRUD operations
def get_compound_by_inchi_key(db: Session, inchikey: str):
    return db.query(models.Compound).filter(models.Compound.inchikey == inchikey).first()


def get_compound_by_hash(db: Session, hash_mol: str):
    """
    Search for compounds in the database by hash_mol.
    """
    if not isinstance(hash_mol, str) or len(hash_mol) != 40:
        raise HTTPException(status_code=400, detail=f"Invalid hash_mol format: {hash_mol}")

    return db.query(models.Compound).filter(models.Compound.hash_mol == hash_mol).all()


def get_compound_by_canonical_smiles(db: Session, canonical_smiles: str):
    return db.query(models.Compound).filter(models.Compound.canonical_smiles == canonical_smiles).first()


def get_compound_by_tautomer(db: Session, tautomer: str):
    return db.query(models.Compound).filter(models.Compound.tautomer == tautomer).all()


def get_compounds(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Compound).options(joinedload(models.Compound.batches)).offset(skip).limit(limit).all()


def create_compound(db: Session, compound: models.CompoundCreate):
    # Create RDKit molecule from SMILES
    mol = Chem.MolFromSmiles(compound.smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")

    # Standardize mol and calculate hashes

    standarized_mol = standardize_mol(mol)
    mol_layers = generate_hash_layers(standarized_mol)
    hash_mol = GetMolHash(mol_layers)
    formula = CalcMolFormula(standarized_mol)
    canonical_smiles = mol_layers[HashLayer.CANONICAL_SMILES]
    hash_canonical_smiles = generate_uuid_from_string(mol_layers[HashLayer.CANONICAL_SMILES])
    hash_tautomer = generate_uuid_from_string(mol_layers[HashLayer.TAUTOMER_HASH])
    hash_no_stereo_smiles = generate_uuid_from_string(mol_layers[HashLayer.NO_STEREO_SMILES])
    hash_no_stereo_tautomer = generate_uuid_from_string(mol_layers[HashLayer.NO_STEREO_TAUTOMER_HASH])
    # sgroup_data = mol_layers[HashLayer.SGROUP_DATA]

    inchi = Chem.MolToInchi(standarized_mol)
    inchikey = Chem.InchiToInchiKey(inchi)

    existing_compound = get_compound_by_hash(db, hash_mol)
    if existing_compound:
        raise HTTPException(status_code=400, detail=f"Compound with hash_mol {hash_mol} already exists")

    # existing_compound = get_compound_by_inchi_key(db, inchikey)
    # if existing_compound:
    #     raise HTTPException(status_code=400, detail=f"Compound with InChIKey {inchikey} already exists")

    # existing_compound = get_compound_by_canonical_smiles(db, canonical_smiles)
    # if existing_compound:
    #     raise HTTPException(status_code=400, detail=f"Compound with canonical SMILES {canonical_smiles} already exists")

    # Create the compound with calculated values
    # NOTE: The following UUID values are placeholders for hash fields
    # and will be replaced with actual computed hash values later.
    db_compound = models.Compound(
        canonical_smiles=canonical_smiles,
        original_molfile=compound.original_molfile,
        inchi=inchi,
        inchikey=inchikey,
        molregno=random.randint(1, 100),
        formula=formula,
        hash_mol=hash_mol,
        hash_tautomer=hash_tautomer,
        hash_canonical_smiles=hash_canonical_smiles,
        hash_no_stereo_smiles=hash_no_stereo_smiles,
        hash_no_stereo_tautomer=hash_no_stereo_tautomer,
        created_at=datetime.now(),
        updated_at=datetime.now(),
        created_by=main.admin_user_id,
        updated_by=main.admin_user_id,
        is_archived=compound.is_archived,
    )

    db.add(db_compound)
    db.commit()
    db.refresh(db_compound)
    return db_compound


def create_compounds_batch(db: Session, smiles_list: List[str]):
    """
    Create multiple compounds from a list of SMILES strings.

    Args:
        db: Database session
        smiles_list: List of SMILES strings

    Returns:
        List of created compounds

    Raises:
        HTTPException: If any SMILES is invalid or already exists
    """
    # Validate all SMILES first
    invalid_smiles = []
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid_smiles.append({"index": i, "smiles": smiles})

    if invalid_smiles:
        raise HTTPException(status_code=400, detail=f"Invalid SMILES strings found: {invalid_smiles}")

    # Calculate InChIKeys for all compounds
    compounds_data = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi)

        compounds_data.append(
            {
                "smiles": smiles,
                "canonical_smiles": canonical_smiles,
                "inchi": inchi,
                "inchikey": inchikey,
            }
        )

    # Check if any compounds already exist
    existing_compounds = []
    for i, compound in enumerate(compounds_data):
        existing = get_compound_by_inchi_key(db, compound["inchikey"])
        if existing:
            existing_compounds.append(
                {
                    "index": i,
                    "smiles": compound["smiles"],
                    "inchikey": compound["inchikey"],
                }
            )

    if existing_compounds:
        raise HTTPException(
            status_code=400, detail=f"Some compounds already exist in the database: {existing_compounds}"
        )

    # Create all compounds
    # NOTE: The following UUID values are placeholders for hash fields
    # and will be replaced with actual computed hash values later.
    created_compounds = []
    for compound_data in compounds_data:
        db_compound = models.Compound(
            canonical_smiles=compound_data["canonical_smiles"],
            inchi=compound_data["inchi"],
            inchikey=compound_data["inchikey"],
            molregno=random.randint(1, 100),
            formula=CalcMolFormula(mol),
            hash_mol=uuid.uuid4(),
            hash_tautomer=uuid.uuid4(),
            hash_canonical_smiles=uuid.uuid4(),
            hash_no_stereo_smiles=uuid.uuid4(),
            hash_no_stereo_tautomer=uuid.uuid4(),
            created_by=main.admin_user_id,
            updated_by=main.admin_user_id,
            created_at=datetime.now(),
            updated_at=datetime.now(),
            is_archived=False,
        )
        db.add(db_compound)
        created_compounds.append(db_compound)

    db.commit()

    # Refresh all compounds to get their IDs
    for compound in created_compounds:
        db.refresh(compound)

    return created_compounds


# def update_compound(db: Session, compound_id: int, compound: schemas.CompoundUpdate):
#     db_compound = db.query(models.Compound).filter(models.Compound.id == compound_id).first()

#     update_data = compound.model_dump(exclude_unset=True)
#     for key, value in update_data.items():
#         setattr(db_compound, key, value)

#     db_compound.updated_at = datetime.now()
#     db.add(db_compound)
#     db.commit()
#     db.refresh(db_compound)
#     return db_compound


# Batch CRUD operations
def get_batch(db: Session, batch_id: int):
    return db.query(models.Batch).filter(models.Batch.id == batch_id).first()


def get_batches(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Batch).offset(skip).limit(limit).all()


def get_batches_by_compound(db: Session, compound_id: int, skip: int = 0, limit: int = 100):
    return db.query(models.Batch).filter(models.Batch.compound_id == compound_id).offset(skip).limit(limit).all()


def create_batch(db: Session, batch: models.BatchBase):
    db_batch = models.Batch(
        compound_id=batch.compound_id,
        notes=batch.notes,
        created_by=main.admin_user_id,
        updated_by=main.admin_user_id,
        created_at=datetime.now(),
        batch_regno=random.randint(1, 100),
    )

    db.add(db_batch)
    db.commit()
    db.refresh(db_batch)
    return db_batch


# def update_batch(db: Session, batch_id: int, batch: schemas.BatchUpdate):
#     db_batch = db.query(models.Batch).filter(models.Batch.id == batch_id).first()

#     update_data = batch.model_dump(exclude_unset=True)
#     for key, value in update_data.items():
#         setattr(db_batch, key, value)

#     db_batch.updated_at = datetime.now()
#     db.add(db_batch)
#     db.commit()
#     db.refresh(db_batch)
#     return db_batch


def delete_batch(db: Session, batch_id: int):
    db_batch = db.query(models.Batch).filter(models.Batch.id == batch_id).first()
    db.delete(db_batch)
    db.commit()
    return db_batch


def create_semantic_type(db: Session, semantic_type: models.SemanticTypeBase):
    db_semantic_type = models.SemanticType(name=semantic_type.name, description=semantic_type.description)

    db.add(db_semantic_type)
    db.commit()
    db.refresh(db_semantic_type)
    return db_semantic_type


# Property CRUD operations
def get_property(db: Session, property_id: int):
    return db.query(models.Property).filter(models.Property.id == property_id).first()


def get_properties(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Property).offset(skip).limit(limit).all()


def create_property(db: Session, property: models.PropertyBase):
    db_property = models.Property(
        name=property.name,
        value_type=property.value_type,
        property_class=property.property_class,
        unit=property.unit,
        semantic_type_id=property.semantic_type_id,
        scope=property.scope,
        created_by=main.admin_user_id,
        updated_by=main.admin_user_id,
    )

    db.add(db_property)
    db.commit()
    db.refresh(db_property)
    return db_property


# def update_property(db: Session, property_id: int, property: schemas.PropertyUpdate):
#     db_property = db.query(models.Property).filter(models.Property.id == property_id).first()

#     update_data = property.model_dump(exclude_unset=True)
#     for key, value in update_data.items():
#         setattr(db_property, key, value)

#     db_property.updated_at = datetime.now()
#     db.add(db_property)
#     db.commit()
#     db.refresh(db_property)
#     return db_property


def delete_property(db: Session, property_id: int):
    db_property = db.query(models.Property).filter(models.Property.id == property_id).first()
    db.delete(db_property)
    db.commit()
    return db_property


# AssayType CRUD operations
def get_assay_type(db: Session, assay_type_id: int):
    return db.query(models.AssayType).filter(models.AssayType.id == assay_type_id).first()


def get_assay_types(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.AssayType).offset(skip).limit(limit).all()


def create_assay_type(db: Session, assay_type: models.AssayTypeCreate):
    # Create the assay type
    current_time = datetime.now(timezone.utc)
    db_assay_type = models.AssayType(
        name=assay_type.name,
        description=assay_type.description,
        created_at=current_time,
        updated_at=current_time,
        created_by=main.admin_user_id,
        updated_by=main.admin_user_id,
    )
    db.add(db_assay_type)
    db.commit()
    db.refresh(db_assay_type)

    # Add properties to the assay type if provided
    if assay_type.property_ids:
        properties = db.query(models.Property).filter(models.Property.id.in_(assay_type.property_ids)).all()
        for prop in properties:
            db_assay_type.properties.append(prop)

    # Add property requirements if provided
    for req in assay_type.property_requirements:
        # Skip if property not in the assay type's properties
        if req["property_id"] not in assay_type.property_ids:
            continue

        property_req = models.AssayTypeProperty(
            assay_type_id=db_assay_type.id,
            property_id=req["property_id"],
            required=req.get("required", False),
        )
        db.add(property_req)

    # Add property details if provided
    for detail in assay_type.property_details:
        # Skip if property not in the assay type's properties
        if detail["property_id"] not in assay_type.property_ids:
            continue

        # Get the property to determine its value type
        property = db.query(models.Property).filter(models.Property.id == detail["property_id"]).first()
        if not property:
            continue

        property_detail = models.AssayTypeDetail(
            assay_type_id=db_assay_type.id,
            property_id=detail["property_id"],
        )

        # Set the appropriate value based on property type
        if property.value_type == "datetime" and "value_datetime" in detail:
            property_detail.value_datetime = detail["value_datetime"]
        elif property.value_type in ("int", "double") and "value_num" in detail:
            property_detail.value_num = detail["value_num"]
        elif property.value_type == "string" and "value_string" in detail:
            property_detail.value_string = detail["value_string"]

        db.add(property_detail)

    db.commit()
    db.refresh(db_assay_type)
    return db_assay_type


# Assay CRUD operations
def get_assay(db: Session, assay_id: int):
    # Get the assay
    assay = db.query(models.Assay).filter(models.Assay.id == assay_id).first()

    if assay:
        # Get properties associated with this assay through assay details
        assay_details = db.query(models.AssayDetail).filter(models.AssayDetail.assay_id == assay_id).all()
        property_ids = [detail.property_id for detail in assay_details]

        # If no properties from assay details, get them from the assay type
        if not property_ids and assay.assay_type:
            assay_type_properties = (
                db.query(models.AssayTypeProperty)
                .filter(models.AssayTypeProperty.assay_type_id == assay.assay_type_id)
                .all()
            )
            property_ids = [prop.property_id for prop in assay_type_properties]

        # Get the property objects
        if property_ids:
            properties = db.query(models.Property).filter(models.Property.id.in_(property_ids)).all()
            # Add properties to assay
            assay.properties = properties
        else:
            assay.properties = []

    return assay


def get_assays(db: Session, skip: int = 0, limit: int = 100):
    # Get assays with pagination
    assays = db.query(models.Assay).offset(skip).limit(limit).all()

    # For each assay, add its properties
    for assay in assays:
        # Get properties associated with this assay through assay details
        assay_details = db.query(models.AssayDetail).filter(models.AssayDetail.assay_id == assay.id).all()
        property_ids = [detail.property_id for detail in assay_details]

        # If no properties from assay details, get them from the assay type
        if not property_ids and assay.assay_type:
            assay_type_properties = (
                db.query(models.AssayTypeProperty)
                .filter(models.AssayTypeProperty.assay_type_id == assay.assay_type_id)
                .all()
            )
            property_ids = [prop.property_id for prop in assay_type_properties]

        # Get the property objects
        if property_ids:
            properties = db.query(models.Property).filter(models.Property.id.in_(property_ids)).all()
            # Add properties to assay
            assay.properties = properties
        else:
            assay.properties = []

    return assays


def create_assay(db: Session, assay: models.AssayCreate):
    # Create the assay
    db_assay = models.Assay(
        name=assay.name,
        description=assay.description,
        assay_type_id=assay.assay_type_id,
        created_by=main.admin_user_id,
        updated_by=main.admin_user_id,
        created_at=datetime.now(),
    )
    db.add(db_assay)
    db.commit()
    db.refresh(db_assay)

    # Initialize empty properties list
    db_assay.properties = []

    # Add assay details for properties if provided
    if assay.property_ids:
        # Get the properties
        properties = db.query(models.Property).filter(models.Property.id.in_(assay.property_ids)).all()

        # Get the assay type to check which properties are expected
        assay_type = db.query(models.AssayType).filter(models.AssayType.id == assay.assay_type_id).first()

        # Only add properties that are part of the assay type
        valid_property_ids = set()
        if assay_type and assay_type.properties:
            valid_property_ids = {prop.id for prop in assay_type.properties}

        # Collect valid properties to add to assay.properties
        valid_properties = []

        # Create assay_details entries for each valid property
        for prop in properties:
            if prop.id in valid_property_ids or not valid_property_ids:  # If valid_property_ids is empty, accept all
                # Initialize empty detail
                assay_detail = models.AssayDetail(
                    assay_id=db_assay.id,
                    property_id=prop.id,
                )
                db.add(assay_detail)
                valid_properties.append(prop)

        db.commit()
        db.refresh(db_assay)

        # Set the properties list
        db_assay.properties = valid_properties

    return db_assay


def get_compounds_ex(db: Session, query_params: models.CompoundQueryParams):
    """
    Get compounds with optional filtering parameters.

    Args:
        db: Database session
        query_params: Query parameters including substructure, skip, and limit

    Returns:
        List of compounds matching the query parameters
    """
    # If substructure is provided, use substructure search
    if query_params.substructure:
        # Validate the substructure SMILES
        mol = Chem.MolFromSmiles(query_params.substructure)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid substructure SMILES string")

        # SQL query using the RDKit cartridge substructure operator '@>'
        sql = text(f"""
            SELECT c.* FROM {models.DB_SCHEMA}.compounds c
            JOIN rdk.mols ON rdk.mols.id = c.id
            WHERE rdk.mols.m@>'{query_params.substructure}'
            ORDER BY c.id
            OFFSET :skip LIMIT :limit
        """)

        result = db.execute(sql, {"skip": query_params.skip, "limit": query_params.limit})
        compounds = []
        for row in result:
            compound = models.Compound()
            for column, value in row._mapping.items():
                setattr(compound, column, value)
            compounds.append(compound)

        return compounds

    else:
        # If no substructure provided, use regular get_compounds function
        return get_compounds(db, skip=query_params.skip, limit=query_params.limit)


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


# AssayResult CRUD operations
def get_assay_result(db: Session, assay_result_id: int):
    return db.query(models.AssayResult).filter(models.AssayResult.id == assay_result_id).first()


def get_assay_results(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.AssayResult).offset(skip).limit(limit).all()


def get_batch_assay_results(db: Session, batch_id: int):
    """Get all assay results for a specific batch"""
    results = db.query(models.AssayResult).filter(models.AssayResult.batch_id == batch_id).all()

    # Group results by assay_id
    grouped_results = {}
    for result in results:
        assay_id = result.assay_id
        if assay_id not in grouped_results:
            # Get the assay name
            assay = db.query(models.Assay).filter(models.Assay.id == assay_id).first()
            assay_name = assay.name if assay else "Unknown Assay"

            grouped_results[assay_id] = {
                "assay_id": assay_id,
                "batch_id": batch_id,
                "assay_name": assay_name,
                "measurements": {},
            }

        # Get property name and type
        property = db.query(models.Property).filter(models.Property.id == result.property_id).first()
        property_name = property.name if property else f"Property-{result.property_id}"
        property_type = property.value_type if property else "double"

        # Get value based on property type
        value = None
        if property_type in ("int", "double"):
            value = result.value_num
        elif property_type == "string":
            value = result.value_string
        elif property_type == "bool":
            value = result.value_bool

        # If we have a qualifier other than "=" (0), include it in the result
        if result.value_qualifier != 0:
            grouped_results[assay_id]["measurements"][property_name] = {
                "qualifier": result.value_qualifier,
                "value": value,
            }
        else:
            grouped_results[assay_id]["measurements"][property_name] = value

    return list(grouped_results.values())


def create_assay_result(db: Session, assay_result: models.AssayResultBase):
    """Create a single assay result entry for a specific property"""
    # Validate assay exists
    assay = db.query(models.Assay).filter(models.Assay.id == assay_result.assay_id).first()
    if assay is None:
        raise HTTPException(status_code=404, detail=f"Assay with ID {assay_result.assay_id} not found")

    # Validate batch exists
    batch = db.query(models.Batch).filter(models.Batch.id == assay_result.batch_id).first()
    if batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {assay_result.batch_id} not found")

    # Validate property exists and get its type
    property = db.query(models.Property).filter(models.Property.id == assay_result.property_id).first()
    if property is None:
        raise HTTPException(status_code=404, detail=f"Property with ID {assay_result.property_id} not found")

    # Check if the property is associated with the assay through assay_details
    assay_detail = (
        db.query(models.AssayDetail)
        .filter(
            models.AssayDetail.assay_id == assay_result.assay_id,
            models.AssayDetail.property_id == assay_result.property_id,
        )
        .first()
    )

    # If not found in assay_details, check if the property is associated with the assay's assay_type
    if not assay_detail:
        assay_type = db.query(models.AssayType).filter(models.AssayType.id == assay.assay_type_id).first()
        if assay_type:
            # Check if property is in assay type properties
            assay_type_property = (
                db.query(models.AssayTypeProperty)
                .filter(
                    models.AssayTypeProperty.assay_type_id == assay_type.id,
                    models.AssayTypeProperty.property_id == assay_result.property_id,
                )
                .first()
            )

            if not assay_type_property:
                raise HTTPException(
                    status_code=400,
                    detail=f"Property with ID {assay_result.property_id} is not associated with assay type '{assay_type.name}'",
                )

    # Create the assay result with the appropriate value field
    db_assay_result = models.AssayResult(
        assay_id=assay_result.assay_id,
        batch_id=assay_result.batch_id,
        property_id=assay_result.property_id,
        value_qualifier=assay_result.value_qualifier,
    )

    # Set the appropriate value based on property type
    if property.value_type in ("int", "double") and assay_result.value_num is not None:
        db_assay_result.value_num = assay_result.value_num
    elif property.value_type == "string" and assay_result.value_string is not None:
        db_assay_result.value_string = assay_result.value_string
    elif property.value_type == "bool" and assay_result.value_bool is not None:
        db_assay_result.value_bool = assay_result.value_bool
    else:
        raise HTTPException(
            status_code=400,
            detail=f"Property value of type {property.value_type} missing or not supported for AssayResult",
        )

    db.add(db_assay_result)
    db.commit()
    db.refresh(db_assay_result)
    return db_assay_result


def create_batch_assay_results(db: Session, batch_results: models.BatchAssayResultsCreate):
    """Create multiple assay results for a batch in a single transaction"""
    # Validate assay exists
    assay = db.query(models.Assay).filter(models.Assay.id == batch_results.assay_id).first()
    if assay is None:
        raise HTTPException(status_code=404, detail=f"Assay with ID {batch_results.assay_id} not found")

    # Validate batch exists
    batch = db.query(models.Batch).filter(models.Batch.id == batch_results.batch_id).first()
    if batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {batch_results.batch_id} not found")

    # Get the assay type
    assay_type = db.query(models.AssayType).filter(models.AssayType.id == assay.assay_type_id).first()
    if assay_type is None:
        raise HTTPException(status_code=404, detail=f"AssayType for Assay with ID {batch_results.assay_id} not found")

    # Get all properties from the assay type
    properties = {}
    for prop in assay_type.properties:
        properties[prop.name] = {
            "id": prop.id,
            "value_type": prop.value_type,
        }

    # Validate all property names in measurements
    for prop_name in batch_results.measurements.keys():
        if prop_name not in properties:
            raise HTTPException(
                status_code=400, detail=f"Property '{prop_name}' is not associated with assay type '{assay_type.name}'"
            )

    # Create assay results for each property
    created_results = []
    processed_measurements = {}

    for prop_name, measurement in batch_results.measurements.items():
        property_id = properties[prop_name]["id"]
        property_type = properties[prop_name]["value_type"]

        # Prepare the assay result with common fields
        db_assay_result = models.AssayResult(
            assay_id=batch_results.assay_id,
            batch_id=batch_results.batch_id,
            property_id=property_id,
            value_qualifier=0,  # Default to equals
        )

        # Handle complex measurement (dict with qualifier and value)
        if isinstance(measurement, dict) and "value" in measurement:
            # Set qualifier if provided
            if "qualifier" in measurement:
                db_assay_result.value_qualifier = measurement["qualifier"]

            value = measurement["value"]

            # Set the value based on property type
            if property_type in ("int", "double"):
                db_assay_result.value_num = float(value)
                processed_measurements[prop_name] = {
                    "qualifier": db_assay_result.value_qualifier,
                    "value": value,
                }
            elif property_type == "string":
                db_assay_result.value_string = str(value)
                processed_measurements[prop_name] = {
                    "qualifier": db_assay_result.value_qualifier,
                    "value": value,
                }
            elif property_type == "bool":
                db_assay_result.value_bool = bool(value)
                processed_measurements[prop_name] = {
                    "qualifier": db_assay_result.value_qualifier,
                    "value": value,
                }
        # Handle simple value (backward compatibility)
        else:
            if property_type in ("int", "double"):
                db_assay_result.value_num = float(measurement)
                processed_measurements[prop_name] = measurement
            elif property_type == "string":
                db_assay_result.value_string = str(measurement)
                processed_measurements[prop_name] = measurement
            elif property_type == "bool":
                db_assay_result.value_bool = bool(measurement)
                processed_measurements[prop_name] = measurement

        db.add(db_assay_result)
        created_results.append(db_assay_result)

    db.commit()

    # Refresh all results to get their IDs
    for result in created_results:
        db.refresh(result)

    # Return in grouped format
    return {
        "assay_id": batch_results.assay_id,
        "batch_id": batch_results.batch_id,
        "assay_name": assay.name,
        "measurements": processed_measurements,
    }


# BatchDetail CRUD operations
def get_batch_detail(db: Session, batch_detail_id: int):
    return db.query(models.BatchDetail).filter(models.BatchDetail.id == batch_detail_id).first()


def get_batch_details_by_batch(db: Session, batch_id: int, skip: int = 0, limit: int = 100):
    return db.query(models.BatchDetail).filter(models.BatchDetail.batch_id == batch_id).offset(skip).limit(limit).all()


def create_batch_detail(db: Session, batch_detail: models.BatchDetailBase):
    # Get the property to determine its value type
    property = db.query(models.Property).filter(models.Property.id == batch_detail.property_id).first()
    if not property:
        raise HTTPException(status_code=404, detail=f"Property with ID {batch_detail.property_id} not found")

    db_batch_detail = models.BatchDetail(
        batch_id=batch_detail.batch_id,
        property_id=batch_detail.property_id,
        value_qualifier=batch_detail.value_qualifier,
        created_by=main.admin_user_id,
        updated_by=main.admin_user_id,
    )

    # Set the value based on the property type
    if batch_detail.value_datetime is not None and property.value_type == "datetime":
        db_batch_detail.value_datetime = batch_detail.value_datetime
    elif batch_detail.value_num is not None and property.value_type in ("int", "double"):
        db_batch_detail.value_num = batch_detail.value_num
    elif batch_detail.value_string is not None and property.value_type == "string":
        db_batch_detail.value_string = batch_detail.value_string

    db.add(db_batch_detail)
    db.commit()
    db.refresh(db_batch_detail)
    return db_batch_detail


def bulk_create_if_not_exists(
    db: Session,
    model_cls: Type,
    base_model_cls: Type,
    items: List[Any],
    *,
    name_attr: str = "name",
    validate: bool = True,
) -> List[Dict[str, Any]]:
    input_names = [getattr(item, name_attr) for item in items]
    existing_names = {
        name
        for (name,) in db.query(getattr(model_cls, name_attr))
        .filter(getattr(model_cls, name_attr).in_(input_names))
        .all()
    }

    to_insert = []
    for item in items:
        item_name = getattr(item, name_attr)
        if item_name not in existing_names:
            validated = base_model_cls.model_validate(item) if validate else item
            data = validated.model_dump()
            data.update(
                {
                    "created_by": main.admin_user_id,
                    "updated_by": main.admin_user_id,
                }
            )
            to_insert.append(data)

    if not to_insert:
        return []

    stmt = insert(model_cls).values(to_insert).returning(model_cls)
    result = db.execute(stmt).fetchall()
    db.commit()
    return [model_cls.model_validate(row[0]) for row in result]


def create_synonym_types(db: Session, synonym_types: list[models.SynonymTypeBase]) -> list[dict]:
    return bulk_create_if_not_exists(db, models.SynonymType, models.SynonymTypeBase, synonym_types)


def get_properties_by_scope(scope: enums.ScopeClass, db: Session) -> List[models.Property]:
    return db.query(models.Property).filter(models.Property.scope == scope).all()


def get_synonyms_by_level(level: enums.SynonymLevel, db: Session) -> List[models.SynonymType]:
    return db.query(models.SynonymType).filter(models.SynonymType.synonym_level == level).all()


def create_properties(db: Session, properties: list[models.PropertyBase]) -> list[dict]:
    return bulk_create_if_not_exists(db, models.Property, models.PropertyBase, properties)


def enrich_addition(add: models.AdditionBase) -> models.AdditionBase:
    smiles, molfile, formula, mw = add.smiles, add.molfile, add.formula, add.molecular_weight

    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if not mol and molfile:
        try:
            mol = Chem.MolFromMolBlock(molfile, sanitize=True)
        except Exception:
            mol = None

    if not mol:
        return add

    add.smiles = smiles or Chem.MolToSmiles(mol)
    add.molfile = molfile or Chem.MolToMolBlock(mol)
    add.formula = formula or rdMolDescriptors.CalcMolFormula(mol)
    add.molecular_weight = mw or Descriptors.MolWt(mol)

    return add


def create_additions(db: Session, additions: list[models.AdditionBase]) -> list[dict]:
    enriched_additions = [enrich_addition(add) for add in additions]
    return bulk_create_if_not_exists(db, models.Addition, models.AdditionBase, enriched_additions, validate=False)


def get_additions(db: Session, role: enums.AdditionsRole | None = None) -> List[models.AdditionBase]:
    query = db.query(models.Addition)
    if role is None:
        return query.all()
    return query.filter_by(role=role).all()


def get_addition_by_id(db: Session, addition_id: int) -> models.Addition:
    db_addition = db.get(models.Addition, addition_id)
    if db_addition is None:
        raise HTTPException(status_code=404, detail="Addition not found")
    return db_addition


def update_addition_by_id(db: Session, addition_id: int, addition_update: models.AdditionUpdate):
    db_addition = get_addition_by_id(db, addition_id=addition_id)
    update_data = addition_update.dict(exclude_unset=True)
    for key, value in update_data.items():
        setattr(db_addition, key, value)

    db_addition.updated_at = datetime.now()
    db.add(db_addition)
    db.commit()
    db.refresh(db_addition)
    return db_addition


def delete_addition_by_id(db: Session, addition_id: int):
    db_addition = get_addition_by_id(db, addition_id=addition_id)
    db_addition.deleted_at = datetime.now()
    db_addition.is_active = False
    db.commit()
    db.refresh(db_addition)
    return db_addition


def read_compounds(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Compound).offset(skip).limit(limit).all()


def get_compound_by_id(db: Session, compound_id: int):
    return db.query(models.Compound).filter(models.Compound.id == compound_id).first()


def delete_compound(db: Session, compound_id: int):
    db_compound = db.get(models.Compound, compound_id)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")
    db_compound.deleted_at = datetime.now()

    db.commit()
    db.refresh(db_compound)
    return db_compound


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
