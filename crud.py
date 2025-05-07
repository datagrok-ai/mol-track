from sqlalchemy.orm import Session, joinedload
from fastapi import HTTPException
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from typing import List, Dict, Any, Optional, Union
from sqlalchemy import text
from datetime import datetime, timezone
import yaml
from chemistry_utils import standardize_mol,generate_molhash
from rdkit.Chem.RegistrationHash import HashLayer 

# Handle both package imports and direct execution
try:
    # When imported as a package (for tests)
    from . import models, schemas
except ImportError:
    # When run directly
    import models, schemas

# Compound CRUD operations
def get_compound(db: Session, compound_id: int):
    return db.query(models.Compound).filter(models.Compound.id == compound_id).first()

def get_compound_by_inchi_key(db: Session, inchikey: str):
    return db.query(models.Compound).filter(models.Compound.inchikey == inchikey).first()

def get_compound_by_molhash(db: Session, molhash: str):
    return db.query(models.Compound).filter(models.Compound.molhash == molhash).first()

def get_compound_by_canonical_smiles(db: Session, canonical_smiles: str):
    return db.query(models.Compound).filter(models.Compound.canonical_smiles == canonical_smiles).first()

def get_compounds(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Compound).options(joinedload(models.Compound.batches)).offset(skip).limit(limit).all()

def create_compound(db: Session, compound: schemas.CompoundCreate):
    # Create RDKit molecule from SMILES
    mol = Chem.MolFromSmiles(compound.smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    # Standardize mol and calculate hashes

    standarized_mol = standardize_mol(mol)
    molhash, mol_layers = generate_molhash(standarized_mol)


    formula = mol_layers[HashLayer.FORMULA]
    canonical_smiles = mol_layers[HashLayer.CANONICAL_SMILES]
    tautomer = mol_layers[HashLayer.TAUTOMER_HASH]
    no_stereo_smiles = mol_layers[HashLayer.NO_STEREO_SMILES]
    no_stereo_tautomer = mol_layers[HashLayer.NO_STEREO_TAUTOMER_HASH]
    sgroup_data = mol_layers[HashLayer.NO_STEREO_TAUTOMER_HASH]

    
    # canonical_smiles = Chem.MolToSmiles(standarized_mol, isomericSmiles=True, canonical=True)
    # Calculate canonical SMILES, InChI, and InChIKey
    inchi = Chem.MolToInchi(standarized_mol)
    inchikey = Chem.InchiToInchiKey(inchi)
    
    
    # Check if compound with this InChIKey already exists

    existing_compound = get_compound_by_molhash(db, molhash)
    if existing_compound:
        raise HTTPException(status_code=400, detail=f"Compound with molhash {molhash} already exists")

    existing_compound = get_compound_by_inchi_key(db, inchikey)
    if existing_compound:
        raise HTTPException(status_code=400, detail=f"Compound with InChIKey {inchikey} already exists")

    existing_compound = get_compound_by_canonical_smiles(db, canonical_smiles)
    if existing_compound:
        raise HTTPException(status_code=400, detail=f"Compound with canonical SMILES {canonical_smiles} already exists")
    


    # Create the compound with calculated values
    db_compound = models.Compound(
        canonical_smiles=canonical_smiles,
        original_molfile=compound.original_molfile,
        molhash = molhash,
        formula = formula,
        tautomer = tautomer,
        no_stereo_smiles = no_stereo_smiles,
        no_stereo_tautomer = no_stereo_tautomer,
        sgroup_data = sgroup_data, 
        inchi=inchi,
        inchikey=inchikey,
        created_at=datetime.now(),
        updated_at=datetime.now(),
        is_archived=compound.is_archived
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
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid SMILES strings found: {invalid_smiles}"
        )
    
    # Calculate InChIKeys for all compounds
    compounds_data = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi)
        
        compounds_data.append({
            "smiles": smiles,
            "canonical_smiles": canonical_smiles,
            "inchi": inchi,
            "inchikey": inchikey
        })
    
    # Check if any compounds already exist
    existing_compounds = []
    for i, compound in enumerate(compounds_data):
        existing = get_compound_by_inchi_key(db, compound["inchikey"])
        if existing:
            existing_compounds.append({
                "index": i, 
                "smiles": compound["smiles"], 
                "inchikey": compound["inchikey"]
            })
    
    if existing_compounds:
        raise HTTPException(
            status_code=400, 
            detail=f"Some compounds already exist in the database: {existing_compounds}"
        )
    
    # Create all compounds
    created_compounds = []
    for compound_data in compounds_data:
        db_compound = models.Compound(
            canonical_smiles=compound_data["canonical_smiles"],
            inchi=compound_data["inchi"],
            inchikey=compound_data["inchikey"],
            created_at=datetime.now(),
            updated_at=datetime.now(),
            is_archived=False
        )
        db.add(db_compound)
        created_compounds.append(db_compound)
    
    db.commit()
    
    # Refresh all compounds to get their IDs
    for compound in created_compounds:
        db.refresh(compound)
    
    return created_compounds

def update_compound(db: Session, compound_id: int, compound: schemas.CompoundUpdate):
    db_compound = db.query(models.Compound).filter(models.Compound.id == compound_id).first()
    
    update_data = compound.model_dump(exclude_unset=True)
    for key, value in update_data.items():
        setattr(db_compound, key, value)
    
    db_compound.updated_at = datetime.now()
    db.add(db_compound)
    db.commit()
    db.refresh(db_compound)
    return db_compound

def delete_compound(db: Session, compound_id: int):
    db_compound = db.query(models.Compound).filter(models.Compound.id == compound_id).first()
    db.delete(db_compound)
    db.commit()
    return db_compound

# Batch CRUD operations
def get_batch(db: Session, batch_id: int):
    return db.query(models.Batch).filter(models.Batch.id == batch_id).first()

def get_batches(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Batch).offset(skip).limit(limit).all()

def get_batches_by_compound(db: Session, compound_id: int, skip: int = 0, limit: int = 100):
    return db.query(models.Batch).filter(models.Batch.compound_id == compound_id).offset(skip).limit(limit).all()

def create_batch(db: Session, batch: schemas.BatchCreate):
    db_batch = models.Batch(
        compound_id=batch.compound_id,
        batch_number=batch.batch_number,
        amount=batch.amount,
        amount_unit=batch.amount_unit,
        purity=batch.purity,
        notes=batch.notes,
        expiry_date=batch.expiry_date,
        created_at=datetime.now()
    )
    db.add(db_batch)
    db.commit()
    db.refresh(db_batch)
    return db_batch

def update_batch(db: Session, batch_id: int, batch: schemas.BatchUpdate):
    db_batch = db.query(models.Batch).filter(models.Batch.id == batch_id).first()
    
    update_data = batch.model_dump(exclude_unset=True)
    for key, value in update_data.items():
        setattr(db_batch, key, value)
    
    db_batch.updated_at = datetime.now()
    db.add(db_batch)
    db.commit()
    db.refresh(db_batch)
    return db_batch

def delete_batch(db: Session, batch_id: int):
    db_batch = db.query(models.Batch).filter(models.Batch.id == batch_id).first()
    db.delete(db_batch)
    db.commit()
    return db_batch

# Property CRUD operations
def get_property(db: Session, property_id: int):
    return db.query(models.Property).filter(models.Property.id == property_id).first()

def get_properties(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Property).offset(skip).limit(limit).all()

def create_property(db: Session, property: schemas.PropertyCreate):
    db_property = models.Property(
        name=property.name,
        value_type=property.value_type,
        property_class=property.property_class,
        unit=property.unit
    )
    db.add(db_property)
    db.commit()
    db.refresh(db_property)
    return db_property

def update_property(db: Session, property_id: int, property: schemas.PropertyUpdate):
    db_property = db.query(models.Property).filter(models.Property.id == property_id).first()
    
    update_data = property.model_dump(exclude_unset=True)
    for key, value in update_data.items():
        setattr(db_property, key, value)
    
    db_property.updated_at = datetime.now()
    db.add(db_property)
    db.commit()
    db.refresh(db_property)
    return db_property

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

def create_assay_type(db: Session, assay_type: schemas.AssayTypeCreate):
    # Create the assay type
    current_time = datetime.now(timezone.utc)
    db_assay_type = models.AssayType(
        name=assay_type.name,
        description=assay_type.description,
        created_on=current_time,
        updated_on=current_time
    )
    db.add(db_assay_type)
    db.commit()
    db.refresh(db_assay_type)
    
    # Add properties to the assay type if provided
    if assay_type.property_ids:
        properties = db.query(models.Property).filter(models.Property.id.in_(assay_type.property_ids)).all()
        db_assay_type.properties.extend(properties)
        db.commit()
        db.refresh(db_assay_type)
    
    return db_assay_type

def update_assay_type(db: Session, assay_type_id: int, assay_type: schemas.AssayTypeUpdate):
    db_assay_type = db.query(models.AssayType).filter(models.AssayType.id == assay_type_id).first()
    
    # Update basic fields
    update_data = {k: v for k, v in assay_type.model_dump(exclude_unset=True).items() if k != 'property_ids'}
    for key, value in update_data.items():
        setattr(db_assay_type, key, value)
    
    # Update properties if provided
    if assay_type.property_ids is not None:
        # Clear existing properties
        db_assay_type.properties = []
        # Add new properties
        properties = db.query(models.Property).filter(models.Property.id.in_(assay_type.property_ids)).all()
        db_assay_type.properties.extend(properties)
    
    db_assay_type.updated_on = datetime.now()
    db.add(db_assay_type)
    db.commit()
    db.refresh(db_assay_type)
    return db_assay_type

def delete_assay_type(db: Session, assay_type_id: int):
    db_assay_type = db.query(models.AssayType).filter(models.AssayType.id == assay_type_id).first()
    db.delete(db_assay_type)
    db.commit()
    return db_assay_type

# Assay CRUD operations
def get_assay(db: Session, assay_id: int):
    return db.query(models.Assay).filter(models.Assay.id == assay_id).first()

def get_assays(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Assay).offset(skip).limit(limit).all()

def create_assay(db: Session, assay: schemas.AssayCreate):
    # Create the assay
    db_assay = models.Assay(
        name=assay.name,
        description=assay.description,
        assay_type_id=assay.assay_type_id,
        created_at=datetime.now()
    )
    db.add(db_assay)
    db.commit()
    db.refresh(db_assay)
    
    # Add properties to the assay if provided
    if assay.property_ids:
        properties = db.query(models.Property).filter(models.Property.id.in_(assay.property_ids)).all()
        db_assay.properties.extend(properties)
        db.commit()
        db.refresh(db_assay)
    
    return db_assay

def update_assay(db: Session, assay_id: int, assay: schemas.AssayUpdate):
    db_assay = db.query(models.Assay).filter(models.Assay.id == assay_id).first()
    
    # Update basic fields
    update_data = {k: v for k, v in assay.model_dump(exclude_unset=True).items() if k != 'property_ids'}
    for key, value in update_data.items():
        setattr(db_assay, key, value)
    
    # Update properties if provided
    if assay.property_ids is not None:
        # Clear existing properties
        db_assay.properties = []
        # Add new properties
        properties = db.query(models.Property).filter(models.Property.id.in_(assay.property_ids)).all()
        db_assay.properties.extend(properties)
    
    db_assay.updated_at = datetime.now()
    db.add(db_assay)
    db.commit()
    db.refresh(db_assay)
    return db_assay

def delete_assay(db: Session, assay_id: int):
    db_assay = db.query(models.Assay).filter(models.Assay.id == assay_id).first()
    if db_assay is None:
        return None
    
    # Save the ID for return
    assay_id = db_assay.id
    
    db.delete(db_assay)
    db.commit()
    return {"id": assay_id}

def get_compounds_ex(db: Session, query_params: schemas.CompoundQueryParams):
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
                "measurements": {}
            }
        
        # Get property name
        property_name = db.query(models.Property).filter(models.Property.id == result.property_id).first().name
        grouped_results[assay_id]["measurements"][property_name] = result.result_value
    
    return list(grouped_results.values())

def create_assay_result(db: Session, assay_result: schemas.AssayResultCreate):
    """Create a single assay result entry for a specific property"""
    # Validate assay exists
    assay = db.query(models.Assay).filter(models.Assay.id == assay_result.assay_id).first()
    if assay is None:
        raise HTTPException(status_code=404, detail=f"Assay with ID {assay_result.assay_id} not found")
    
    # Validate batch exists
    batch = db.query(models.Batch).filter(models.Batch.id == assay_result.batch_id).first()
    if batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {assay_result.batch_id} not found")
    
    # Validate property exists
    property = db.query(models.Property).filter(models.Property.id == assay_result.property_id).first()
    if property is None:
        raise HTTPException(status_code=404, detail=f"Property with ID {assay_result.property_id} not found")
    
    # Check if the property is associated with the assay
    if property not in assay.properties:
        raise HTTPException(
            status_code=400, 
            detail=f"Property with ID {assay_result.property_id} is not associated with assay '{assay.name}'"
        )
    
    # Create the assay result
    db_assay_result = models.AssayResult(
        assay_id=assay_result.assay_id,
        batch_id=assay_result.batch_id,
        property_id=assay_result.property_id,
        result_value=assay_result.result_value
    )
    
    db.add(db_assay_result)
    db.commit()
    db.refresh(db_assay_result)
    return db_assay_result

def create_batch_assay_results(db: Session, batch_results: schemas.BatchAssayResultsCreate):
    """Create multiple assay results for a batch in a single transaction"""
    # Validate assay exists
    assay = db.query(models.Assay).filter(models.Assay.id == batch_results.assay_id).first()
    if assay is None:
        raise HTTPException(status_code=404, detail=f"Assay with ID {batch_results.assay_id} not found")
    
    # Validate batch exists
    batch = db.query(models.Batch).filter(models.Batch.id == batch_results.batch_id).first()
    if batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {batch_results.batch_id} not found")
    
    # Get all property names and IDs for the assay
    properties = {prop.name: prop.id for prop in assay.properties}
    
    # Validate all property names in measurements
    for prop_name in batch_results.measurements.keys():
        if prop_name not in properties:
            raise HTTPException(
                status_code=400, 
                detail=f"Property '{prop_name}' is not associated with assay '{assay.name}'"
            )
    
    # Create assay results for each property
    created_results = []
    for prop_name, result_value in batch_results.measurements.items():
        property_id = properties[prop_name]
        db_assay_result = models.AssayResult(
            assay_id=batch_results.assay_id,
            batch_id=batch_results.batch_id,
            property_id=property_id,
            result_value=result_value
        )
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
        "measurements": batch_results.measurements
    }

def update_assay_result(db: Session, assay_result_id: int, assay_result: schemas.AssayResultUpdate):
    """Update a single assay result"""
    db_assay_result = db.query(models.AssayResult).filter(models.AssayResult.id == assay_result_id).first()
    if db_assay_result is None:
        raise HTTPException(status_code=404, detail=f"Assay result with ID {assay_result_id} not found")
    
    # Update the result value
    db_assay_result.result_value = assay_result.result_value
    
    db.add(db_assay_result)
    db.commit()
    db.refresh(db_assay_result)
    return db_assay_result

def update_batch_assay_results(db: Session, assay_id: int, batch_id: int, measurements: Dict[str, float]):
    """Update multiple assay results for a batch in a single transaction"""
    # Validate assay exists
    assay = db.query(models.Assay).filter(models.Assay.id == assay_id).first()
    if assay is None:
        raise HTTPException(status_code=404, detail=f"Assay with ID {assay_id} not found")
    
    # Validate batch exists
    batch = db.query(models.Batch).filter(models.Batch.id == batch_id).first()
    if batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {batch_id} not found")
    
    # Get all property names and IDs for the assay
    properties = {prop.name: prop.id for prop in assay.properties}
    
    # Validate all property names in measurements
    for prop_name in measurements.keys():
        if prop_name not in properties:
            raise HTTPException(
                status_code=400, 
                detail=f"Property '{prop_name}' is not associated with assay '{assay.name}'"
            )
    
    # Get existing results for this batch/assay
    existing_results = db.query(models.AssayResult).filter(
        models.AssayResult.assay_id == assay_id,
        models.AssayResult.batch_id == batch_id
    ).all()
    
    # Create a lookup by property_id
    existing_by_property = {result.property_id: result for result in existing_results}
    
    # Update existing or create new results
    updated_results = []
    for prop_name, result_value in measurements.items():
        property_id = properties[prop_name]
        
        if property_id in existing_by_property:
            # Update existing result
            existing_result = existing_by_property[property_id]
            existing_result.result_value = result_value
            db.add(existing_result)
            updated_results.append(existing_result)
        else:
            # Create new result
            db_assay_result = models.AssayResult(
                assay_id=assay_id,
                batch_id=batch_id,
                property_id=property_id,
                result_value=result_value
            )
            db.add(db_assay_result)
            updated_results.append(db_assay_result)
    
    db.commit()
    
    # Return in grouped format
    return {
        "assay_id": assay_id,
        "batch_id": batch_id,
        "assay_name": assay.name,
        "measurements": measurements
    }

def delete_assay_result(db: Session, assay_result_id: int):
    """Delete a single assay result"""
    db_assay_result = db.query(models.AssayResult).filter(models.AssayResult.id == assay_result_id).first()
    if db_assay_result is None:
        raise HTTPException(status_code=404, detail=f"Assay result with ID {assay_result_id} not found")
    
    assay_result_id = db_assay_result.id
    db.delete(db_assay_result)
    db.commit()
    return {"id": assay_result_id}

def delete_batch_assay_results(db: Session, assay_id: int, batch_id: int):
    """Delete all assay results for a batch/assay combination"""
    # Validate assay exists
    assay = db.query(models.Assay).filter(models.Assay.id == assay_id).first()
    if assay is None:
        raise HTTPException(status_code=404, detail=f"Assay with ID {assay_id} not found")
    
    # Validate batch exists
    batch = db.query(models.Batch).filter(models.Batch.id == batch_id).first()
    if batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {batch_id} not found")
    
    # Delete all results for this batch/assay
    results = db.query(models.AssayResult).filter(
        models.AssayResult.assay_id == assay_id,
        models.AssayResult.batch_id == batch_id
    ).all()
    
    if not results:
        raise HTTPException(
            status_code=404, 
            detail=f"No assay results found for batch {batch_id} and assay {assay_id}"
        )
    
    for result in results:
        db.delete(result)
    
    db.commit()
    
    return {
        "assay_id": assay_id,
        "batch_id": batch_id,
        "deleted_count": len(results)
    }

# BatchDetail CRUD operations
def get_batch_detail(db: Session, batch_detail_id: int):
    return db.query(models.BatchDetail).filter(models.BatchDetail.id == batch_detail_id).first()

def get_batch_details_by_batch(db: Session, batch_id: int, skip: int = 0, limit: int = 100):
    return db.query(models.BatchDetail).filter(models.BatchDetail.batch_id == batch_id).offset(skip).limit(limit).all()

def create_batch_detail(db: Session, batch_detail: schemas.BatchDetailCreate):
    db_batch_detail = models.BatchDetail(
        batch_id=batch_detail.batch_id,
        property_id=batch_detail.property_id,
        result_value=batch_detail.result_value
    )
    db.add(db_batch_detail)
    db.commit()
    db.refresh(db_batch_detail)
    return db_batch_detail

def update_batch_detail(db: Session, batch_detail_id: int, batch_detail: schemas.BatchDetailUpdate):
    db_batch_detail = db.query(models.BatchDetail).filter(models.BatchDetail.id == batch_detail_id).first()
    if db_batch_detail:
        update_data = batch_detail.model_dump(exclude_unset=True)
        for key, value in update_data.items():
            setattr(db_batch_detail, key, value)
        db.add(db_batch_detail)
        db.commit()
        db.refresh(db_batch_detail)
    return db_batch_detail

def delete_batch_detail(db: Session, batch_detail_id: int):
    db_batch_detail = db.query(models.BatchDetail).filter(models.BatchDetail.id == batch_detail_id).first()
    if db_batch_detail:
        db.delete(db_batch_detail)
        db.commit()
    return db_batch_detail




