from sqlalchemy.orm import Session
from fastapi import HTTPException
from rdkit import Chem
from typing import List

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

def get_compounds(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Compound).offset(skip).limit(limit).all()

def create_compound(db: Session, compound: schemas.CompoundCreate):
    # Create RDKit molecule from SMILES
    mol = Chem.MolFromSmiles(compound.smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    # Calculate canonical SMILES, InChI, and InChIKey
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.InchiToInchiKey(inchi)
    
    # Check if compound with this InChIKey already exists
    existing_compound = get_compound_by_inchi_key(db, inchikey)
    if existing_compound:
        raise HTTPException(status_code=400, detail=f"Compound with InChIKey {inchikey} already exists")
    
    # Create the compound with calculated values
    db_compound = models.Compound(
        canonical_smiles=canonical_smiles,
        original_molfile=compound.original_molfile,
        inchi=inchi,
        inchikey=inchikey,
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

def get_compound_batches(db: Session, compound_id: int):
    return db.query(models.Batch).filter(models.Batch.compound_id == compound_id).all()

def create_batch(db: Session, batch: schemas.BatchCreate):
    db_batch = models.Batch(
        compound_id=batch.compound_id,
        batch_number=batch.batch_number,
        amount=batch.amount,
        amount_unit=batch.amount_unit,
        purity=batch.purity,
        vendor=batch.vendor,
        catalog_id=batch.catalog_id,
        acquisition_date=batch.acquisition_date,
        expiry_date=batch.expiry_date,
        storage_location=batch.storage_location,
        notes=batch.notes,
        created_by=batch.created_by
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
    
    db.add(db_property)
    db.commit()
    db.refresh(db_property)
    return db_property

def delete_property(db: Session, property_id: int):
    db_property = db.query(models.Property).filter(models.Property.id == property_id).first()
    db.delete(db_property)
    db.commit()
    return db_property

# Assay CRUD operations
def get_assay(db: Session, assay_id: int):
    return db.query(models.Assay).filter(models.Assay.id == assay_id).first()

def get_assays(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Assay).offset(skip).limit(limit).all()

def create_assay(db: Session, assay: schemas.AssayCreate):
    # Create the assay
    db_assay = models.Assay(
        name=assay.name,
        description=assay.description
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
    
    db.add(db_assay)
    db.commit()
    db.refresh(db_assay)
    return db_assay

def delete_assay(db: Session, assay_id: int):
    db_assay = db.query(models.Assay).filter(models.Assay.id == assay_id).first()
    db.delete(db_assay)
    db.commit()
    return db_assay 