from fastapi import FastAPI, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List

# Handle both package imports and direct execution
try:
    # When imported as a package (for tests)
    from . import models, schemas, crud
    from .database import SessionLocal, engine
except ImportError:
    # When run directly
    import models, schemas, crud
    from database import SessionLocal, engine

#models.Base.metadata.create_all(bind=engine)

app = FastAPI(title="MolTrack API", description="API for managing chemical compounds and batches")

# Dependency
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

# Compounds endpoints
@app.post("/compounds/", response_model=schemas.Compound)
def create_compound(compound: schemas.CompoundCreate, db: Session = Depends(get_db)):
    return crud.create_compound(db=db, compound=compound)

@app.post("/compounds/batch/", response_model=List[schemas.Compound])
def create_compounds_batch(batch: schemas.CompoundBatchCreate, db: Session = Depends(get_db)):
    """
    Create multiple compounds from a list of SMILES strings.
    
    All SMILES must be valid and not already exist in the database.
    If any SMILES is invalid or already exists, the entire batch will fail.
    """
    return crud.create_compounds_batch(db=db, smiles_list=batch.compounds)

@app.get("/compounds/", response_model=List[schemas.Compound])
def read_compounds(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    compounds = crud.get_compounds(db, skip=skip, limit=limit)
    return compounds

@app.get("/compounds/{compound_id}", response_model=schemas.Compound)
def read_compound(compound_id: int, db: Session = Depends(get_db)):
    db_compound = crud.get_compound(db, compound_id=compound_id)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")
    return db_compound

@app.put("/compounds/{compound_id}", response_model=schemas.Compound)
def update_compound(compound_id: int, compound: schemas.CompoundUpdate, db: Session = Depends(get_db)):
    db_compound = crud.get_compound(db, compound_id=compound_id)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")
    return crud.update_compound(db=db, compound_id=compound_id, compound=compound)

@app.delete("/compounds/{compound_id}", response_model=schemas.Compound)
def delete_compound(compound_id: int, db: Session = Depends(get_db)):
    db_compound = crud.get_compound(db, compound_id=compound_id)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")
    return crud.delete_compound(db=db, compound_id=compound_id)

# Batches endpoints
@app.post("/batches/", response_model=schemas.Batch)
def create_batch(batch: schemas.BatchCreate, db: Session = Depends(get_db)):
    return crud.create_batch(db=db, batch=batch)

@app.get("/batches/", response_model=List[schemas.Batch])
def read_batches(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    batches = crud.get_batches(db, skip=skip, limit=limit)
    return batches

@app.get("/batches/{batch_id}", response_model=schemas.Batch)
def read_batch(batch_id: int, db: Session = Depends(get_db)):
    db_batch = crud.get_batch(db, batch_id=batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail="Batch not found")
    return db_batch

@app.put("/batches/{batch_id}", response_model=schemas.Batch)
def update_batch(batch_id: int, batch: schemas.BatchUpdate, db: Session = Depends(get_db)):
    db_batch = crud.get_batch(db, batch_id=batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail="Batch not found")
    return crud.update_batch(db=db, batch_id=batch_id, batch=batch)

@app.delete("/batches/{batch_id}", response_model=schemas.Batch)
def delete_batch(batch_id: int, db: Session = Depends(get_db)):
    db_batch = crud.get_batch(db, batch_id=batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail="Batch not found")
    return crud.delete_batch(db=db, batch_id=batch_id)

# Properties endpoints
@app.post("/properties/", response_model=schemas.Property)
def create_property(property: schemas.PropertyCreate, db: Session = Depends(get_db)):
    return crud.create_property(db=db, property=property)

@app.get("/properties/", response_model=List[schemas.Property])
def read_properties(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    properties = crud.get_properties(db, skip=skip, limit=limit)
    return properties

@app.get("/properties/{property_id}", response_model=schemas.Property)
def read_property(property_id: int, db: Session = Depends(get_db)):
    db_property = crud.get_property(db, property_id=property_id)
    if db_property is None:
        raise HTTPException(status_code=404, detail="Property not found")
    return db_property

# Assay endpoints
@app.post("/assays/", response_model=schemas.Assay)
def create_assay(assay: schemas.AssayCreate, db: Session = Depends(get_db)):
    # Validate that all property IDs exist
    if assay.property_ids:
        for property_id in assay.property_ids:
            db_property = crud.get_property(db, property_id=property_id)
            if db_property is None:
                raise HTTPException(status_code=404, detail=f"Property with ID {property_id} not found")
    
    return crud.create_assay(db=db, assay=assay)

@app.get("/assays/", response_model=List[schemas.Assay])
def read_assays(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    assays = crud.get_assays(db, skip=skip, limit=limit)
    return assays

@app.get("/assays/{assay_id}", response_model=schemas.Assay)
def read_assay(assay_id: int, db: Session = Depends(get_db)):
    db_assay = crud.get_assay(db, assay_id=assay_id)
    if db_assay is None:
        raise HTTPException(status_code=404, detail="Assay not found")
    return db_assay

@app.put("/assays/{assay_id}", response_model=schemas.Assay)
def update_assay(assay_id: int, assay: schemas.AssayUpdate, db: Session = Depends(get_db)):
    db_assay = crud.get_assay(db, assay_id=assay_id)
    if db_assay is None:
        raise HTTPException(status_code=404, detail="Assay not found")
    
    # Validate that all property IDs exist if provided
    if assay.property_ids:
        for property_id in assay.property_ids:
            db_property = crud.get_property(db, property_id=property_id)
            if db_property is None:
                raise HTTPException(status_code=404, detail=f"Property with ID {property_id} not found")
    
    return crud.update_assay(db=db, assay_id=assay_id, assay=assay)

@app.delete("/assays/{assay_id}", response_model=schemas.Assay)
def delete_assay(assay_id: int, db: Session = Depends(get_db)):
    db_assay = crud.get_assay(db, assay_id=assay_id)
    if db_assay is None:
        raise HTTPException(status_code=404, detail="Assay not found")
    return crud.delete_assay(db=db, assay_id=assay_id) 