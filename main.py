from fastapi import FastAPI, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List, Dict, Any

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
    print(db.bind.url);
    print("Database connection successful");
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
def read_compounds(
    query: schemas.CompoundQueryParams = Depends(),
    db: Session = Depends(get_db)
):
    """
    Get a list of compounds with optional filtering by substructure.
    
    - **substructure**: Optional SMILES pattern to search for substructures
    - **skip**: Number of records to skip (for pagination)
    - **limit**: Maximum number of records to return (for pagination)
    """
    compounds = crud.get_compounds_ex(db, query_params=query)
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

# AssayType endpoints
@app.post("/assay-types/", response_model=schemas.AssayType)
def create_assay_type(assay_type: schemas.AssayTypeCreate, db: Session = Depends(get_db)):
    # Validate that all property IDs exist
    if assay_type.property_ids:
        for property_id in assay_type.property_ids:
            db_property = crud.get_property(db, property_id=property_id)
            if db_property is None:
                raise HTTPException(status_code=404, detail=f"Property with ID {property_id} not found")
    
    return crud.create_assay_type(db=db, assay_type=assay_type)

@app.get("/assay-types/", response_model=List[schemas.AssayType])
def read_assay_types(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    assay_types = crud.get_assay_types(db, skip=skip, limit=limit)
    return assay_types

@app.get("/assay-types/{assay_type_id}", response_model=schemas.AssayType)
def read_assay_type(assay_type_id: int, db: Session = Depends(get_db)):
    db_assay_type = crud.get_assay_type(db, assay_type_id=assay_type_id)
    if db_assay_type is None:
        raise HTTPException(status_code=404, detail="Assay type not found")
    return db_assay_type

@app.put("/assay-types/{assay_type_id}", response_model=schemas.AssayType)
def update_assay_type(assay_type_id: int, assay_type: schemas.AssayTypeUpdate, db: Session = Depends(get_db)):
    db_assay_type = crud.get_assay_type(db, assay_type_id=assay_type_id)
    if db_assay_type is None:
        raise HTTPException(status_code=404, detail="Assay type not found")
    
    # Validate that all property IDs exist if provided
    if assay_type.property_ids:
        for property_id in assay_type.property_ids:
            db_property = crud.get_property(db, property_id=property_id)
            if db_property is None:
                raise HTTPException(status_code=404, detail=f"Property with ID {property_id} not found")
    
    return crud.update_assay_type(db=db, assay_type_id=assay_type_id, assay_type=assay_type)

@app.delete("/assay-types/{assay_type_id}", response_model=schemas.AssayType)
def delete_assay_type(assay_type_id: int, db: Session = Depends(get_db)):
    db_assay_type = crud.get_assay_type(db, assay_type_id=assay_type_id)
    if db_assay_type is None:
        raise HTTPException(status_code=404, detail="Assay type not found")
    return crud.delete_assay_type(db=db, assay_type_id=assay_type_id)

# Assay endpoints
@app.post("/assays/", response_model=schemas.Assay)
def create_assay(assay: schemas.AssayCreate, db: Session = Depends(get_db)):
    # Validate that the assay type exists
    db_assay_type = crud.get_assay_type(db, assay_type_id=assay.assay_type_id)
    if db_assay_type is None:
        raise HTTPException(status_code=404, detail=f"Assay type with ID {assay.assay_type_id} not found")
    
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

@app.delete("/assays/{assay_id}")
def delete_assay(assay_id: int, db: Session = Depends(get_db)):
    db_assay = crud.get_assay(db, assay_id=assay_id)
    if db_assay is None:
        raise HTTPException(status_code=404, detail="Assay not found")
    result = crud.delete_assay(db=db, assay_id=assay_id)
    return result

# AssayResult endpoints
@app.post("/assay-results/", response_model=schemas.AssayResult)
def create_assay_result(assay_result: schemas.AssayResultCreate, db: Session = Depends(get_db)):
    """
    Create a single assay result entry for a specific property.
    
    - **assay_id**: The ID of the assay
    - **batch_id**: The ID of the batch
    - **property_id**: The ID of the property
    - **result_value**: The value of the measurement
    """
    return crud.create_assay_result(db=db, assay_result=assay_result)

@app.post("/batch-assay-results/", response_model=schemas.BatchAssayResultsResponse)
def create_batch_assay_results(batch_results: schemas.BatchAssayResultsCreate, db: Session = Depends(get_db)):
    """
    Register multiple measurements for a batch against an assay at once.
    
    - **assay_id**: The ID of the assay
    - **batch_id**: The ID of the batch
    - **measurements**: A dictionary mapping property names to their values
    """
    return crud.create_batch_assay_results(db=db, batch_results=batch_results)

@app.get("/assay-results/", response_model=List[schemas.AssayResult])
def read_assay_results(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """
    Get a list of all individual assay results.
    """
    return crud.get_assay_results(db, skip=skip, limit=limit)

@app.get("/assay-results/{assay_result_id}", response_model=schemas.AssayResult)
def read_assay_result(assay_result_id: int, db: Session = Depends(get_db)):
    """
    Get a specific assay result by ID.
    """
    db_assay_result = crud.get_assay_result(db, assay_result_id=assay_result_id)
    if db_assay_result is None:
        raise HTTPException(status_code=404, detail="Assay result not found")
    return db_assay_result

@app.get("/batches/{batch_id}/assay-results", response_model=List[schemas.BatchAssayResultsResponse])
def read_batch_assay_results(batch_id: int, db: Session = Depends(get_db)):
    """
    Get all assay results for a specific batch, grouped by assay.
    """
    return crud.get_batch_assay_results(db, batch_id=batch_id)

@app.put("/assay-results/{assay_result_id}", response_model=schemas.AssayResult)
def update_assay_result(assay_result_id: int, assay_result: schemas.AssayResultUpdate, db: Session = Depends(get_db)):
    """
    Update a single assay result entry.
    """
    return crud.update_assay_result(db=db, assay_result_id=assay_result_id, assay_result=assay_result)

@app.put("/assays/{assay_id}/batches/{batch_id}/results", response_model=schemas.BatchAssayResultsResponse)
def update_batch_assay_results(
    assay_id: int, 
    batch_id: int, 
    measurements: Dict[str, float], 
    db: Session = Depends(get_db)
):
    """
    Update multiple measurements for a batch and assay at once.
    
    The request body should be a dictionary mapping property names to their values.
    """
    return crud.update_batch_assay_results(
        db=db, 
        assay_id=assay_id, 
        batch_id=batch_id, 
        measurements=measurements
    )

@app.delete("/assay-results/{assay_result_id}")
def delete_assay_result(assay_result_id: int, db: Session = Depends(get_db)):
    """
    Delete a single assay result.
    """
    db_assay_result = crud.get_assay_result(db, assay_result_id=assay_result_id)
    if db_assay_result is None:
        raise HTTPException(status_code=404, detail="Assay result not found")
    return crud.delete_assay_result(db=db, assay_result_id=assay_result_id)

@app.delete("/assays/{assay_id}/batches/{batch_id}/results")
def delete_batch_assay_results(assay_id: int, batch_id: int, db: Session = Depends(get_db)):
    """
    Delete all assay results for a batch/assay combination.
    
    - **assay_id**: The ID of the assay
    - **batch_id**: The ID of the batch
    """
    result = crud.delete_batch_assay_results(db=db, assay_id=assay_id, batch_id=batch_id)
    return result

# BatchDetail endpoints
@app.post("/batch-details/", response_model=schemas.BatchDetail)
def create_batch_detail(batch_detail: schemas.BatchDetailCreate, db: Session = Depends(get_db)):
    """
    Create a batch detail entry.
    
    - **batch_id**: The ID of the batch
    - **property_id**: The ID of the property
    - **result_value**: The value for this property
    """
    # Validate batch exists
    db_batch = crud.get_batch(db, batch_id=batch_detail.batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {batch_detail.batch_id} not found")
    
    # Validate property exists
    db_property = crud.get_property(db, property_id=batch_detail.property_id)
    if db_property is None:
        raise HTTPException(status_code=404, detail=f"Property with ID {batch_detail.property_id} not found")
    
    return crud.create_batch_detail(db=db, batch_detail=batch_detail)

@app.get("/batch-details/{batch_detail_id}", response_model=schemas.BatchDetail)
def read_batch_detail(batch_detail_id: int, db: Session = Depends(get_db)):
    """
    Get a specific batch detail by ID.
    
    - **batch_detail_id**: The ID of the batch detail
    """
    db_batch_detail = crud.get_batch_detail(db, batch_detail_id=batch_detail_id)
    if db_batch_detail is None:
        raise HTTPException(status_code=404, detail="Batch detail not found")
    return db_batch_detail

@app.get("/batches/{batch_id}/details", response_model=List[schemas.BatchDetail])
def read_batch_details(batch_id: int, skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """
    Get all details for a specific batch.
    
    - **batch_id**: The ID of the batch
    - **skip**: Number of records to skip
    - **limit**: Maximum number of records to return
    """
    # Validate batch exists
    db_batch = crud.get_batch(db, batch_id=batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {batch_id} not found")
    
    return crud.get_batch_details_by_batch(db, batch_id=batch_id, skip=skip, limit=limit)

@app.put("/batch-details/{batch_detail_id}", response_model=schemas.BatchDetail)
def update_batch_detail(batch_detail_id: int, batch_detail: schemas.BatchDetailUpdate, db: Session = Depends(get_db)):
    """
    Update a batch detail.
    
    - **batch_detail_id**: The ID of the batch detail
    - **result_value**: The new value
    """
    db_batch_detail = crud.get_batch_detail(db, batch_detail_id=batch_detail_id)
    if db_batch_detail is None:
        raise HTTPException(status_code=404, detail="Batch detail not found")
    
    updated = crud.update_batch_detail(db=db, batch_detail_id=batch_detail_id, batch_detail=batch_detail)
    if updated is None:
        raise HTTPException(status_code=500, detail="Failed to update batch detail")
    
    return updated

@app.delete("/batch-details/{batch_detail_id}", response_model=Dict[str, Any])
def delete_batch_detail(batch_detail_id: int, db: Session = Depends(get_db)):
    """
    Delete a batch detail.
    
    - **batch_detail_id**: The ID of the batch detail
    """
    db_batch_detail = crud.get_batch_detail(db, batch_detail_id=batch_detail_id)
    if db_batch_detail is None:
        raise HTTPException(status_code=404, detail="Batch detail not found")
    
    deleted = crud.delete_batch_detail(db=db, batch_detail_id=batch_detail_id)
    if deleted is None:
        raise HTTPException(status_code=500, detail="Failed to delete batch detail")
    
    return {"id": batch_detail_id, "deleted": True} 