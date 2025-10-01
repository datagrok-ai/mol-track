from fastapi import HTTPException
from sqlalchemy import and_
from sqlalchemy.orm import Session
from sqlalchemy.orm import joinedload
from app.crud.properties import enrich_model
from app import crud, models


def enrich_batch(batch: models.Batch) -> models.BatchResponse:
    batch_resp = enrich_model(batch, models.BatchResponse, "batch_details", "batch_id")
    if batch.compound:
        batch_resp.compound = enrich_model(batch.compound, models.CompoundResponse, "compound_details", "compound_id")
    return batch_resp


def get_batch_by_synonym(db: Session, property_value: str, property_name: str = None, enrich: bool = True):
    if not property_value:
        return None

    filters = [models.Property.semantic_type_id == 1, models.BatchDetail.value_string == property_value]

    if property_name:
        filters.append(models.Property.name == property_name)

    batch = (
        db.query(models.Batch)
        .join(models.Batch.batch_details)
        .join(models.BatchDetail.property)
        .options(joinedload(models.Batch.batch_details).joinedload(models.BatchDetail.property))
        .filter(and_(*filters))
        .first()
    )

    if not batch:
        return None

    return enrich_batch(batch) if enrich else batch


def get_batches(db: Session, skip: int = 0, limit: int = 100):
    batches = db.query(models.Batch).offset(skip).limit(limit).all()
    return [enrich_batch(batch) for batch in batches]


def get_batches_by_compound(db: Session, compound_id: int, skip: int = 0, limit: int = 100):
    return db.query(models.Batch).filter(models.Batch.compound_id == compound_id).offset(skip).limit(limit).all()


def delete_batch_by_synonym(
    db: Session,
    property_value: str,
    property_name: str,
):
    db_batch = get_batch_by_synonym(db, property_value, property_name, False)
    if db_batch is None:
        raise HTTPException(status_code=404, detail="Batch not found")

    assay_results = crud.get_all_assay_results_for_batch(db, db_batch.id)
    if assay_results:
        raise HTTPException(status_code=400, detail="Batch has dependent assay results")

    db.delete(db_batch)
    db.commit()
    return db_batch
