from fastapi import HTTPException
from sqlalchemy.orm import Session
from sqlalchemy.orm import selectinload
from app.crud.properties import enrich_properties
from app import models


def enrich_batch(batch: models.Batch) -> models.BatchResponse:
    return models.BatchResponse(**batch.dict(), properties=enrich_properties(batch, "batch_details", "batch_id"))


def get_batch(db: Session, batch_id: int):
    batch = (
        db.query(models.Batch)
        .options(selectinload(models.Batch.properties).selectinload(models.Property.batch_details))
        .filter(models.Batch.id == batch_id)
        .first()
    )
    if not batch:
        return None
    return enrich_batch(batch)


def get_batches(db: Session, skip: int = 0, limit: int = 100):
    batches = db.query(models.Batch).offset(skip).limit(limit).all()
    return [enrich_batch(batch) for batch in batches]


def get_batches_by_compound(db: Session, compound_id: int, skip: int = 0, limit: int = 100):
    return db.query(models.Batch).filter(models.Batch.compound_id == compound_id).offset(skip).limit(limit).all()


def delete_batch(db: Session, batch_id: int):
    db_batch = db.get(models.Batch, batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail="Batch not found")

    db.query(models.BatchDetail).filter(models.BatchDetail.batch_id == batch_id).delete(synchronize_session=False)
    db.query(models.BatchAddition).filter(models.BatchAddition.batch_id == batch_id).delete(synchronize_session=False)
    db.delete(db_batch)
    db.commit()
    return db_batch
