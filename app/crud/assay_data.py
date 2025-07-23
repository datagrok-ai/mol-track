from typing import List
from sqlalchemy.orm import Session
from app import models


# === Assay-related operations ===
def get_assay(db: Session, assay_id: int):
    return db.query(models.Assay).filter(models.Assay.id == assay_id).first()


def get_assays(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Assay).offset(skip).limit(limit).all()


# === AssayRun-related operations ===
def get_assay_run(db: Session, assay_run_id: int):
    # Get the assay
    assay = db.query(models.AssayRun).filter(models.AssayRun.id == assay_run_id).first()

    if assay:
        # Get properties associated with this assay through assay details
        assay_details = db.query(models.AssayRunDetail).filter(models.AssayRunDetail.assay_run_id == assay_run_id).all()
        property_ids = [detail.property_id for detail in assay_details]

        # If no properties from assay details, get them from the assay type
        if not property_ids and assay.assay:
            assay_properties = (
                db.query(models.AssayProperty).filter(models.AssayProperty.assay_id == assay.assay_id).all()
            )
            property_ids = [prop.property_id for prop in assay_properties]

        # Get the property objects
        if property_ids:
            properties = db.query(models.Property).filter(models.Property.id.in_(property_ids)).all()
            # Add properties to assay
            assay.properties = properties
        else:
            assay.properties = []

    return assay


def get_assay_runs(db: Session, skip: int = 0, limit: int = 100):
    # Get assay_runs with pagination
    assay_runs = db.query(models.AssayRun).offset(skip).limit(limit).all()

    # For each assay, add its properties
    for assay_run in assay_runs:
        # Get properties associated with this assay through assay details
        assay_details = db.query(models.AssayRunDetail).filter(models.AssayRunDetail.assay_run_id == assay_run.id).all()
        property_ids = [detail.property_id for detail in assay_details]

        # If no properties from assay details, get them from the assay type
        if not property_ids and assay_run.assay:
            assay_properties = (
                db.query(models.AssayProperty).filter(models.AssayProperty.assay_id == assay_run.assay_id).all()
            )
            property_ids = [prop.property_id for prop in assay_properties]

        # Get the property objects
        if property_ids:
            properties = db.query(models.Property).filter(models.Property.id.in_(property_ids)).all()
            # Add properties to assay
            assay_run.properties = properties
        else:
            assay_run.properties = []

    return assay_runs


# === AssayResult-related operations ===
def get_assay_result(db: Session, assay_result_id: int):
    return db.query(models.AssayResult).filter(models.AssayResult.id == assay_result_id).first()


def get_assay_results(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.AssayResult).offset(skip).limit(limit).all()


def get_all_assay_results_for_batch(db: Session, batch_id: int) -> List[models.AssayResult]:
    return db.query(models.AssayResult).filter(models.AssayResult.batch_id == batch_id).all()


def get_batch_assay_results(db: Session, batch_id: int):
    """Get all assay results for a specific batch"""
    results = get_all_assay_results_for_batch(db, batch_id)

    # Group results by assay_run_id
    grouped_results = {}
    for result in results:
        assay_run_id = result.assay_run_id
        if assay_run_id not in grouped_results:
            # Get the assay name
            assay_run = db.query(models.AssayRun).filter(models.AssayRun.id == assay_run_id).first()
            assay = db.query(models.Assay).filter(models.Assay.id == assay_run.assay_id).first()
            assay_name = assay.name if assay else "Unknown Assay"

            grouped_results[assay_run_id] = {
                "assay_run_id": assay_run_id,
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
            grouped_results[assay_run_id]["measurements"][property_name] = {
                "qualifier": result.value_qualifier,
                "value": value,
            }
        else:
            grouped_results[assay_run_id]["measurements"][property_name] = value

    return list(grouped_results.values())
