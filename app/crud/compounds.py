from sqlalchemy.orm import Session
from fastapi import HTTPException
from rdkit import Chem
from sqlalchemy import text
from sqlalchemy.orm import joinedload
from sqlalchemy import and_
from app.crud.properties import enrich_model
from app import crud, models
from app.utils.admin_utils import admin
from app.utils import type_casting_utils


def get_compound_by_hash(db: Session, hash_mol: str):
    """
    Search for compounds in the database by hash_mol.
    """
    if not isinstance(hash_mol, str) or len(hash_mol) != 40:
        raise HTTPException(status_code=400, detail=f"Invalid hash_mol format: {hash_mol}")

    return db.query(models.Compound).filter(models.Compound.hash_mol == hash_mol).all()


def read_compounds(db: Session, skip: int = 0, limit: int = 100):
    compounds = db.query(models.Compound).offset(skip).limit(limit).all()
    return [enrich_model(c, models.CompoundResponse, "compound_details", "compound_id") for c in compounds]


def get_compound_by_synonym(db: Session, property_value: str, property_name: str = None, enrich: bool = True):
    if not property_value:
        return None

    filters = [models.Property.semantic_type_id == 1, models.CompoundDetail.value_string == property_value]

    if property_name:
        filters.append(models.Property.name == property_name)

    compound = (
        db.query(models.Compound)
        .join(models.Compound.compound_details)
        .join(models.CompoundDetail.property)
        .options(joinedload(models.Compound.compound_details).joinedload(models.CompoundDetail.property))
        .filter(and_(*filters))
        .first()
    )

    if not compound:
        return None

    return enrich_model(compound, models.CompoundResponse, "compound_details", "compound_id") if enrich else compound


def update_compound(db: Session, property_value: str, property_name: str, update_data: models.CompoundUpdate):
    db_compound = get_compound_by_synonym(db, property_value=property_value, property_name=property_name, enrich=False)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")

    if db_compound.batches and update_data.canonical_smiles is not None:
        raise HTTPException(
            status_code=400,
            detail="Structure updates are not allowed for compounds with batches attached.",
        )

    for field in ("original_molfile", "is_archived", "canonical_smiles"):
        if (value := getattr(update_data, field)) is not None:
            setattr(db_compound, field, value)

    for detail in update_data.properties or []:
        db_property = db.get(models.Property, detail.property_id)
        if db_property is None:
            raise HTTPException(status_code=404, detail=f"Property with id {detail.property_id} not found")

        expected_type = db_property.value_type
        field_name = type_casting_utils.value_type_to_field.get(expected_type)
        cast_func = type_casting_utils.value_type_cast_map.get(expected_type)

        if not field_name or not cast_func:
            raise HTTPException(status_code=400, detail=f"Unsupported value type '{expected_type}'")

        existing_details = (
            db.query(models.CompoundDetail)
            .filter(
                models.CompoundDetail.compound_id == db_compound.id,
                models.CompoundDetail.property_id == detail.property_id,
            )
            .all()
        )

        new_values = detail.value if isinstance(detail.value, list) else [detail.value]
        old_values = [getattr(d, field_name) for d in existing_details]

        to_add = set(new_values) - set(old_values)
        to_remove = set(old_values) - set(new_values)

        if not isinstance(detail.value, list):
            existing_detail = existing_details[0] if existing_details else None
            if not existing_detail:
                raise HTTPException(
                    status_code=400,
                    detail=f"No existing CompoundDetail found for property_id {detail.property_id}",
                )
            try:
                setattr(existing_detail, field_name, cast_func(detail.value))
            except Exception:
                raise HTTPException(
                    status_code=400,
                    detail=f"Invalid value '{detail.value}' for type '{expected_type}'",
                )
            existing_detail.updated_by = admin.admin_user_id
            continue

        for val in to_add:
            try:
                cast_val = cast_func(val)
            except Exception:
                raise HTTPException(
                    status_code=400,
                    detail=f"Invalid value '{val}' for type '{expected_type}'",
                )
            db.add(
                models.CompoundDetail(
                    compound_id=db_compound.id,
                    property_id=detail.property_id,
                    **{field_name: cast_val},
                    created_by=admin.admin_user_id,
                    updated_by=admin.admin_user_id,
                )
            )

        for d in existing_details:
            if getattr(d, field_name) in to_remove:
                db.delete(d)

    db.add(db_compound)
    db.commit()
    db.refresh(db_compound)
    return db_compound


def delete_compound(db: Session, property_value: str, property_name: str):
    db_compound = get_compound_by_synonym(db, property_value=property_value, property_name=property_name, enrich=False)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")

    batches = crud.get_batches_by_compound(db, compound_id=db_compound.id)
    if batches:
        raise HTTPException(status_code=400, detail="Compound has dependent batches")

    db.delete(db_compound)
    db.commit()
    return db_compound


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
        return read_compounds(db, skip=query_params.skip, limit=query_params.limit)
