from sqlalchemy.orm import Session
from fastapi import HTTPException
from rdkit import Chem
from sqlalchemy import select, text
from sqlalchemy.orm import joinedload
from sqlalchemy import and_
from app.crud.properties import enrich_properties
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


def enrich_compound(compound: models.Compound) -> models.CompoundResponse:
    return models.CompoundResponse(
        **compound.dict(), properties=enrich_properties(compound, "compound_details", "compound_id")
    )


def read_compounds(db: Session, skip: int = 0, limit: int = 100):
    compounds = db.query(models.Compound).offset(skip).limit(limit).all()
    return [enrich_compound(c) for c in compounds]


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

    return enrich_compound(compound) if enrich else compound


def update_compound(db: Session, property_value: str, property_name: str, update_data: models.CompoundUpdate):
    db_compound = get_compound_by_synonym(db, property_value=property_value, property_name=property_name, enrich=False)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")

    if db_compound.batches and update_data.smiles is not None:
        raise HTTPException(
            status_code=400, detail="Structure updates are not allowed for compounds with batches attached."
        )

    for field in ["original_molfile", "is_archived", "canonical_smiles"]:
        value = getattr(update_data, field)
        if value is not None:
            setattr(db_compound, field, value)

    if update_data.properties:
        for detail in update_data.properties:
            db_property = db.get(models.Property, detail.property_id)
            if db_property is None:
                raise HTTPException(status_code=404, detail=f"Property with id {detail.property_id} not found")

            expected_type = db_property.value_type
            if (
                expected_type not in type_casting_utils.value_type_to_field
                or expected_type not in type_casting_utils.value_type_cast_map
            ):
                raise HTTPException(status_code=400, detail=f"Unsupported value type '{expected_type}'")

            try:
                cast_value = type_casting_utils.value_type_cast_map[expected_type](detail.value)
            except Exception:
                raise HTTPException(
                    status_code=400, detail=f"Invalid value '{detail.value}' for type '{expected_type}'"
                )

            stmt = select(models.CompoundDetail).where(
                models.CompoundDetail.compound_id == db_compound.id,
                models.CompoundDetail.property_id == detail.property_id,
            )
            existing_detail = db.execute(stmt).scalars().first()
            if not existing_detail:
                raise HTTPException(
                    status_code=400, detail=f"No existing CompoundDetail found for property_id {detail.property_id}"
                )

            setattr(existing_detail, type_casting_utils.value_type_to_field[expected_type], cast_value)
            existing_detail.updated_by = admin.admin_user_id

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
