from sqlalchemy.orm import Session
from fastapi import HTTPException
from typing import List, Optional
from app.crud.properties import bulk_create_if_not_exists
from app import models

from app.utils import enums
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from datetime import datetime


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


def get_batch_addition_for_addition(db: Session, addition_id: int) -> Optional[models.BatchAddition]:
    return db.query(models.BatchAddition).filter(models.BatchAddition.addition_id == addition_id).first()


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
    db.delete(db_addition)
    db.commit()
    return db_addition
