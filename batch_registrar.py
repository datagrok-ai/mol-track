from typing import Dict, Optional
from fastapi import HTTPException
from pytest import Session
from sqlalchemy import select
from compound_registrar import CompoundRegistrar
import crud
import models


class BatchRegistrar(CompoundRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.additions_map = self.load_additions()

    def load_additions(self) -> Dict[str, int]:
        result = self.db.execute(select(models.Addition)).all()
        return {self.normalize_key(a.name): a.id for (a,) in result}

    def process_row(self, row):
        grouped = self.group_data(row)
        compound = self.create_compound(grouped.get("compounds", {}))
        self.create_synonyms(compound.id, grouped.get("compounds_synonyms", {}))
        self.create_properties(compound.id, grouped.get("properties", {}))
        batch = self.create_batch(compound.id)
        self.create_batch_synonyms(batch.id, grouped.get("batches_synonyms", {}))
        self.create_batch_additions(batch.id, grouped.get("batches_additions", {}))

    def create_batch(self, compound_id: int):
        batch_data = models.BatchBase(compound_id=compound_id)
        return crud.create_batch(self.db, batch=batch_data)

    def create_batch_synonyms(self, batch_id: int, synonyms: Dict[str, str]):
        for synonym_name, value in synonyms.items():
            norm_synonym_name = self.normalize_key(synonym_name)
            type_id = self.synonym_type_map.get(norm_synonym_name)
            if type_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown synonym type: {synonym_name}")

            synonym_input = models.BatchSynonymBase(batch_id=batch_id, synonym_type_id=type_id, synonym_value=value)
            crud.create_batch_synonym(self.db, batch_synonym=synonym_input)

    def create_batch_additions(self, batch_id: int, batch_additions: Dict[str, str]):
        for addition_name, value in batch_additions.items():
            norm_addition_name = self.normalize_key(addition_name)
            addition_id = self.additions_map.get(norm_addition_name)
            if addition_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown addition: {addition_name}")
            addition_input = models.BatchAdditionBase(
                batch_id=batch_id, addition_id=addition_id, addition_equivalent=value
            )
            crud.create_batch_addition(self.db, batch_addition=addition_input)
