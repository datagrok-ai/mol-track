from pydantic import BaseModel, Field
from typing import Optional, List
from datetime import datetime, date
from enum import Enum

# Property schemas
class ValueType(str, Enum):
    NUMBER = "NUMBER"
    STRING = "STRING"
    BOOL = "BOOL"

class PropertyClass(str, Enum):
    CALCULATED = "CALCULATED"
    MEASURED = "MEASURED"
    PREDICTED = "PREDICTED"

class PropertyBase(BaseModel):
    name: str
    value_type: ValueType
    property_class: PropertyClass
    unit: Optional[str] = None

class PropertyCreate(PropertyBase):
    pass

class PropertyUpdate(BaseModel):
    name: Optional[str] = None
    value_type: Optional[ValueType] = None
    property_class: Optional[PropertyClass] = None
    unit: Optional[str] = None

class Property(PropertyBase):
    id: int
    created_at: datetime

    class Config:
        orm_mode = True
        from_attributes = True
        populate_by_name = True

# Assay schemas
class AssayBase(BaseModel):
    name: str
    description: Optional[str] = None

class AssayCreate(AssayBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay

class AssayUpdate(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    property_ids: Optional[List[int]] = None  # Optional list of property IDs to update

class Assay(AssayBase):
    id: int
    created_at: datetime
    properties: List[Property] = []

    class Config:
        orm_mode = True
        from_attributes = True
        populate_by_name = True

# Batch schemas
class BatchBase(BaseModel):
    batch_number: str
    amount: Optional[float] = None
    amount_unit: Optional[str] = None
    purity: Optional[float] = None
    vendor: Optional[str] = None
    catalog_id: Optional[str] = None
    acquisition_date: Optional[date] = None
    expiry_date: Optional[date] = None
    storage_location: Optional[str] = None
    notes: Optional[str] = None
    created_by: Optional[int] = None

class BatchCreate(BatchBase):
    compound_id: int

class BatchUpdate(BaseModel):
    batch_number: Optional[str] = None
    amount: Optional[float] = None
    amount_unit: Optional[str] = None
    purity: Optional[float] = None
    vendor: Optional[str] = None
    catalog_id: Optional[str] = None
    acquisition_date: Optional[date] = None
    expiry_date: Optional[date] = None
    storage_location: Optional[str] = None
    notes: Optional[str] = None

class Batch(BatchBase):
    id: int
    compound_id: int
    created_at: datetime

    class Config:
        orm_mode = True
        from_attributes = True

# Compound schemas
class CompoundBase(BaseModel):
    canonical_smiles: str
    original_molfile: Optional[str] = None
    inchi: str
    inchikey: str
    is_archived: Optional[bool] = False

class CompoundCreate(BaseModel):
    smiles: str
    original_molfile: Optional[str] = None
    is_archived: Optional[bool] = False

class CompoundBatchCreate(BaseModel):
    compounds: List[str]  # List of SMILES strings

class CompoundUpdate(BaseModel):
    canonical_smiles: Optional[str] = None
    original_molfile: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    is_archived: Optional[bool] = None

class Compound(CompoundBase):
    id: int
    created_at: datetime
    updated_at: datetime
    batches: List[Batch] = []

    class Config:
        orm_mode = True
        from_attributes = True 