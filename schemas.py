from pydantic import BaseModel, Field, validator
from typing import Optional, List, Dict, Any, Union
from datetime import datetime, date
from enum import Enum

# Property schemas
class ValueType(str, Enum):
    INT = "int"
    DOUBLE = "double"
    BOOL = "bool"
    DATETIME = "datetime"
    STRING = "string"

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

# BatchDetail schemas
class BatchDetailBase(BaseModel):
    batch_id: int
    property_id: int
    result_value: str

class BatchDetailCreate(BatchDetailBase):
    pass

class BatchDetailUpdate(BaseModel):
    result_value: Optional[str] = None

class BatchDetail(BatchDetailBase):
    id: int
    
    class Config:
        orm_mode = True
        from_attributes = True
        populate_by_name = True

# AssayType schemas
class AssayTypeBase(BaseModel):
    name: str
    description: Optional[str] = None

class AssayTypeCreate(AssayTypeBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay type

class AssayTypeUpdate(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    property_ids: Optional[List[int]] = None  # Optional list of property IDs to update

class AssayType(AssayTypeBase):
    id: int
    created_on: Optional[datetime] = None
    updated_on: Optional[datetime] = None
    properties: List[Property] = []

    class Config:
        orm_mode = True
        from_attributes = True
        populate_by_name = True

# Assay schemas
class AssayBase(BaseModel):
    name: str
    description: Optional[str] = None
    assay_type_id: int

class AssayCreate(AssayBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay

class AssayUpdate(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    assay_type_id: Optional[int] = None
    property_ids: Optional[List[int]] = None  # Optional list of property IDs to update

class Assay(AssayBase):
    id: int
    created_at: datetime
    properties: List[Property] = []
    assay_type: Optional[AssayType] = None

    class Config:
        orm_mode = True
        from_attributes = True
        populate_by_name = True

# AssayResult schemas
class AssayResultBase(BaseModel):
    batch_id: int
    assay_id: int
    property_id: int
    result_value: float

class AssayResultCreate(AssayResultBase):
    pass

# For updating a single result value
class AssayResultUpdate(BaseModel):
    result_value: float

# Schema for submitting multiple measurements at once for a batch/assay combination
class BatchAssayResultsCreate(BaseModel):
    assay_id: int
    batch_id: int
    measurements: Dict[str, float]  # Map of property name to result value

class AssayResult(AssayResultBase):
    id: int

    class Config:
        orm_mode = True
        from_attributes = True
        populate_by_name = True

# Schema for returning grouped results for a batch
class BatchAssayResultsResponse(BaseModel):
    assay_id: int
    batch_id: int
    assay_name: str
    measurements: Dict[str, Union[float, str, bool]]

    class Config:
        orm_mode = True

# Batch schemas
class BatchBase(BaseModel):
    batch_number: str
    amount: Optional[float] = None
    amount_unit: Optional[str] = None
    purity: Optional[float] = None
    notes: Optional[str] = None
    expiry_date: Optional[date] = None
    created_by: Optional[int] = None

class BatchCreate(BatchBase):
    compound_id: int

class BatchUpdate(BaseModel):
    batch_number: Optional[str] = None
    amount: Optional[float] = None
    amount_unit: Optional[str] = None
    purity: Optional[float] = None
    notes: Optional[str] = None
    expiry_date: Optional[date] = None

class Batch(BatchBase):
    id: int
    compound_id: int
    created_at: datetime
    created_by: Optional[int] = None
    batch_details: List[BatchDetail] = []

    class Config:
        orm_mode = True
        from_attributes = True

# Compound schemas
class CompoundBase(BaseModel):
    canonical_smiles: Optional[str] = None
    original_molfile: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None

    @validator('inchi', 'inchikey', always=True)
    def set_inchi(cls, v, values):
        if v is not None:
            return v
        # In a real implementation, these would be calculated
        # based on canonical_smiles
        if 'canonical_smiles' in values and values['canonical_smiles'] is not None:
            # Placeholder logic
            if 'inchi' in values:
                return 'InChI=1S/' + values['canonical_smiles']
            else:
                return 'INCHIKEY' + values['canonical_smiles']
        return v

class CompoundCreate(CompoundBase):
    smiles: str
    is_archived: Optional[bool] = False

class CompoundBatchCreate(BaseModel):
    compounds: List[str]  # List of canonical SMILES strings

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
    is_archived: bool
    batches: List[Batch] = []

    class Config:
        orm_mode = True
        from_attributes = True

# Query parameters for compound search
class CompoundQueryParams(BaseModel):
    substructure: Optional[str] = None
    skip: int = 0
    limit: int = 100

    class Config:
        from_attributes = True 