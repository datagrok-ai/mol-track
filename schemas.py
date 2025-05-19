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
        from_attributes = True
        populate_by_name = True

# BatchDetail schemas
class BatchDetailBase(BaseModel):
    batch_id: int
    property_id: int
    value_qualifier: Optional[int] = 0  # 0 for "=", 1 for "<", 2 for ">"
    value_datetime: Optional[datetime] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

class BatchDetailCreate(BatchDetailBase):
    pass

class BatchDetailUpdate(BaseModel):
    value_qualifier: Optional[int] = None
    value_datetime: Optional[datetime] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

class BatchDetail(BatchDetailBase):
    id: int
    
    class Config:
        from_attributes = True
        populate_by_name = True

# AssayType property requirement schemas
class AssayTypePropertyBase(BaseModel):
    assay_type_id: int
    property_id: int
    required: bool = False

class AssayTypePropertyCreate(AssayTypePropertyBase):
    pass

class AssayTypeProperty(AssayTypePropertyBase):
    class Config:
        from_attributes = True
        populate_by_name = True

# AssayType detail schemas for metadata
class AssayTypeDetailBase(BaseModel):
    assay_type_id: int
    property_id: int
    value_datetime: Optional[datetime] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

class AssayTypeDetailCreate(AssayTypeDetailBase):
    pass

class AssayTypeDetailUpdate(BaseModel):
    value_datetime: Optional[datetime] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

class AssayTypeDetail(AssayTypeDetailBase):
    class Config:
        from_attributes = True
        populate_by_name = True

# AssayDetail schemas for assay-specific metadata
class AssayDetailBase(BaseModel):
    assay_id: int
    property_id: int
    value_datetime: Optional[datetime] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

class AssayDetailCreate(AssayDetailBase):
    pass

class AssayDetailUpdate(BaseModel):
    value_datetime: Optional[datetime] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

class AssayDetail(AssayDetailBase):
    class Config:
        from_attributes = True
        populate_by_name = True

# Update AssayType schemas to include details and requirements
class AssayTypeBase(BaseModel):
    name: str
    description: Optional[str] = None

class AssayTypeCreate(AssayTypeBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay type
    property_requirements: List[Dict[str, Any]] = []  # List of property requirements
    property_details: List[Dict[str, Any]] = []  # List of property metadata

class AssayTypeUpdate(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    property_ids: Optional[List[int]] = None  # Optional list of property IDs to update
    property_requirements: Optional[List[Dict[str, Any]]] = None
    property_details: Optional[List[Dict[str, Any]]] = None

class AssayType(AssayTypeBase):
    id: int
    created_on: Optional[datetime] = None
    updated_on: Optional[datetime] = None
    properties: List[Property] = []
    assay_type_details: List[AssayTypeDetail] = []
    property_requirements: List[AssayTypeProperty] = []

    class Config:
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
    assay_type: Optional[AssayType] = None
    assay_details: List[AssayDetail] = []
    properties: List[Property] = []

    class Config:
        from_attributes = True
        populate_by_name = True

# AssayResult schemas
class AssayResultBase(BaseModel):
    batch_id: int
    assay_id: int
    property_id: int
    value_qualifier: Optional[int] = 0  # 0 for "=", 1 for "<", 2 for ">"
    value_num: Optional[float] = None
    value_string: Optional[str] = None
    value_bool: Optional[bool] = None

class AssayResultCreate(AssayResultBase):
    pass

# For updating a single result value
class AssayResultUpdate(BaseModel):
    value_qualifier: Optional[int] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None
    value_bool: Optional[bool] = None

# Schema for submitting multiple measurements at once for a batch/assay combination
class BatchAssayResultsCreate(BaseModel):
    assay_id: int
    batch_id: int
    measurements: Dict[str, Union[float, str, bool, Dict[str, Any]]]  # Map of property name to result value or object

class AssayResult(AssayResultBase):
    id: int

    class Config:
        from_attributes = True
        populate_by_name = True

# Extended response model with backward compatibility field
class AssayResultResponse(AssayResult):
    # Add a computed field for backward compatibility
    result_value: Optional[Union[float, str, bool]] = None

    @validator('result_value', always=True)
    def compute_result_value(cls, v, values):
        """Compute result_value from the appropriate typed value field"""
        if 'value_num' in values and values['value_num'] is not None:
            return values['value_num']
        elif 'value_string' in values and values['value_string'] is not None:
            return values['value_string']
        elif 'value_bool' in values and values['value_bool'] is not None:
            return values['value_bool']
        return None

# Schema for returning grouped results for a batch
class BatchAssayResultsResponse(BaseModel):
    assay_id: int
    batch_id: int
    assay_name: str
    measurements: Dict[str, Union[float, str, bool, Dict[str, Any]]]

    class Config:
        from_attributes = True

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
        from_attributes = True

# Query parameters for compound search
class CompoundQueryParams(BaseModel):
    substructure: Optional[str] = None
    skip: int = 0
    limit: int = 100

    class Config:
        from_attributes = True

# SynonymType schemas
class SynonymTypeBase(BaseModel):
    synonym_level: str  # 'batch' or 'compound'
    name: str
    pattern: str
    description: str
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None

    @validator('synonym_level')
    def validate_synonym_level(cls, value):
        allowed_values = ['BATCH', 'COMPOUND']
        if value.upper() not in allowed_values:
            raise ValueError(f"synonym_level must be one of {allowed_values}")
        return value

class SynonymType(SynonymTypeBase):
    id: int

    class Config:
        from_attributes = True

class SynonymTypeCreate(SynonymTypeBase):
    pass

# Compound synonym schemas
class CompoundSynonymBase(BaseModel):
    compound_id: int  # ID of the compound
    synonym_type_id: int  # ID of the synonym type
    synonym_value: str
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None

class CompoundSynonym(CompoundSynonymBase):
    id: int

    class Config:
        from_attributes = True

class CompoundSynonymCreate(CompoundSynonymBase):
    pass

class CompoundSynonymUpdate(BaseModel):
    compound_id: Optional[int] = None
    synonym_type_id: Optional[int] = None
    synonym_value: Optional[str] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None

# Batch synonym schemas
class BatchSynonymBase(BaseModel):
    batch_id: int  # ID of the batch
    synonym_type_id: int  # ID of the synonym type
    synonym_value: str
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None

class BatchSynonym(BatchSynonymBase):
    id: int

    class Config:
        from_attributes = True

class BatchSynonymCreate(BatchSynonymBase):
    pass

class BatchSynonymUpdate(BaseModel):
    batch_id: Optional[int] = None
    synonym_type_id: Optional[int] = None
    synonym_value: Optional[str] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
