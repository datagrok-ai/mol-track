import uuid
from datetime import date, datetime
from enum import Enum
from typing import Any, Dict, List, Literal, Optional, Union

from pydantic import BaseModel, Field, validator
from rdkit import Chem

from chemistry_utils import (
    generate_hash_layers,
    generate_uuid_hash_mol,
    standardize_mol,
)


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


class ScopeClass(str, Enum):
    BATCH = "BATCH"
    COMPOUND = "COMPOUND"
    ASSAY = "ASSAY"
    SYSTEM = "SYSTEM"


class PropertyBase(BaseModel):
    name: str
    value_type: ValueType
    property_class: PropertyClass
    unit: Optional[str] = None
    semantic_type_id: Optional[int]
    scope: ScopeClass


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
    property_ids: List[
        int
    ] = []  # List of property IDs to associate with this assay type
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
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None
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
    measurements: Dict[
        str, Union[float, str, bool, Dict[str, Any]]
    ]  # Map of property name to result value or object


class AssayResult(AssayResultBase):
    id: int

    class Config:
        from_attributes = True
        populate_by_name = True


# Extended response model with backward compatibility field
class AssayResultResponse(AssayResult):
    # Add a computed field for backward compatibility
    result_value: Optional[Union[float, str, bool]] = None

    @validator("result_value", always=True)
    def compute_result_value(cls, v, values):
        """Compute result_value from the appropriate typed value field"""
        if "value_num" in values and values["value_num"] is not None:
            return values["value_num"]
        elif "value_string" in values and values["value_string"] is not None:
            return values["value_string"]
        elif "value_bool" in values and values["value_bool"] is not None:
            return values["value_bool"]
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
    notes: Optional[str] = None


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
    batch_details: List[BatchDetail] = []

    class Config:
        from_attributes = True


# Compound schemas
class CompoundBase(BaseModel):
    canonical_smiles: Optional[str] = None
    original_molfile: Optional[str] = None
    inchi: Optional[str] = None
    hash_mol: Optional[str] = None
    formula: Optional[str] = None
    hash_tautomer: Optional[uuid.UUID] = None
    hash_canonical_smiles: Optional[uuid.UUID] = None
    hash_no_stereo_smiles: Optional[uuid.UUID] = None
    hash_no_stereo_tautomer: Optional[uuid.UUID] = None
    sgroup_data: Optional[str] = None
    inchikey: Optional[str] = None
    is_archived: Optional[bool] = None

    @validator("inchi", "inchikey", always=True)
    def set_inchi(cls, v, values):
        if v is not None:
            return v
        # In a real implementation, these would be calculated
        # based on canonical_smiles
        if "canonical_smiles" in values and values["canonical_smiles"] is not None:
            # Placeholder logic
            if "inchi" in values:
                return "InChI=1S/" + values["canonical_smiles"]
            else:
                return "INCHIKEY" + values["canonical_smiles"]
        return v


class CompoundCreate(CompoundBase):
    smiles: str
    is_archived: Optional[bool] = None


class CompoundBatchCreate(BaseModel):
    compounds: List[str]  # List of canonical SMILES strings


class CompoundUpdate(BaseModel):
    canonical_smiles: Optional[str] = None
    original_molfile: Optional[str] = None
    inchi: Optional[str] = None
    hash_mol: Optional[str] = None
    formula: Optional[str] = None
    hash_tautomer: Optional[uuid.UUID] = None
    hash_canonical_smiles: Optional[uuid.UUID] = None
    hash_no_stereo_smiles: Optional[uuid.UUID] = None
    hash_no_stereo_tautomer: Optional[uuid.UUID] = None
    sgroup_data: Optional[str] = None
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


class SubstructureQuery(BaseModel):
    pattern: str  # The SMILES pattern to search for
    limit: Optional[int] = 10  # Pagination limit (optional, defaults to 10)
    skip: Optional[int] = 0  # Pagination skip (optional, defaults to 0)


class ExactQuery(BaseModel):
    fields: List[str]


class CompoundQueryParams(BaseModel):
    search_type: str  # Either "substructure" or "exact"
    substructure: Optional[SubstructureQuery] = None
    exact: Optional[ExactQuery] = None


class SearchMethod(BaseModel):
    substructure: Optional[str] = "substructure"
    exact: Optional[str] = "exact"
    similarity: Optional[str] = "similarity"


class CompoundSearchMethod(BaseModel):
    search_method: SearchMethod
    query_smiles: str


class SubstructureSearchParameters(BaseModel):
    skip: int = 0
    limit: int = 100


class ExactSearchParameters(BaseModel):
    fields: List[str]


class SemanticTypeCreate(BaseModel):
    name: str
    description: Optional[str] = None


class SemanticType(SemanticTypeCreate):
    id: int

    class Config:
        orm_mode = True


class ExactSearchParameters(BaseModel):
    field: str
    value: Optional[str] = None


# ---- Search Schemas ----


class CompoundSearchRequest(BaseModel):
    search_method: str  # or use an Enum for validation
    query_smiles: str
    search_parameters: Optional[Dict[str, Any]] = None


class BaseSearchModel(BaseModel):
    search_method: Optional[str] = None  # e.g., "exact", "substructure", "similarity"
    filters: Optional[Dict[str, Any]] = (
        None  # Additional filters (e.g., project, source, date range)
    )
    skip: int = 0  # Pagination: number of records to skip
    limit: int = 100  # Pagination: maximum number of records to return


class SearchCompoundsModel(BaseSearchModel):
    query_smiles: Optional[str] = None  # SMILES string for structure-based searches
    standardization_steps: Optional[List[str]] = None  # Optional standardization steps


class ExactSearchModel(BaseModel):
    query_smiles: str  # SMILES string for the molecule
    standardization_steps: Optional[List[str]] = None  # Optional standardization steps
    hash_mol: Optional[uuid.UUID] = (
        None  # UUID hash generated from the standardized SMILES
    )

    @validator("hash_mol", always=True, pre=True)
    def validate_or_generate_hash(cls, v, values):
        """
        Validate or generate a UUID hash from the standardized SMILES.
        """
        query_smiles = values.get("query_smiles")
        if not query_smiles:
            raise ValueError("query_smiles is required to generate a hash.")

        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(query_smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {query_smiles}")

        # Standardize the molecule
        standardized_mol = standardize_mol(mol)

        # Generate molecular layers
        layers = generate_hash_layers(standardized_mol)

        # Generate the hash if not provided - this is a placeholder
        # this would be GetMolHash
        if v is None:
            return generate_uuid_hash_mol(layers)
        return v


class SearchCompoundStructure(BaseModel):
    search_type: Literal[
        "substructure", "tautomer", "stereo", "similarity"
    ]  # Type of structure search
    query_smiles: str  # SMILES string for the structure search
    search_parameters: Optional[Dict[str, Any]] = Field(
        default_factory=dict,
        description="Additional parameters for the search (e.g., similarity threshold, tautomer options)",
    )

    @validator("query_smiles")
    def validate_smiles(cls, v):
        """
        Validate that the provided SMILES string is valid.
        """
        mol = Chem.MolFromSmiles(v)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {v}")
        return v


class QueryCondition(BaseModel):
    table: Literal["batch", "compounds", "assays"]  # Specify the tables to query
    field: str  # Field/column to filter on
    operator: Literal[
        "=", "!=", ">", "<", ">=", "<=", "LIKE", "IN"
    ]  # expand for supported by rdkit cartridge like @>?
    value: Optional[Any] = None  #
    query_smiles: Optional[str] = None
    columns: Optional[List[str]] = None  # List of columns to return for table


class ComplexQueryRequest(BaseModel):
    conditions: List[QueryCondition]
    logic: Literal["AND", "OR"] = "AND"
