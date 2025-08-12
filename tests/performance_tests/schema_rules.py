from app.utils.enums import EntityType, ValueType


ASSAY_RUN_DATE = "Assay Run Date"
EPA_BATCH_ID = "epa batch id"
SMILES = "smiles"
ASSAY_NAME = "name"
SCHEMA_GENERATING_RULES = {
    EntityType.COMPOUND: {
        "identity_columns": {
            "property_names": [SMILES],
            "property_mapping": [None],
            "property_type": [ValueType.string],
        },
        "alias": "compound_details",
    },
    EntityType.BATCH: {
        "identity_columns": {
            "property_names": [SMILES, EPA_BATCH_ID],
            "property_mapping": [None, "batch_details"],
            "property_type": [ValueType.string, ValueType.string],
        },
        "alias": "batch_details",
    },
    EntityType.ASSAY: {
        "identity_columns": {
            "property_names": [ASSAY_NAME],
            "property_mapping": ["assay"],
            "property_type": [ValueType.string],
        },
        "alias": "assay",
    },
    EntityType.ASSAY_RUN: {
        "identity_columns": {
            "property_names": [ASSAY_NAME, ASSAY_RUN_DATE],
            "property_mapping": ["assay", "assay_run_details"],
            "property_type": [ValueType.string, ValueType.datetime],
        },
        "alias": "assay_run_details",
    },
    EntityType.ASSAY_RESULT: {
        "identity_columns": {
            "property_names": [ASSAY_NAME, ASSAY_RUN_DATE, EPA_BATCH_ID],
            "property_mapping": ["assay", "assay_run_details", "batch_details"],
            "property_type": [ValueType.string, ValueType.datetime, ValueType.string],
        },
        "alias": "assay_results",
    },
}
