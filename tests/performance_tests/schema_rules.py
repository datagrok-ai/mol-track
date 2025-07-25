ASSAY_RUN_DATE = "Assay Run Date"
EPA_BATCH_ID = "epa batch id"
SMILES = "smiles"
ASSAY_NAME = "name"
SCHEMA_GENERATING_RULES = {
    "compound": {
        "num_properties": 5,
        "identity_columns": {"property_names": [SMILES], "property_mapping": [""], "property_type": ["string"]},
        "alias": "compound_details",
    },
    "batch": {
        "num_properties": 10,
        "identity_columns": {
            "property_names": [SMILES, EPA_BATCH_ID],
            "property_mapping": ["", "batch_details"],
            "property_type": ["string", "string"],
        },
        "alias": "batch_details",
    },
    "assay": {
        "num_properties": 2,
        "identity_columns": {
            "property_names": [ASSAY_NAME],
            "property_mapping": ["assay"],
            "property_type": ["string"],
        },
        "alias": "assay",
    },
    "assay_run": {
        "num_properties": 10,
        "identity_columns": {
            "property_names": [ASSAY_NAME, ASSAY_RUN_DATE],
            "property_mapping": ["assay", "assay_run_details"],
            "property_type": ["string", "datetime"],
        },
        "alias": "assay_run_details",
    },
    "assay_result": {
        "num_properties": 5,
        "identity_columns": {
            "property_names": [ASSAY_NAME, ASSAY_RUN_DATE, EPA_BATCH_ID],
            "property_mapping": ["assay", "assay_run_details", "batch_details"],
            "property_type": ["string", "datetime", "string"],
        },
        "alias": "assay_results",
    },
}
