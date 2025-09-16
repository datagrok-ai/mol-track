import re


def create_column_mapping(output, level):
    """
    Create a mapping from output column names to actual field names in the response.

    The server returns field names like 'compounds_details_epa_compound_id'
    but the output might specify 'epa_compound_id' or 'compounds.details.epa_compound_id'.
    """
    col_map = {}

    # Derive prefixes directly from level
    prefix = level
    details_prefix = f"{level}_details"

    # Define direct fields for each level (fields that belong to the main table)
    direct_fields = {
        "compounds": [
            "id",
            "canonical_smiles",
            "inchi",
            "inchikey",
            "formula",
            "molregno",
            "hash_mol",
            "hash_tautomer",
            "hash_canonical_smiles",
            "hash_no_stereo_smiles",
            "hash_no_stereo_tautomer",
            "created_at",
            "updated_at",
            "is_archived",
            "structure",
        ],
        "batches": ["id", "compound_id", "batch_regno", "notes", "created_at", "updated_at"],
        "assay_results": ["id", "batch_id", "assay_run_id", "created_at", "updated_at", "created_by", "updated_by"],
    }

    level_direct_fields = direct_fields.get(level, [])

    for col in output:
        # Handle different input formats
        if "." in col:
            # Convert dot notation to underscore: 'compounds.details.epa compound id' -> 'compounds_details_epa_compound_id'
            # First replace dots with underscores, then replace spaces with underscores
            mapped_col = col.replace(".", "_").replace(" ", "_")
        else:
            # For simple column names, determine the appropriate prefix
            if col.startswith(f"{prefix}_"):
                # Already has the correct prefix
                mapped_col = col
            elif col in level_direct_fields:
                # Direct field from the main table
                mapped_col = f"{prefix}_{col}"
            else:
                # Assume it's a details field
                mapped_col = f"{details_prefix}_{col.replace(' ', '_')}"

        col_map[col] = mapped_col

    return col_map


def find_field_in_response(item, col_name, level):
    """
    Try to find a field in the response data using multiple possible field names.
    This is a fallback method when the standard mapping doesn't work.
    """
    if not item:
        return ""

    # Derive prefixes directly from level
    prefix = level
    details_prefix = f"{level}_details"

    # Compile regex pattern once for efficiency
    normalize_pattern = re.compile(r"[.\s]+")

    # Create normalized version of column name
    normalized_col = normalize_pattern.sub("_", col_name)

    # Create variations with exact case first
    variations = [
        col_name,  # Exact match
        normalized_col,  # Normalized version
        f"{prefix}_{normalized_col}",  # With prefix
        f"{details_prefix}_{normalized_col}",  # With details prefix
    ]

    # Try each variation with exact case first
    for variation in variations:
        if variation in item:
            return str(item[variation])

    # If no exact match found, try case-insensitive matching
    # Create a case-insensitive lookup dictionary
    item_lower = {k.lower(): k for k in item.keys()}

    # Try case-insensitive variations
    for variation in variations:
        variation_lower = variation.lower()
        if variation_lower in item_lower:
            # Return the value using the original case key
            return str(item[item_lower[variation_lower]])

    # If nothing found, return empty string
    return ""
