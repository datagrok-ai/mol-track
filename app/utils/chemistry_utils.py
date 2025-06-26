import uuid
import os
import yaml

from rdkit import Chem
from rdkit.Chem import RegistrationHash
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.RegistrationHash import HashLayer

STANDARDIZER_CONFIG_FILE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "setup", "molecule_standarizer_operations.yaml")
)


def standardize_mol(
    mol: Chem.Mol,
    standardizer_config_file: str = STANDARDIZER_CONFIG_FILE,
) -> Chem.Mol:
    """
    Standardizes a given RDKit molecule using operations defined in a YAML configuration file.

    The operations are dynamically executed in the order defined in the file, but only if they are enabled.

    Args:
        mol (Chem.Mol): The molecule to standardize.
        standardizer_config_file (str): Path to the YAML configuration file with operation definitions.

    Returns:
        Chem.Mol: The standardized molecule after performing all configured operations.
    """

    with open(standardizer_config_file, "r") as f:
        config = yaml.safe_load(f)

    # Apply only the enabled operations in the order of declaration in the yml file.
    for operation in config.get("operations", []):
        operation_type = operation.get("type")
        is_enabled = operation.get("enable", True)

        if not is_enabled:
            continue

        if not operation_type:
            raise ValueError(f"Operation type is missing in the configuration file:{standardizer_config_file}.")

        mol = apply_standardizer_operation(mol, operation_type)

    return mol


def apply_standardizer_operation(mol: Chem.Mol, operation_type: str) -> Chem.Mol:
    """
    Applies a specific operation to the molecule based on the operation type.

    Args:
        mol (Chem.Mol): The molecule to modify.
        operation_type (str): The type of standardization operation to perform.

    Returns:
        Chem.Mol: The transformed molecule.
    """
    operation_map = {
        "Cleanup": rdMolStandardize.Cleanup,
        "FragmentParent": rdMolStandardize.FragmentParent,
        "RemoveHs": Chem.RemoveHs,
        "Uncharger": lambda mol: rdMolStandardize.Uncharger().uncharge(mol),
    }

    if operation_type not in operation_map:
        raise ValueError(f"Unknown operation type: {operation_type}")

    return operation_map[operation_type](mol)


def generate_hash_layers(mol: Chem.Mol) -> dict:
    """
    Generate layers for a given molecule.

    This function calculates the layers using the `RegistrationHash` module.

    Args:
        mol: An RDKit molecule object (`rdkit.Chem.Mol`) for which the layers
                  will be generated.

    Returns:
        dict: A dictionary containing the layers used to compute the MolHash.
    """

    return RegistrationHash.GetMolLayers(mol, enable_tautomer_hash_v2=True)


def generate_uuid_from_string(input_string: str) -> uuid.UUID:
    """
    Generate a UUID hash for a given input string, for hashing different molecule layers.

    Args:
        input_string (str): The input string to hash.

    Returns:
        uuid.UUID: The UUID hash of the input string, ready for PostgreSQL UUID type.
    """
    return uuid.uuid5(uuid.NAMESPACE_DNS, input_string)


def calculate_tautomer_hash(mol: Chem.Mol) -> str:
    """
    Calculate the tautomer hash for a given molecule.
    """
    return generate_uuid_from_string(generate_hash_layers(mol)[HashLayer.TAUTOMER_HASH])


def calculate_no_stereo_smiles_hash(mol: Chem.Mol) -> str:
    """Calculate the no-stereo SMILES hash for a given molecule."""

    return generate_uuid_from_string(generate_hash_layers(mol)[HashLayer.NO_STEREO_SMILES])


def calculate_no_stereo_tautomer_hash(mol: Chem.Mol) -> str:
    """
    Calculate the no-stereo tautomer hash for a given molecule.
    """
    return generate_uuid_from_string(generate_hash_layers(mol)[HashLayer.NO_STEREO_TAUTOMER_HASH])
