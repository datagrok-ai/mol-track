import uuid

import yaml
from rdkit import Chem
from rdkit.Chem import RegistrationHash
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.RegistrationHash import HashLayer



def standardize_mol(
    mol: Chem.Mol,
    standardizer_config_file: str = "molecule_standarizer_operations.yaml",
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
    # Generate molecular layers
    layers = RegistrationHash.GetMolLayers(mol, enable_tautomer_hash_v2=True)

    return layers


def generate_uuid_from_string(input_string: str) -> uuid.UUID:
    """
    Generate a UUID hash for a given input string, for hashing different molecule layers.

    Args:
        input_string (str): The input string to hash.

    Returns:
        uuid.UUID: The UUID hash of the input string, ready for PostgreSQL UUID type.
    """
    uuid_hash = uuid.uuid5(uuid.NAMESPACE_DNS, input_string)
    return uuid_hash

def generate_uuid_hash_mol(layers: dict) ->  uuid.UUID:
    """
    Generate a UUID hash for a given molecule's layers.

    Args:
        layers (dict): The layers of the molecule.

    Returns:
        str: The UUID hash of the molecule.
    """
    # Convert the layers to a string representation
    sorted_layers = tuple(sorted(layers.items(), key=lambda item: item[0].name))
    layers_str = str(sorted_layers)
    # Generate a UUID based on the string representation of the layers
    uuid_hash = uuid.uuid5(uuid.NAMESPACE_DNS, layers_str)
    return uuid_hash

def calculate_tautomer_hash(mol: Chem.Mol) -> str:
    """
    Calculate the tautomer hash for a given molecule.
    """
    layers = generate_hash_layers(mol)
    return generate_uuid_from_string(layers[HashLayer.TAUTOMER_HASH])

def calculate_no_stereo_smiles_hash(mol: Chem.Mol) -> str:
    """
    Calculate the no-stereo SMILES hash for a given molecule.
    """
    layers = generate_hash_layers(mol)
    return generate_uuid_from_string(layers[HashLayer.NO_STEREO_SMILES])

def calculate_no_stereo_tautomer_hash(mol: Chem.Mol) -> str:
    """
    Calculate the no-stereo tautomer hash for a given molecule.
    """
    # Generate a no-stereo version of the molecule
    layers = generate_hash_layers(mol)
    return generate_uuid_from_string(layers[HashLayer.NO_STEREO_TAUTOMER_HASH])

