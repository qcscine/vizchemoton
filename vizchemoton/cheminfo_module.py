'''
Enric Petrus, August 2025. Cheminformatics helper functions.
'''

# Standard Library Imports
from collections import Counter
import json
import copy
import random

# Third-Party Library Imports
import yaml
import numpy as np
from xyz2mol import xyz2mol
from rdkit.Chem import MolToSmiles, MolFromSmiles, Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem import GetPeriodicTable
from rdkit import Chem

def get_cartesian_descriptors(xyz):
    """
    Custom function to extract basic cartesian descriptors for the post-clustering step.
    """
    atom_list = [a for a,_ in xyz]
    pt = GetPeriodicTable()
    atom_masses = np.array([pt.GetAtomicWeight(a) for a in atom_list])
    coords = [b for _,b in xyz]
    counts = Counter(atom_list)  # e.g., ['C', 'H', 'H', 'O']
    total_atoms = sum(counts.values())
    heavy_atoms = sum(counts[a] for a in counts if a != 'H')

    coords = np.array(coords)
    if atom_masses is None:
        atom_masses = np.ones(len(coords))  # default: equal mass

    total_mass = np.sum(atom_masses)
    center_of_mass = np.sum(coords * atom_masses[:, None], axis=0) / total_mass

    # Radius of gyration
    rg = np.sqrt(np.sum(atom_masses[:, None] * (coords - center_of_mass)**2) / total_mass)

    # Bounding box
    bbox = np.ptp(coords, axis=0)  # peak-to-peak along each axis
    descripcart = [total_atoms, heavy_atoms, rg, bbox[0] * bbox[1] * bbox[2]]
    return descripcart

def _convert_xyz_to_smiles(elements, coordinates, charge):
    """
    Convert xyz file to smiles using the external xyz2mol library.
    """
    # deprecated - now implemented in rdkit
    data = {'smiles': None}
    try:
        molformat = xyz2mol(elements, coordinates, charge, use_huckel=False,
                            embed_chiral=False, allow_charged_fragments=True)
        if len(molformat) != 0:
            smiles = MolToSmiles(molformat[0])
            m = MolFromSmiles(smiles)
            smiles = MolToSmiles(m)
            data = {'smiles': smiles}
            return data
        else:
            return data
    except:  # filter out cases where a SMILES is not feasible
        print("WARNING! Aggregate could not be converted to SMILES format")
        return data


def convert_struct_to_smile(centroid):
    """
    Convert structure instance to a smile.
    """
    conv2angs = 0.529177  # conversion of bohrs to anstrongs
    elements = [atom.value for atom in centroid.get_atoms().elements]
    coordinates = [[cj * conv2angs for cj in ci]
                   for ci in centroid.get_atoms().positions.tolist()]
    charge = centroid.get_charge()
    dsmiles = _convert_xyz_to_smiles(
        elements, coordinates, charge)
    return dsmiles

def is_valid_smiles(smiles):
    """Check if a SMILES string is valid using RDKit."""
    return Chem.MolFromSmiles(smiles) is not None

def query_pubchem(smiles_list, delay=0.5, verbose=True):
    """
    Check if SMILES are in PubChem.

    Args:
        smiles_list (list): List of SMILES strings.
        delay (float): Delay between API requests (default 0.5s).

    Returns:
        dict: {SMILES: True/False} indicating whether the compound exists in PubChem.
    """
    results = {}
    totalsmi = len(smiles_list)
    for ind, smi in enumerate(smiles_list):
        if not is_valid_smiles(smi):
            results[smi] = None  # Invalid SMILES, cannot be in PubChem
            continue
        time.sleep(delay)  # Avoid PubChem rate limits
        try:
            compounds = get_compounds(smi, 'smiles')
            if len(compounds) > 0:  # True if found
                results[smi] = compounds[0].cid
        except BadRequestError:
            results[smi] = None  # PubChem rejected the request
        if verbose: print(results[smi], smi)

    return results

def get_bio_properties(smiles):
    """
    Compute MW, logP, and TPSA for a list of SMILES.
    
    Returns a Pandas DataFrame.
    """
    mol = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    dprop = {"MolWt": mw, "LogP": logp, "TPSA": tpsa}
    return dprop
