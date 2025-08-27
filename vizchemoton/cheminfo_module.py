'''
Enric Petrus, August 2025. Cheminformatics helper functions.
'''

# Standard Library Imports
from collections import Counter
import json
import copy
import random
import time

# Third-Party Library Imports
import yaml
import numpy as np
from xyz2mol import xyz2mol
from rdkit.Chem import MolToSmiles, MolFromSmiles, Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem import GetPeriodicTable
from rdkit import Chem
from chembl_webresource_client.new_client import new_client
from pubchempy import get_compounds, BadRequestError

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

def get_pubchem_cid(smiles, delay=0.5, verbose=True):
    """
    Check if SMILES are in PubChem.

    Args:
        smiles_list (list): List of SMILES strings.
        delay (float): Delay between API requests (default 0.5s).

    Returns:
        dict: {SMILES: True/False} indicating whether the compound exists in PubChem.
    """
    time.sleep(delay)  # Avoid PubChem rate limits
    dpub = {"cid": None}
    try:
        compounds = get_compounds(smiles, 'smiles')
        if len(compounds) > 0:  # True if found
            dpub["cid"] = compounds[0].cid
    except BadRequestError:
        dpub["cid"] = None  # PubChem rejected the request
    if verbose: print("#### Querying PuChem. CID is " + str(dpub["cid"]))

    return dpub

def get_ChEMBL_id():
    """
    Check if InChIKey is in ChEMBL database using their Python API.
    """
    mol = Chem.MolFromSmiles(smiles)
    inchikey = Chem.inchi.MolToInchiKey(mol)
    # Query ChEMBL by InChIKey
    molecule = new_client.molecule
    results = molecule.filter(molecule_structures__standard_inchi_key=inchikey)
    if 'molecule_chembl_id' in results.keys():
        return 
    

def get_bio_properties(smiles):
    """
    Compute MW, logP, and TPSA for a list of SMILES.
    
    Returns a Pandas DataFrame.
    """
    mol = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    dprop = {"molwt": mw, "logp": logp, "tpsa": tpsa}
    return dprop

def pubchem_node_check(graph,compounds):
    """
    Checks whether the nodes in the graph have PubChem IDs, to state colors:
    0 - not present, 1 - some species present, 2 - all species present
    """
    pubchem_mapping = {v["crn_id"]:v["pubchem"] for k,v in compounds.items()}
    for nd in graph.nodes(data=True):
        pubchem_info = pubchem_mapping.get(nd[0],None)
        if not isinstance(pubchem_info,list):
            pubchem_info = [pubchem_info]
            
        if not pubchem_info:
            rnk = 0
        elif isinstance(pubchem_info,int):
            rnk = 2
        elif isinstance(pubchem_info,list):
            if all(pubchem_info):
                rnk = 2
            elif any(pubchem_info):
                rnk = 1
            else:
                rnk = 0
        
        nd[1]["pubchemRank"] = rnk 
        nd[1]["pubchemInfo"] = "//".join([str(pchm) for pchm in pubchem_info])

    return None
