'''
Enric Petrus, August 2025. Cheminformatics helper functions.
'''

# Standard Library Imports
from collections import Counter,defaultdict
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

# Local Imports
from vizchemoton.html_module import (format_value_list)

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
    """
    Check if a SMILES string is valid using RDKit.
    """
    return Chem.MolFromSmiles(smiles) is not None

def _get_inchikey_from_smiles(smiles):
    """
    Transform SMILES to InChIKey.
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.inchi.MolToInchiKey(mol)

def get_public_database_id(name, smiles):
    """
    Wrapper
    """

    ddb = {name: None}
    if name == "pubchem":
        from .cheminfo_module import get_pubchem_cid
        ddb = get_pubchem_cid(smiles, delay=0.5, verbose=True)
        return ddb
    elif name == "chembl":
        from .cheminfo_module import get_chembl_id
        ddb = get_chembl_id(smiles)
        return ddb
    elif name == "chebi":
        from .cheminfo_module import get_chebi_id
        ddb = get_chebi_id(smiles)

    return ddb


def get_pubchem_cid(smiles, delay=0.5, verbose=True):
    """
    Check if SMILES are in PubChem.

    Args:
        smiles_list (list): List of SMILES strings.
        delay (float): Delay between API requests (default 0.5s).

    Returns:
        dict: {SMILES: True/False} indicating whether the compound exists in PubChem.
    """
    from pubchempy import get_compounds
    inchikey = _get_inchikey_from_smiles(smiles)
    time.sleep(delay)  # Avoid PubChem rate limits
    dpub = {"id": None}
    try:
        compounds = get_compounds(inchikey, 'inchikey')
        if len(compounds) > 0:  # True if found
            dpub["id"] = compounds[0].cid
    except Exception:
        dpub["id"] = "Error"  # PubChem rejected the request
    strtmp = "#### Querying PubChem. {s} has id = {b}"
    if verbose: print(strtmp.format(s=smiles, b=str(dpub["id"])))

    return dpub

def get_chembl_id(smiles, verbose=True):
    """
    Check if InChIKey is in ChEMBL database using their Python API.
    """
    from chembl_webresource_client.new_client import new_client
    inchikey = _get_inchikey_from_smiles(smiles)
    dchembl = {"id": None}
    try:
        molecule = new_client.molecule
        results = molecule.filter(molecule_structures__standard_inchi_key=inchikey)
        if len(results) > 0:
            idchembl = results[0]['molecule_chembl_id']
            number = int(''.join(filter(str.isdigit, idchembl)))
            dchembl["id"] = number
    except Exception:
        dchembl["id"] = "Error"
    strtmp = "#### Querying ChEMBL. {s} has id = {b}"
    if verbose: print(strtmp.format(s=smiles, b=str(dchembl["id"])))
    return dchembl

def get_chemspider_id(smiles, apikey, verbose=True):
    """
    Check if InChIKey is in ChemSpider database using their Python API.
    """
    from chemspipy import ChemSpider
    inchikey = _get_inchikey_from_smiles(smiles)
    try:
        cs = ChemSpider(apikey)
        results = cs.search(inchikey)
        dchemspi = {"id": None}
        if len(results) > 0:
            idchemspi = results[0].csid
            dchemspi["id"] = idchemspi
    except Exception:
        dchemspi["id"] = "Error"
    strtmp = "#### Querying ChemSpider. {s} has id = {b}"
    if verbose: print(strtmp.format(s=smiles, b=str(dchemspi["id"])))
    return dchemspi

def get_chebi_id(smiles, verbose=True):
   """
   Check if InChIKey is in ChEBI database using their Python API.
   """
   from libchebipy import search
   inchikey = _get_inchikey_from_smiles(smiles)
   dchebi = {"id": None}
   try:
       # search() returns a list of ChebiEntity objects
       entities = search(inchikey)
       if entities:
           idchebi = entities[0].get_id()  # take the first match
           number = int(''.join(filter(str.isdigit, idchebi)))
           dchebi["id"] = number
   except Exception as e:
       dchebi["id"] = "Error"
   strtmp = "#### Querying ChEBI. {s} has id = {b}"
   if verbose: print(strtmp.format(s=smiles, b=str(dchebi["id"])))
   return dchebi


def _get_rdkit_descriptor(mol, name, modules):
    """
    General function to extract properties from RDKit.
    """
    for module in modules:
        if hasattr(module, name):
            func = getattr(module, name)
            return func(mol)
    raise ValueError(f"Descriptor '{name}' not found in RDKit modules.")

def get_rdkit_properties(smiles, rdkitprop):
    """
    Returns the RDKit properties demanded in the config.yaml
    """
    modules = [Descriptors, Crippen, rdMolDescriptors]
    mol = Chem.MolFromSmiles(smiles)
    dprop = {k:_get_rdkit_descriptor(mol, k, modules) for k in rdkitprop}
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
        nd[1]["pubchemInfo"] = [str(pchm) for pchm in pubchem_info]
        nd[1]["pubchemInfoStr"] = "//".join(nd[1]["pubchemInfo"])

    return None

def compute_cheminf_props(graph,prop_keys):
    id_to_props = {}
    for nd in graph.nodes(data=True):
        id_list = nd[0].split("+")
        smiles_list = nd[1]["smiles"]
        current_elements = defaultdict(list)
        for spid,smi in zip(id_list,smiles_list):
            if smi == "None":
                current = {k:None for k in prop_keys}
            elif spid in id_to_props.keys(): 
                current = id_to_props[spid]
            else:
                current = get_rdkit_properties(smi,prop_keys)
                id_to_props[spid] = current
                
            for k in prop_keys:
                current_elements[k].append(current[k])
                # Apply these to the graph
                nd[1][k] = current_elements[k]
        # And prepare string formatting too
        for k in prop_keys:
            nd[1][k + "Str"] = format_value_list(current_elements[k])
    return None
