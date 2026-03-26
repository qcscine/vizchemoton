'''
Enric Petrus, August 2025. Collection of text functions to load or dump 
chemical reaction network data. 
Diego Garay-Ruiz, November 2023. Collection of helper functions to link
amk-tools and grrm-tools, generating interactive
HTML dashboards to visualize GRRM-generated reaction networks.
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
import networkx as nx
from rdkit.Chem import GetPeriodicTable

# Local Imports
from .cheminfo_module import (get_public_database_id, _convert_xyz_to_smiles,
                              get_rdkit_properties, get_cartesian_descriptors)

def vizchemoton_header():
    """
    Prints the VizChemoton ASCII art banner to the console.
    
    This function serves as the visual entry point for the CLI tool, 
    signaling the successful initialization of the package.
    """
    ascii_text = r"""
    __      ___      _____ _                          _
    \ \    / (_)    / ____| |                        | |
     \ \  / / _ ___| |    | |__   ___ _ __ ___   ___ | |_ ___  _ __
      \ \/ / | |_  / |    | '_ \ / _ \ '_ ` _ \ / _ \| __/ _ \| '_ \
       \  /  | |/ /| |____| | | |  __/ | | | | | (_) | || (_) | | | |
        \/   |_/___|\_____|_| |_|\___|_| |_| |_|\___/ \__\___/|_| |_|
    """
    print(ascii_text)


def load_config(config_file="config.yaml"):
    """
    Loads runtime parameters from a YAML configuration file.

    Args:
        config_file (str): Path to the configuration file. Defaults to "config.yaml".

    Returns:
        dict: Parsed configuration parameters.

    Raises:
        FileNotFoundError: If the specified config file does not exist.
        yaml.YAMLError: If the file contains invalid YAML syntax.
    """
    with open(config_file, "r") as file:
        return yaml.safe_load(file)


def custom_json_dump(obj, indent=2, level=0):
    """
    Recursively serializes a Python object to a JSON-formatted string with 
    selective compaction.

    Unlike standard `json.dumps`, this function preserves readability by:
    1. Pretty-printing dictionaries and lists of dictionaries.
    2. Collapsing simple primitive lists (like atomic coordinates or lists of 
       strings) into a single line to keep files compact.

    Args:
        obj (any): The Python object to serialize.
        indent (int): Number of spaces for indentation. Defaults to 2.
        level (int): Current recursion depth (used internally).

    Returns:
        str: A formatted JSON string.
        
    Example:
        Coordinates are kept on one line: `"xyz": [1.0, 2.0, 3.0]`
        Complex structures are expanded: 
        `"atoms": [
            {"symbol": "C", "id": 1},
            {"symbol": "O", "id": 2}
        ]`
    """
    space = ' ' * (indent * level)
    space_next = ' ' * (indent * (level + 1))

    if isinstance(obj, dict):
        items = []
        for k, v in obj.items():
            dumped_v = custom_json_dump(v, indent, level + 1)
            items.append(f'{space_next}"{k}": {dumped_v}')
        return '{\n' + ',\n'.join(items) + '\n' + space + '}'

    elif isinstance(obj, list):
        if not obj:
            return '[]'
        # Check if list of dicts → pretty print
        if all(isinstance(i, dict) for i in obj):
            items = [custom_json_dump(i, indent, level + 1) for i in obj]
            return '[\n' + ',\n'.join(space_next + item for item in items) + '\n' + space + ']'
        else:
            # Otherwise compact it
            return json.dumps(obj, separators=(',', ':'))

    else:
        return json.dumps(obj)


def write_compound_reactions_files(
        html_reactions,
        html_compounds,
        reaction_file,
        compound_file,
        verbose=True):
    """
    Serializes reaction connectivity and compound metadata to disk.

    This helper saves the elementary reactions as a line-delimited CSV and 
    the compound/transition state metadata as a formatted JSON file.

    Args:
        html_reactions (list of tuples): Elementary reactions defined as 
            (reactant_idx, product_idx, transition_state_idx).
        html_compounds (dict): Metadata mapping for all indices, including 
            coordinates, energy, and identifiers.
        reaction_file (str): Output path for the reaction connectivity file.
        compound_file (str): Output path for the JSON compound metadata.
        verbose (bool): If True, logs file creation status to the console.

    Returns:
        None: Writes data directly to the specified file paths.
    """ 
    if verbose:
        print(
            "## Writing {f1} and {f2} files".format(
                f1=reaction_file,
                f2=compound_file))
    # Open a file in write mode
    with open(reaction_file, 'w') as f:
        # Loop through the list and write each tuple to the file
        for item in html_reactions:
            r, p, ts = item
            f.write("{r},{p},{ts}\n".format(r=r, p=p, ts=ts))
    
    with open(compound_file, "w") as f:
        f.write(custom_json_dump(html_compounds, indent=2))


def read_compound_reactions_files(reaction_file, compounds_file, verbose=True):
    """
    Loads reaction connectivity and compound metadata from external files.

    Args:
        reaction_file (str): Path to the reaction CSV file (r,p,ts).
        compounds_file (str): Path to the JSON metadata file.
        verbose (bool): If True, logs the reading status to the console.

    Returns:
        tuple: A pair containing:
            - reaction_tuples (list of lists): The raw reaction connectivity.
            - compounds (dict): The parsed JSON metadata dictionary.

    Note:
        The `reaction_tuples` are returned as lists of strings from the 
        file read; further casting to integers may be required depending 
        on downstream graph construction logic.
    """
    if verbose:
        print(
            "## Reading {f1} and {f2} files".format(
                f1=reaction_file,
                f2=compounds_file))
    with open(reaction_file, "r") as freac:
        reaction_tuples = [line.strip().split(",")
                           for line in freac.readlines()]
    with open(compounds_file, "r") as fcomp:
        compounds = json.load(fcomp)

    return reaction_tuples, compounds


def review_compound_file(compounds_file, verbose=True, checkpoint_every=100):
    """
    Scans and repairs missing database identifiers in a compound metadata file.

    This function identifies entries where previous API queries failed (marked 
    as `False` or containing `False` in a list). It re-queries PubChem, ChEMBL, 
    and ChEBI to fill these gaps. To prevent data loss during long-running 
    network tasks, it implements an alternating "A/B" checkpointing strategy.
    
    Args:
        compounds_file (str): Path to the JSON file containing compound metadata.
        verbose (bool): If True, logs progress and API re-attempts to the console.
        checkpoint_every (int): Frequency of metadata serialization to 
            prevent data loss (number of compounds processed).

    Returns:
        dict: The updated compounds dictionary with repaired identifiers.

    Note:
        - Creates `.checkpointA`, `.checkpointB`, and final `.reviewed` files.
        - The "A/B" alternating logic ensures that if the system crashes 
          while writing a checkpoint, a valid backup from the previous 
          cycle remains available.
        - Handles both single SMILES and aggregate SMILES (delimited by '//').
    """
    if verbose:
        print(f"## Reviewing {compounds_file} file")

    with open(compounds_file, "r") as fcomp:
        compounds = json.load(fcomp)

    lkeys = ["pubchem", "chembl", "chebi"]
    total = len(compounds)
    last_checkpoint = time.time()
    flagcheckpoint = "A"
    for i, k1 in enumerate(compounds):
        cmp = compounds[k1]
        for k2 in lkeys:
            if cmp[k2] is False:
                smiles = compounds[k1]['smiles']
                if (smiles is not None) and ("None" not in smiles):
                    ddb = get_public_database_id(k2, smiles)
                    compounds[k1][k2] = ddb["id"]
                else:
                    compounds[k1][k2] = None
            elif isinstance(cmp[k2], list) and (False in cmp[k2]):
                lsmiles = compounds[k1]['smiles'].split("//")
                ltmp = []
                for smiles in lsmiles:
                    if (smiles != 'None') and (smiles is not None):
                        ddb = get_public_database_id(k2, smiles)
                        ltmp.append(ddb["id"])
                    else:
                        ltmp.append(None)
                compounds[k1][k2] = ltmp

        # Save checkpoint periodically
        if (i + 1) % checkpoint_every == 0:
            # Alternate checkpoints to have backup in case it crashes when editing
            # the checkpoint file.
            if flagcheckpoint == "A":
                checkpoint_file = compounds_file + ".checkpoint" + flagcheckpoint
                flagcheckpoint = "B"
            elif flagcheckpoint == "B":
                checkpoint_file = compounds_file + ".checkpoint" + flagcheckpoint
                flagcheckpoint = "A"
            
            with open(checkpoint_file, "w") as fcheckpoint:
                print(compounds[k1][k2])
                fcheckpoint.write(custom_json_dump(compounds, indent=2))
            
            if verbose:
                elapsed = time.time() - last_checkpoint
                print(f"Checkpoint saved after {i+1}/{total} compounds (elapsed: {elapsed:.1f}s)")
                last_checkpoint = time.time()

    # Final save
    reviewed_file = compounds_file + ".reviewed"
    with open(reviewed_file, "w") as f:
        f.write(custom_json_dump(compounds, indent=2))
    if verbose:
        print(f"Review completed. Output saved to {reviewed_file}")

    return compounds

def _add_smiles_to_compounds(xyz, charge):
    """
    Top-level wrapper for converting Cartesian coordinates to a SMILES dictionary.

    This function handles the conversion of atomic symbols to atomic numbers and 
    scales coordinates from Bohr to Angstroms ($0.529177$ factor) before passing 
    them to the inference engine.

    Args:
        xyz (list of tuple): List containing (symbol, [x, y, z]) coordinates in Bohr.
        charge (int): The total formal charge of the species.

    Returns:
        dict: A dictionary containing the inferred SMILES string.
            Example: {'smiles': 'CCO'}
    """
    ptable = GetPeriodicTable()
    elements = [ptable.GetAtomicNumber(a) for a,b in xyz]
    coordinates = [[b2* 0.529177 for b2 in b1] for a,b1 in xyz]
    dsmiles = convert_xyz_to_smiles(elements, coordinates, charge)
    return dsmiles

def _add_rdkit_properties(dsmiles, rdkitprop):
    """
    Appends calculated RDKit descriptors to an existing SMILES entry.

    Args:
        dsmiles (dict): Dictionary containing the 'smiles' key.
        rdkitprop (list of str): List of descriptor names to compute 
            (e.g., ['MolWt', 'LogP']).

    Returns:
        dict: A dictionary of descriptors. Returns a dictionary populated with 
              `None` values if the input SMILES is missing or invalid.
    """
    dprop = {k:None for k in rdkitprop}
    if dsmiles['smiles'] != None:
        dprop = get_rdkit_properties(dsmiles['smiles'], rdkitprop)
    
    return dprop

def _add_public_db_ids(dsmiles, databases):
    """
    Cross-references a SMILES string against selected public chemical databases.

    Args:
        dsmiles (dict): Dictionary containing the 'smiles' key.
        databases (dict): Configuration mapping database names to booleans, 
            e.g., {'pubchem': True, 'chembl': False}.

    Returns:
        dict: A dictionary of retrieved IDs (e.g., {'pubchem': 2244, ...}). 
              If a database is disabled in config or SMILES is missing, 
              the corresponding value remains `None`.
    """
    dpublidbs = {"pubchem": None, "chembl": None, "chebi": None}
    for name in databases.keys():
        if databases[name] and dsmiles['smiles'] != None:
            db_id = get_public_database_id(name, dsmiles['smiles'])["id"]
            dpublidbs[name] = db_id
    return dpublidbs


def upgrade_compound_file(compounds_file, rdkitprop, databases, verbose=True):
    """
    Upgrades a compound metadata file by calculating structural descriptors and 
    fetching missing database identifiers.

    This function iterates through a JSON compound file and fills in missing 
    information for each entry, including SMILES inference from XYZ coordinates, 
    RDKit physical properties, and cross-references to public databases. It 
    specifically handles "aggregate" species (mixtures/complexes) by processing 
    lists of coordinates and averaging Cartesian descriptors.

    Args:
        compounds_file (str): Path to the JSON file containing compound data.
        rdkitprop (list of str): RDKit descriptors to calculate (e.g., 'MolLogP').
        databases (dict): Configuration mapping database names to booleans for 
            API lookups (e.g., {'pubchem': True}).
        verbose (bool): If True, prints the progress of the file review.

    Returns:
        dict: The fully enriched compounds dictionary.

    Note:
        - The function saves a new file with the '.upgraded' extension.
        - For aggregate species (lists), SMILES are joined using the "//" delimiter.
        - Cartesian descriptors for aggregates are calculated as the mean 
          of the individual components.
    """
    if verbose:
            print("## Reviewing {f} file".format(
                    f=compounds_file))
    with open(compounds_file, "r") as fcomp:
        compounds = json.load(fcomp)
    # identify which compounds have errors
    for cnt, k1 in enumerate(compounds):
        #if cnt > 5:
        #    continue
        cmp = compounds[k1]
        if cmp == {}: # empty dict - artifact of submethods
            continue
        if isinstance(cmp["method"], str):
            dsmiles = _add_smiles_to_compounds(cmp["xyz"], cmp["charge"])
            smiles = dsmiles['smiles']
            dprop = _add_rdkit_properties(dsmiles, rdkitprop)
            for p in dprop:
                compounds[k1][p] = dprop[p]
            dpublidbs = _add_public_db_ids(dsmiles, databases)
            for d in dpublidbs:
                compounds[k1][d] = dpublidbs[d]
            xyzdes = get_cartesian_descriptors(cmp['xyz'])
        elif isinstance(cmp["method"], list):
            lsmiles, dpublidbs, lxyzdes = [], [], []
            for p in dprop:
                compounds[k1][p] = []
            for d in databases:
                compounds[k1][d] = []
            for xyz,charge in zip(cmp["xyz"], cmp["charge"]):
                dsmiles = _add_smiles_to_compounds(xyz, charge)
                lsmiles.append(dsmiles['smiles'])
                dpublidbs = _add_public_db_ids(dsmiles, databases)
                dprop = _add_rdkit_properties(dsmiles, rdkitprop)
                for p in dprop:
                    compounds[k1][p].append(dprop[p])
                for d in databases:
                    compounds[k1][d].append(dpublidbs[d])
                #dpublidbs = _add_public_db_ids(dsmiles, databases)
                xyzdesi = get_cartesian_descriptors(xyz)
                lxyzdes.append(xyzdesi)
            smiles = "//".join([str(o) for o in lsmiles])
            _sima, _simb = lxyzdes
            xyzdes = [np.mean(s) for s in zip(_sima, _simb)]
        compounds[k1]["xyzdes"] = xyzdes
        compounds[k1]["smiles"] = smiles
            
    with open(compounds_file+'.upgraded', "w") as f:
        f.write(custom_json_dump(compounds, indent=2))

    return compounds


def simplify_compounds(compounds):
    """
    Filters a compound dictionary to remove Transition State (TS) entries.

    Transition states often lack well-defined cheminformatic properties (like 
    standard SMILES or LogP). This function ensures downstream analysis 
    only processes stable intermediates.

    Args:
        compounds (dict): Dictionary mapping indices to compound metadata. 
            Expects a 'crn_id' key in the values to identify TS entries.

    Returns:
        dict: A "clean" dictionary containing only stable nodes (where 'ts' 
              is not in the 'crn_id').
    """
    compounds_clean = {k:v for k,v in compounds.items() if "ts" not in v["crn_id"]}
    return compounds_clean 

def process_compound_dbs(compounds,dblist=["pubchem","chebi","chembl"]):
    """
    Deconstructs aggregate nodes into a unique mapping of individual species.

    Reaction network nodes often represent aggregates (mixtures) of several 
    molecules. This function splits these aggregates using the '//' and '+' 
    delimiters, maps individual MongoDB IDs to their respective SMILES, 
    and aligns them with their public database identifiers.

    

    Args:
        compounds (dict): The master dictionary of reaction network compounds.
        dblist (list of str): The database fields to extract for each unique 
            species (e.g., ["pubchem", "chembl"]).

    Returns:
        dict: A mapping of unique MongoDB IDs to a flattened species dictionary:
            { 
              'mongodb_id': {
                  'smiles': str, 
                  'crn_id': str, 
                  'pubchem': int|None, ...
              }
            }
    """
    mongoid_mapping = {}
    compounds_clean = simplify_compounds(compounds)
    for cmp in compounds_clean.values():
        mongoid_list = cmp["mongodb_id"].split("//")
        smilist = str(cmp["smiles"]).split("//")
        cid_list = cmp["crn_id"].split("+")
        
        db_presence = {}
        for dbii in dblist:
            if isinstance(cmp[dbii],list):
                db_presence[dbii] = cmp[dbii]
            else:
                db_presence[dbii] = [cmp[dbii]]
                
        for ii,mid in enumerate(mongoid_list):
            if mid in mongoid_mapping:
                continue
                
            entry = {dbii:db_presence[dbii][ii] for dbii in dblist}
            entry.update({"smiles":smilist[ii],"crn_id":cid_list[ii]})
            mongoid_mapping[mid] = entry
    return mongoid_mapping

def count_matches(mongoid_mapping,error_flags=["None","Error","False",None,False]):
    """
    Calculates the frequency of valid metadata entries across the species mapping.

    This function audits the results of database cross-referencing and property 
    calculations. It identifies "valid" entries by excluding a customizable 
    list of error flags, providing a clear picture of data coverage for 
    each field (e.g., how many species successfully found a PubChem ID).

    Args:
        mongoid_mapping (dict): The species-level dictionary produced by 
            `process_compound_dbs`, where keys are MongoDB IDs.
        error_flags (list, optional): Values to be treated as missing or failed 
            data. Defaults to ["None", "Error", "False", None, False].

    Returns:
        dict: A summary dictionary where keys are field names (e.g., 'pubchem', 
            'smiles') and values are the integer counts of valid entries.

    Note:
        The function automatically determines fields to count by inspecting the 
        first entry in the mapping. It also includes a total count of unique 
        'mongodb_id' entries.
    """
    counts = {}
    # select a given entry in the dict to determine the keys to count automatically
    kr = next(iter(mongoid_mapping.keys()))
    field_list = mongoid_mapping[kr].keys()

    counts["mongodb_id"] = len(mongoid_mapping) 
    for field in field_list:
        values = [v[field] for v in mongoid_mapping.values() if v[field] not in error_flags]
        counts[field] = len(values)
    return counts

def flatten_xyz_structure(xyz_data):
    """
    Flattens nested xyz structure into a simple list:
    [(atom, x, y, z), ...]
    """
    flattened = []

    for block in xyz_data:
        element = block[0]
        coords = block[1]
        flattened.append((element, coords[0], coords[1], coords[2]))

    return flattened


def convert_bohr_to_angstrom(atoms, BOHR_TO_ANGSTROM=0.529177):
    """
    Converts coordinates from Bohr to Angstrom.
    """
    converted = []
    for element, x, y, z in atoms:
        converted.append(
            (
                element,
                x * BOHR_TO_ANGSTROM,
                y * BOHR_TO_ANGSTROM,
                z * BOHR_TO_ANGSTROM,
            )
        )
    return converted


def write_xyz(filename, atoms):
    """
    Writes atoms to XYZ file format.
    """
    with open(filename, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write("Converted from Bohr to Angstrom\n")
        for element, x, y, z in atoms:
            f.write(f"{element:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")

