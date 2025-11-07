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

# Third-Party Library Imports
import yaml
import numpy as np
import networkx as nx
from rdkit.Chem import GetPeriodicTable

# Local Imports
from .cheminfo_module import (get_public_database_id, convert_xyz_to_smiles,
                              get_rdkit_properties, get_cartesian_descriptors)

def vizchemoton_header():
    """
    Plain text function to signal the start of VizChemoton
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
    Reads the configuration yaml file where the input parameters are defined.

    Input:
    - config_file (str): string with the name of the config file

    Output:
    - yaml (dict): dictionary with the input parameters
    """
    with open(config_file, "r") as file:
        return yaml.safe_load(file)


def custom_json_dump(obj, indent=2, level=0):
    """
    Recursively dump JSON with indent, compacting lists (like 'xyz') to a single line
    unless they contain dicts.
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
    Helper function to write reaction and compound files parsed from Chemoton.

    Input:
    - html_reactions (list): list of tuples of integers of the form [n1,n2,ts]
    specifying the indices of nodes and transition states from the set of
    compounds to define all elementary reactions in the network.
    - html_compounds (dict): dictionary mapping node/ts indices to the
    different computed fields that are available.

    Output:
    - reaction_file (str): path to the reactions file
    - compounds_file (str): path to the compounds file
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
    Helper function to read reaction and compound files parsed from Chemoton.

    Input:
    - reaction_file (str): path to the reactions file
    - compounds_file (str): path to the compounds file

    Output:
    - reaction_tuples (list): list of tuples of integers of the form
    [n1,n2,ts] specifying the indices of nodes and transition states from the
    set of compounds to define all elementary reactions in the network.
    - compounds (dict): dictionary mapping node/ts indices to the different
    computed fields that are available.
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


import json
import time

def review_compound_file(compounds_file, verbose=True, checkpoint_every=50):
    """
    Helper function which reviews the compound file in search for Error messages
    product of timeouts while querying the APIs of the Public Databases.
    Periodically writes checkpoints to avoid losing progress.
    """
    if verbose:
        print(f"## Reviewing {compounds_file} file")

    with open(compounds_file, "r") as fcomp:
        compounds = json.load(fcomp)

    lkeys = ["pubchem", "chebi", "chembl"]
    total = len(compounds)
    last_checkpoint = time.time()

    for i, k1 in enumerate(compounds):
        cmp = compounds[k1]
        for k2 in lkeys:
            if cmp[k2] == 'Error' or cmp[k2] == None:
                smiles = compounds[k1]['smiles']
                if smiles is not None:
                    ddb = get_public_database_id(k2, smiles)
                    compounds[k1][k2] = ddb["id"]
            elif isinstance(cmp[k2], list) and ('Error' in cmp[k2] or None in cmp[k2]):
                lsmiles = compounds[k1]['smiles'].split("//")
                ltmp = []
                for smiles in lsmiles:
                    if smiles != 'None':
                        ddb = get_public_database_id(k2, smiles)
                        ltmp.append(ddb["id"])
                compounds[k1][k2] = ltmp

        # Save checkpoint periodically
        if (i + 1) % checkpoint_every == 0:
            checkpoint_file = compounds_file + ".checkpoint"
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
    
    ptable = GetPeriodicTable()
    elements = [ptable.GetAtomicNumber(a) for a,b in xyz]
    coordinates = [[b2* 0.529177 for b2 in b1] for a,b1 in xyz]
    dsmiles = convert_xyz_to_smiles(elements, coordinates, charge)
    return dsmiles

def _add_rdkit_properties(dsmiles, rdkitprop):
    
    dprop = {k:None for k in rdkitprop}
    if dsmiles['smiles'] != None:
        dprop = get_rdkit_properties(dsmiles['smiles'], rdkitprop)
    
    return dprop

def _add_public_db_ids(dsmiles, databases):

    dpublidbs = {"pubchem": None, "chembl": None, "chebi": None, "chemspider": None}
    for name in databases.keys():
        if databases[name] and dsmiles['smiles'] != None:
            db_id = get_public_database_id(name, dsmiles['smiles'])["id"]
            dpublidbs[name] = db_id
    return dpublidbs


def upgrade_compound_file(compounds_file, rdkitprop, databases, verbose=True):
    """
    Helper function which reviews the compound file in search for Error messages
    product of timeouts while querying the APIs of the Public Databases.
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



