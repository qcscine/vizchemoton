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

# Local Imports
from .cheminfo_module import get_public_database_id

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


def review_compound_file(compounds_file, verbose=True):
    """
    Helper function which reviews the compound file in search for Error messages
    product of timeouts while querying the APIs of the Public Databases.
    """
    if verbose:
            print("## Reviewing {f} file".format(
                    f=compounds_file))
    with open(compounds_file, "r") as fcomp:
        compounds = json.load(fcomp)

    lkeys = ["pubchem", "chebi", "chembl"] #, "chemspider"]
    # identify which compounds have errors
    for k1 in compounds:
        cmp = compounds[k1]
        for k2 in lkeys:  # update json file
            if isinstance(cmp[k2], str) and cmp[k2] == 'Error':
                smiles = compounds[k1]['smiles']
                ddb = get_public_database_id(k2, smiles)
                compounds[k1][k2] = ddb["id"]
            elif isinstance(cmp[k2], list) and 'Error' in cmp[k2]:
                lsmiles = compounds[k1]['smiles'].split("//")
                ltmp = list()
                for smiles in lsmiles:
                    ddb = get_public_database_id(k2, smiles)
                    ltmp.append(ddb["id"])
                compounds[k1][k2] = ltmp

    with open(compounds_file+'.reviewed', "w") as f:
        f.write(custom_json_dump(compounds, indent=2))
 
    return compounds






