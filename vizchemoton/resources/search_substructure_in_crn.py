"""
Basic script to perform substructure search given a custom SMILES/SMARTS and 
a compound JSON file generated with VizChemoton.
"""

import argparse
import json
import sys
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")  # optional

# Define positional and optional arguments
parser = argparse.ArgumentParser(description="Substructure search in a SMILES JSON dataset.")
parser.add_argument(
    "substructure",
    nargs="?",
    default="C=O",
    help="SMILES/SMARTS query string (default: C=O)"
)
parser.add_argument(
    "json_path",
    nargs="?",
    default="compounds_tme_dft.json",
    help="Path to JSON compounds file (default: compounds_tme_dft.json)"
)

args = parser.parse_args()

substructure = args.substructure
pathcompjson = args.json_path

# Load JSON file safely
try:
    with open(pathcompjson, "r") as f:
        compounds = json.load(f)
except FileNotFoundError:
    print(f"Error: File '{pathcompjson}' not found.")
    sys.exit(1)

compoundsmod = {v["smiles"]: v for v in compounds.values()}
keys = compoundsmod.keys()

# Look for substructure
query = Chem.MolFromSmarts(substructure)
if query is None:
    print(f"Error: Invalid SMARTS pattern '{substructure}'")
    sys.exit(1)

cntmatch = 0
for smi2 in keys:
    if not smi2 or "False" in str(smi2) or "None" in str(smi2):
        continue
    for fragment in smi2.split("."):
        for smik2 in fragment.split("//"):
            mol = Chem.MolFromSmiles(smik2)
            if mol is None:
                print(f"{smik2}: invalid SMILES")
                continue
            if mol.HasSubstructMatch(query):
                crnid = compoundsmod[smi2]["crn_id"]
                cntmatch += 1
                print(f"Found {substructure} in CRN ID = {crnid}")
print("---------------------------------------------------")
print(f"Substructure found in a total of {cntmatch} instances")
print("---------------------------------------------------")
