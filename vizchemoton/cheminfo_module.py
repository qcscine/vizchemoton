'''
Enric Petrus, August 2025. Cheminformatics helper functions.
'''

# Standard Library Imports
from collections import Counter,defaultdict
import json
import copy
import random
import time
import requests
import ast
from datetime import datetime

# Third-Party Library Imports
import yaml
import numpy as np
from xyz2mol import xyz2mol
from rdkit.Chem import MolToSmiles, MolFromSmiles, Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem import GetPeriodicTable
from rdkit import Chem
import networkx as nx

#Local Imports
from vizchemoton.html_module import (format_value_list)
import scine_utilities as su
import scine_database as db
import scine_molassembler as masm

def get_cartesian_descriptors(xyz):
    """
    Calculate fundamental geometric and mass-based descriptors from Cartesian coordinates.

    This function computes structural metrics used for post-clustering analysis, 
    including atom counts, Radius of Gyration ($R_g$), and the bounding box volume.

    Args:
        xyz (list of tuple): A list where each element is a tuple of 
            (str: atomic_symbol, np.array: [x, y, z]). 
            Example: [('C', [0, 0, 0]), ('H', [0, 0, 1.08])]

    Returns:
        list: A list containing the following numerical descriptors:
            - total_atoms (int): Total count of all atoms.
            - heavy_atoms (int): Count of non-hydrogen atoms.
            - rg (float): Mass-weighted Radius of Gyration in Angstroms.
            - bbox_volume (float): Volume of the axis-aligned bounding box.

    Notes:
        The Radius of Gyration is calculated as:
        $$R_g = \sqrt{\frac{\sum_{i} m_i (\mathbf{r}_i - \mathbf{r}_{cm})^2}{\sum_{i} m_i}}$$
        where $m_i$ is the atomic weight and $\mathbf{r}_i$ are the coordinates.
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

def _convert_xyz_to_smiles(centroid):
    """
    Infers molecular connectivity (SMILES) from 3D Cartesian coordinates.

    Uses the `xyz2mol` algorithm to determine bond orders based on atomic distances 
    and total charge. It handles the conversion from Bohr to Angstroms and 
    canonicalizes the resulting SMILES string via RDKit.

    

    Args:
        centroid (Object): A molecular object containing:
            - .get_atoms().elements: List of atomic symbols or values.
            - .get_atoms().positions: NumPy array or list of coordinates in Bohr.
            - .get_charge(): Integer representing the total molecular charge.

    Returns:
        dict: A dictionary containing the 'smiles' key. 
            - If successful: {'smiles': 'C1=CC=CC=C1'} (canonical SMILES).
            - If conversion fails: {'smiles': None}.

    Note:
        - Bohr to Angstrom conversion factor used: $0.529177$.
        - This implementation is currently marked as deprecated in favor of 
          native RDKit XYZ-to-mol functionality.
        - Failure to infer a valid structure (e.g., non-physical distances) 
          will trigger a warning and return `None`.
    """
    # deprecated - now implemented in rdkit
    conv2angs = 0.529177  # conversion of bohrs to anstrongs
    elements = [atom.value for atom in centroid.get_atoms().elements]
    coordinates = [[cj * conv2angs for cj in ci]
                   for ci in centroid.get_atoms().positions.tolist()]
    charge = centroid.get_charge()
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

def _convert_scine_bo_to_smiles(centroid, properties, tmpfile, dsmiles):
    """
    Converts SCINE-derived bond orders and atom collections into a SMILES string.

    This function extracts a sparse bond order matrix from a SCINE database object,
    reconstructs the molecular topology, writes it to a temporary file, and 
    uses RDKit to generate a canonical SMILES string.

    Args:
        centroid (db.Structure): A SCINE database Structure object containing 
            atomic positions and bond order properties.
        properties (db.Collection): The SCINE properties collection used to 
            retrieve the bond order data.
        tmpfile (str): Path to a temporary file (e.g., .mol or .pdb) used for 
            intermediary topology storage.
        dsmiles (dict): A fallback dictionary (e.g., {'smiles': None}) to return 
            if the conversion fails early.

    Returns:
        dict: A dictionary containing the key 'smiles'. 
            - Returns the inferred SMILES string on success.
            - Returns the input `dsmiles` or `{'smiles': None}` on failure.

    Raises:
        RuntimeError: Caught internally if bond order properties are inaccessible.
        
    Note:
        Requires the `scine_utilities` (as `su`) and `scine_database` (as `db`) 
        modules, as well as RDKit's `Chem` module.
    """
    atomcollection = centroid.get_atoms()
    bonds = ast.literal_eval(centroid.get_graph("masm_idx_map"))
    if not centroid.has_property("bond_orders"):
        return dsmiles
    try:
        sparsitymatrix = centroid.get_property("bond_orders")
    except RuntimeError:
        return dsmiles
    prop_obj = db.Property(sparsitymatrix, properties)
    prop_json = prop_obj.json()
    data = json.loads(prop_json)
    rowidxs = data["data"]["row_idxs"]
    colidxs = data["data"]["col_idxs"]
    values = data["data"]["values"]
    bondcollection = su.BondOrderCollection(len(centroid.get_atoms()))
    for a, b, c in zip(rowidxs, colidxs, values):
        bondcollection.set_order(a, b, c)
    su.io.write_topology(tmpfile, atomcollection, bondcollection)
    rdkitmol = Chem.MolFromMolFile(tmpfile)
    try:
        smiles = Chem.MolToSmiles(rdkitmol)
    except:
        print("WARNING! Aggregate could not be converted to SMILES format")
        smiles = None
    dsmiles = {"smiles": smiles}
    return dsmiles

def _convert_to_smiles_molassembler(centroid, properties, tmpfile, dsmiles):
    """
    TO-DO
    """
    atomcollection = centroid.get_atoms()
    bonds = ast.literal_eval(centroid.get_graph("masm_idx_map"))
    if not centroid.has_property("bond_orders"):
        return dsmiles
    try:
        sparsitymatrix = centroid.get_property("bond_orders")
    except RuntimeError:
        return dsmiles
    prop_obj = db.Property(sparsitymatrix, properties)
    prop_json = prop_obj.json()
    data = json.loads(prop_json)
    rowidxs = data["data"]["row_idxs"]
    colidxs = data["data"]["col_idxs"]
    values = data["data"]["values"]
    bondcollection = su.BondOrderCollection(len(centroid.get_atoms()))
    for a, b, c in zip(rowidxs, colidxs, values):
        bondcollection.set_order(a, b, c)
    try:
        result = masm.interpret.molecules(atomcollection, bondcollection, masm.interpret.BondDiscretization.RoundToNearest)
        # If interpretation failed to yield a single molecule, try detecting bonds from scratch
        if len(result.molecules) != 1:
            bondcollection = utils.BondDetector.detect_bonds(atomcollection)
            result = masm.interpret.molecules(atomcollection, bondcollection, masm.interpret.BondDiscretization.RoundToNearest)
        if len(result.molecules) != 1:
            print("WARNING! Molassembler could not interpret a single molecule.")
            smiles = None
        else:
            mol = result.molecules[0]
            print("DEBUG", result.molecules, "centroid id", centroid.get_id())
            smiles = masm.io.experimental.emit_smiles(mol)
    except:
        print("WARNING! Aggregate could not be converted to SMILES format")
        smiles = None
    dsmiles = {"smiles": smiles}
    return dsmiles

def get_canolized_compid(centroid, properties, timestmp):
    """
    Generates a canonical identifier for a structure using Weisfeiler-Lehman graph hashing.

    This function extracts bond orders from a SCINE database object, constructs a 
    NetworkX graph where edges are weighted by bond orders, and computes a 
    topological hash. This serves as a robust 'Component ID' that is invariant 
    to atom indexing (canonicalization).

    Args:
        centroid (db.Structure): A SCINE database Structure object.
        properties (db.Collection): The SCINE properties collection containing 
            the 'bond_orders' sparse matrix.
        timestmp (float/str): A timestamp associated with the calculation 
            (currently unused in the function body).

    Returns:
        str: A hexadecimal string representing the Weisfeiler-Lehman graph hash.
        None: If the structure lacks 'bond_orders' or if the property lookup fails.

    Notes:
        - The graph $G = (V, E)$ is constructed where $V$ are atoms and $E$ are bonds.
        - The `edge_attr="bond"` ensures that bond orders (e.g., 1.0 vs 2.0) 
          result in distinct hashes.
        - This is often more computationally robust than SMILES canonicalization 
          for complex organometallic aggregates.
    """
    atomcollection = centroid.get_atoms()
    bonds = ast.literal_eval(centroid.get_graph("masm_idx_map"))
    if not centroid.has_property("bond_orders"):
        return None
    try:
        sparsitymatrix = centroid.get_property("bond_orders")
    except RuntimeError:
        return None
    prop_obj = db.Property(sparsitymatrix, properties)
    prop_json = prop_obj.json()
    data = json.loads(prop_json)
    rowidxs = data["data"]["row_idxs"]
    colidxs = data["data"]["col_idxs"]
    values = data["data"]["values"]
    x = []
    for a, b, c in zip(rowidxs, colidxs, values):
        x.append((a, b, {"bond": str(c)}))
    G = nx.Graph()
    G.add_edges_from(x)
    cancompid = nx.weisfeiler_lehman_graph_hash(G, edge_attr="bond")
    return cancompid

def convert_struct_to_smiles(centroid, properties, timestmp, multiplicity, smilesmode='hybrid'):
    """
    Orchestrates the conversion of 3D molecular structures to SMILES strings.

    This function provides a unified interface for structure-to-SMILES conversion
    using SCINE bond orders, the xyz2mol distance-based algorithm, or a hybrid 
    fallback approach optimized for singlet zwitterions.

    Args:
        centroid (db.Structure): SCINE database structure object.
        properties (db.Collection): SCINE properties collection for bond order lookup.
        timestmp (str): Unique identifier used to prevent collisions in temporary 
            file creation (e.g., "tmp{timestmp}.mol").
        multiplicity (int): The spin multiplicity of the molecule.
        smilesmode (str, optional): The conversion strategy to use. 
            Options include:
            - 'scine': Uses SCINE bond orders from Molassembler [1].
            - 'xyz2mol': Uses the distance-based algorithm by Kim et al. [2].
            - 'hybrid': Default. Uses 'scine' first, falling back to 'xyz2mol' 
              specifically for singlets when 'scine' fails [3].

    Returns:
        dict: A dictionary containing the key 'smiles'.
            - {'smiles': 'string'} on success.
            - {'smiles': None} if all attempted methods fail.

    References:
        [1] Brunken, C., & Reiher, M. J. Chem. Inf. Model. 2020, 60, 8, 3884–3900.
        [2] Kim, J. H., & Kim, H. S. Bull. Korean Chem. Soc. 2015, 36, 1769-1777.
        [3] Forthcoming: ChemRxiv 2026 (Internal hybrid refinement).
    """
    dsmiles = {"smiles": None}
    tmpfile = "tmp"+timestmp+".mol"
    if smilesmode == 'scine':
        dsmiles = _convert_scine_bo_to_smiles(centroid, properties, tmpfile, dsmiles)
    elif smilesmode == 'xyz2mol':
        dsmiles = _convert_xyz_to_smiles(centroid)
    elif smilesmode == 'hybrid': 
        dsmiles = _convert_scine_bo_to_smiles(centroid, properties, tmpfile, dsmiles)
        if (dsmiles["smiles"] is None) and (multiplicity == 1): 
            # xyz2mol handles better singlet zwitterions
            dsmiles = _convert_xyz_to_smiles(centroid)
    elif smilesmode == 'molassembler':
        dsmiles = _convert_to_smiles_molassembler(centroid, properties, tmpfile, dsmiles)
    return dsmiles

def is_valid_smiles(smiles):
    """
    Validates a SMILES string by attempting to sanitize it with RDKit.

    Args:
        smiles (str): The Simplified Molecular Input Line Entry System string.

    Returns:
        bool: True if RDKit can successfully parse and sanitize the molecule, 
              False otherwise.
    """
    return Chem.MolFromSmiles(smiles) is not None

def _get_inchikey_from_smiles(smiles):
    """
    Converts a SMILES string to a standard InChIKey.

    InChIKeys are fixed-length (27 character) hashes that are ideal for 
    database searching and avoiding rate-limit issues compared to long SMILES.

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        str: The resulting InChIKey (e.g., 'BSYREGRVZAWUOY-UHFFFAOYSA-N').
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.inchi.MolToInchiKey(mol)

def get_public_database_id(name, smiles):
    """
    Routes a SMILES query to specific public chemical databases.

    Acts as a central manager for external API calls, handling PubChem, 
    ChEMBL, and ChEBI lookups.

    Args:
        name (str): The database to query. Supported: 'pubchem', 'chembl', 'chebi'.
        smiles (str): The SMILES string to look up.

    Returns:
        dict: A dictionary containing the retrieved identifier.
            Example: {'id': 2244} or {'id': 'CHEMBL25'}.
            Returns {name: None} if the database name is not recognized.
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
        # deprecated from .cheminfo_module import get_chebi_id
        ddb = get_chebi_id(smiles)

    return ddb


def get_pubchem_cid(smiles, delay=0.5, verbose=True):
    """
    Retrieves the PubChem Compound ID (CID) for a given SMILES string.

    This function converts the SMILES to an InChIKey for more reliable 
    searching via the PubChemPy API and includes a sleep delay to comply 
    with PubChem's PUG-REST rate limits.

    Args:
        smiles (str): The SMILES string of the compound.
        delay (float): Time in seconds to wait before the API call. Default 0.5s.
        verbose (bool): If True, prints the query result to the console.

    Returns:
        dict: A dictionary containing the key 'id'.
            - If found: {'id': 12345} (int)
            - If not found: {'id': None}
            - If API error: {'id': 'Error'}
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
    Retrieves the numerical ChEMBL ID for a given SMILES string.

    Converts the SMILES to an InChIKey and queries the ChEMBL web resource 
    client. Note that ChEMBL IDs are typically returned as strings (e.g., 
    'CHEMBL25'); this function extracts the numeric digits and returns 
    them as an integer.

    Args:
        smiles (str): The SMILES string of the compound.
        verbose (bool): If True, prints the query status and result to the console.

    Returns:
        dict: A dictionary containing the key 'id'.
            - If found: {'id': 25} (int)
            - If not found: {'id': None}
            - If API error: {'id': 'Error'}

    Note:
        Requires the `chembl_webresource_client` package.
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

def get_chebi_id(smiles, verbose=True):
   """
    Placeholder for ChEBI ID retrieval (Currently Disabled).

    This function is a stub for the Chemical Entities of Biological 
    Interest (ChEBI) database lookup. It is currently deactivated 
    due to instability in the external API to prevent execution delays.

    Args:
        smiles (str): The SMILES string of the compound.
        verbose (bool): If True, prints a warning message to the console.

    Returns:
        dict: Always returns {'id': 'Error'} in its current state.
    """
   # to not lose time querying the URL
   dchebi = {"id": "Error"}
   print("Warning!: ChEBI deactivate due to problems with API. Returns Error without querying.")
   return dchebi

def _get_chebi_id(smiles, verbose=True):
   """
    Retrieves the ChEBI ID for a SMILES string via the EBI Elasticsearch API.

    This function bypasses the older `libchebipy` library in favor of a direct 
    REST API call. It searches for the compound using its InChIKey and returns 
    the integer ChEBI identifier.

    

    Args:
        smiles (str): The SMILES string of the compound.
        verbose (bool): If True, prints the query status and resulting ID.

    Returns:
        dict: A dictionary containing the key 'id'.
            - If found: {'id': 12345} (int)
            - If not found: {'id': None}
            - If connection/API failure: {'id': 'Error'}

    Notes:
        - API endpoint: https://www.ebi.ac.uk/chebi/backend/api/public/es_search/
        - Uses a broad search term and then validates the exact InChIKey match 
          within the result results to ensure accuracy.
    """
   # deprecated from libchebipy import search
   inchikey = _get_inchikey_from_smiles(smiles)
   dchebi = {"id": None}
   
   url = 'https://www.ebi.ac.uk/chebi/backend/api/public/es_search/'

   params = {
       'term': inchikey,
       'page': 1,
       'size': int(1e6)
   }

   headers = {
       'accept': '*/*'
   }

   try:
       response = requests.get(url, params=params, headers=headers)
       response.raise_for_status()
       data = response.json()
       for ent in data['results']:
           inchikey_i = ent['_source']['inchikey']
           if inchikey_i == inchikey:
               chebi_id = int(ent['_id'])
               dchebi['id'] = chebi_id
               break
   except requests.exceptions.RequestException as e:
       print(f"Error connecting to ChEBI ES API: {e}")
       dchebi["id"] = "Error"

   strtmp = "#### Querying ChEBI. {s} has id = {b}"
   if verbose: print(strtmp.format(s=smiles, b=str(dchebi["id"])))
   return dchebi


def _get_rdkit_descriptor(mol, name, modules):
    """
    Introspectively retrieves a descriptor value from a list of RDKit modules.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
        name (str): The exact name of the descriptor function (e.g., 'MolLogP').
        modules (list): A list of imported RDKit descriptor modules to search.

    Returns:
        Any: The calculated value of the descriptor.

    Raises:
        ValueError: If the descriptor name cannot be found in any of the 
            provided modules.
    """
    for module in modules:
        if hasattr(module, name):
            func = getattr(module, name)
            return func(mol)
    raise ValueError(f"Descriptor '{name}' not found in RDKit modules.")

def get_rdkit_properties(smiles, rdkitprop):
    """
    Calculates a set of chemical properties using RDKit modules.

    Parses a SMILES string and extracts multiple descriptors by searching through 
    the `Descriptors`, `Crippen`, and `rdMolDescriptors` modules.

    Args:
        smiles (str): The SMILES string of the molecule.
        rdkitprop (list of str): A list of descriptor names as defined in 
            config.yaml (e.g., ['MolLogP', 'HeavyAtomCount', 'TPSA']).

    Returns:
        dict: A mapping of property names to their calculated values.
            Example: {'MolLogP': 2.1, 'TPSA': 40.5}

    Note:
        This function assumes the SMILES string is valid. If RDKit cannot parse 
        the SMILES, it will likely raise an error during the mapping phase.
    """
    modules = [Descriptors, Crippen, rdMolDescriptors]
    mol = Chem.MolFromSmiles(smiles)
    dprop = {k:_get_rdkit_descriptor(mol, k, modules) for k in rdkitprop}
    return dprop

def pubchem_node_check(graph,compounds):
    """
    Checks nodes for PubChem IDs and assigns a representation rank for visualization.

    This is a specialized version of `db_node_check` focused specifically on 
    PubChem data. It updates the graph nodes in-place.

    Args:
        graph (networkx.Graph): The molecular graph or CRN.
        compounds (dict): A dictionary mapping compound keys to metadata, 
            containing 'crn_id' and 'pubchem' keys.

    Returns:
        None: Updates `graph` nodes with 'pubchemRank', 'pubchemInfo', 
              and 'pubchemInfoStr'.
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

def db_node_check(graph,compounds,db="pubchem"):
    """
    Enriches graph nodes with availability data from a specified database.

    Evaluates whether the species associated with a node exist in an external 
    database (e.g., PubChem, ChEMBL) and assigns a numerical rank (0, 1, or 2). 
    This rank is typically used to drive node color gradients in graph visualizations.

    Args:
        graph (networkx.Graph): The networkx graph object to be enriched.
        compounds (dict): Compound metadata dictionary. Expected schema per entry:
            {'crn_id': str, 'db_name': int | list | None}.
        db (str): The key for the database to check (e.g., "pubchem", "chembl").

    Returns:
        None: Updates graph nodes in-place with:
            - {db}Rank: '0' (none), '1' (partial), or '2' (complete).
            - {db}Info: List of retrieved IDs as strings.
            - {db}InfoStr: A '//' delimited string of IDs for tooltip display.
    """
    print("Processing nodes in %s" % db)
    db_mapping = {v["crn_id"]:v[db] for k,v in compounds.items()}
    for nd in graph.nodes(data=True):
        db_info = db_mapping.get(nd[0],None)
        if not isinstance(db_info,list):
            db_info = [db_info]
            
        if not db_info:
            rnk = 0
        elif isinstance(db_info,int):
            rnk = 2
        elif isinstance(db_info,list):
            if all(db_info):
                rnk = 2
            elif any(db_info):
                rnk = 1
            else:
                rnk = 0
        
        nd[1][db + "Rank"] = str(rnk) 
        nd[1][db + "Info"] = [str(pchm) for pchm in db_info]
        nd[1][db + "InfoStr"] = "//".join(nd[1][db + "Info"])

    return None

def add_multiple_dbs(graph,compounds,dblist=["pubchem","chembl","chebi","chemspider"]):
    """
    Iteratively runs database availability checks for a suite of external sources.

    This is a wrapper function to batch-process multiple database rankings 
    onto the same graph.

    Args:
        graph (networkx.Graph): The graph to update.
        compounds (dict): Compound metadata source.
        dblist (list of str): List of database keys to process.

    Returns:
        None: Graph is updated with ranking attributes for every database in dblist.
    """
    for db in dblist:
        db_node_check(graph,compounds,db)
    return None

def compute_cheminf_props(graph,prop_keys):
    """
    Computes and maps RDKit chemical properties to all nodes in a reaction graph.

    This function handles nodes that may contain multiple molecular species 
    (represented by '+' delimited IDs and a corresponding SMILES list). It 
    calculates the requested descriptors for each species, caches results to 
    optimize performance, and stores both raw lists and formatted strings 
    back onto the graph nodes.

    Args:
        graph (networkx.Graph): The graph where nodes contain 'smiles' lists.
        prop_keys (list of str): The RDKit descriptor names to calculate 
            (e.g., ['MolLogP', 'ExactMolWt']).

    Returns:
        None: Updates the graph in-place. Each node receives:
            - {prop_key}: A list of numerical values for each species in the node.
            - {prop_key}Str: A formatted string representation for tooltips.

    Notes:
        - If a SMILES string is "None", the property value is set to `None`.
        - Uses an internal dictionary (`id_to_props`) to cache results, ensuring 
          each unique species ID is only processed by RDKit once.
    """
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
