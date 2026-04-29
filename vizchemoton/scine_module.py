'''
Enric Petrus, December 2024. Added SCINE helper function to link with the
amk-tools generation of HTML files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link
amk-tools and grrm-tools, generating interactive
HTML dashboards to visualize GRRM-generated reaction networks.
'''

# Standard Library Imports
from collections import Counter
import json
import copy
import random
from datetime import datetime

# Third-Party Library Imports
import yaml
import numpy as np
from xyz2mol import xyz2mol
from rdkit.Chem import MolToSmiles, MolFromSmiles, Descriptors
from rdkit.Chem import GetPeriodicTable

# Project-Specific SCINE imports
import scine_utilities as utils
import scine_database as db
from scine_chemoton.gears.pathfinder import Pathfinder as pf
from scine_database.energy_query_functions import (
    get_energy_change,
    get_barriers_for_elementary_step_by_type,
    get_energy_for_structure)
from .cheminfo_module import (get_cartesian_descriptors, convert_struct_to_smiles, 
                              get_rdkit_properties, get_public_database_id,
                              get_canolized_compid)


def get_crn_as_pathfinder(
        ip,
        port,
        db_name,
        dmethod,
        read_pathfinder=False,
        write_pathfinder=False,
        verbose=False):
    """
    Retrieve a chemical reaction network (CRN) as a Pathfinder object.

    This helper function connects to a SCINE database instance and
    optionally reads and/or writes a Pathfinder representation of the
    chemical reaction network (CRN).

    Parameters
    ----------
    ip : str
        IP address of the SCINE database server.
    port : int
        Port number of the SCINE database server.
    db_name : str
        Name of the database containing the reaction network.
    dmethod : str
        Database access method or backend identifier.
    read_pathfinder : bool, optional
        If True, load an existing Pathfinder object from storage
        instead of regenerating it from the database (default: False).
    write_pathfinder : bool, optional
        If True, write the generated Pathfinder object to storage
        for later reuse (default: False).
    verbose : bool, optional
        If True, enable verbose logging output (default: False).

    Returns
    -------
    pathfinder : Pathfinder
        The Pathfinder representation of the chemical reaction network.
    """
    
    # Connect to an active MongoDB
    manager = db.Manager()
    credentials = db.Credentials(ip, int(port), db_name)
    manager.set_credentials(credentials)
    if verbose:
        print("## Connecting to the Mongo-DB")
    manager.connect()
    model1 = db.Model(
        dmethod["method_family"],
        dmethod["method"],
        dmethod["basis_set"])
    model1.program = dmethod["program"]

    # Load Pathfinder and assign NetworkX Digraph
    pathfinder = pf(manager)

    if isinstance(read_pathfinder, str):
        if verbose:
            print("## Reading pathfinder object with name " + read_pathfinder)
        pathfinder.load_graph(read_pathfinder)
    elif isinstance(write_pathfinder, str):
        if verbose:
            print("## Writing pathfinder object with name " + write_pathfinder)
        pathfinder.options.model = model1
        pathfinder.options.graph_handler = "barrier"
        pathfinder.options.use_structure_model = True
        pathfinder.options.structure_model = model1
        pathfinder.build_graph()
        pathfinder.export_graph(write_pathfinder)

    return manager, pathfinder


def _calculate_weight(structure: db.Structure, structures: db.Collection, 
                      dstoich, verbose=False):
    """
    Compute the weight and stoichiometry of a given structure.

    This internal helper function determines the effective weight of a
    `db.Structure` object within a reaction network and returns the
    corresponding stoichiometric information. The calculation is based on
    the provided stoichiometric mapping and may access additional
    structures from the given collection.

    Parameters
    ----------
    structure : db.Structure
        The structure for which the weight and stoichiometry are computed.
    structures : db.Collection
        Collection containing structure objects, typically used to resolve
        references or retrieve related structural data.
    dstoich : dict
        Dictionary containing stoichiometric coefficients or mappings
        required for the weight calculation.
    verbose : bool, optional
        If True, print detailed diagnostic information during the
        calculation (default: False).

    Returns
    -------
    dict
        A dictionary containing:
        - "weight": The computed weight of the structure.
        - "stoichiometry": The associated stoichiometric information.
    """
    molec_dict = {}
    structure.link(structures)
    atoms = structure.get_atoms()
    weight, tmp = 0.0, list()
    for e in atoms.elements:
        weight += utils.ElementInfo.mass(e)
        tmp.append(str(e))
    molec_dict['weight'] = weight
    molec_dict['stoich'] = {d:tmp.count(d) for d in dstoich}
    if verbose: print(molec_dict['stoich'])
    return molec_dict


def check_natoms(reactants, reactants_type, compounds, flasks, structures, 
                 dstoich):
    """
    Filter reactions based on atom count consistency.

    This function applies an atom-balance filter to a set of reactants
    and removes reactions that do not satisfy predefined atom count
    criteria. It is primarily used to exclude chemically invalid or
    unbalanced reactions from the final HTML output.

    Parameters
    ----------
    reactants : list
        List of reactant identifiers involved in the reaction.
    reactants_type : list
        List describing the type of each reactant (e.g., compound, flask).
        Must correspond positionally to `reactants`.
    compounds : db.Collection
        Collection containing compound objects.
    flasks : db.Collection
        Collection containing flask objects.
    structures : db.Collection
        Collection containing structure objects used to retrieve
        atom count information.
    dstoich : dict
        Dictionary containing stoichiometric coefficients for the
        reaction participants.

    Returns
    -------
    bool
        True if the reaction passes the atom-count filter and should be
        included in the final output, False otherwise.
    """
    lhs, rhs = reactants
    lhst, rhst = reactants_type
    atomlist = dstoich.keys()
    condlist = []
    for lhsi, lhsti in zip(lhs, lhst):
        compound_id = lhsi.string()
        if lhsti.name == db.CompoundOrFlask.COMPOUND.name:
            compound = db.Compound(db.ID(compound_id), compounds)
        else:
            compound = db.Flask(db.ID(compound_id), flasks)
        structure = compound.get_centroid()
        molec_dict = _calculate_weight(db.Structure(structure), structures, 
                                       dstoich)
        weight, dstoich_i = molec_dict['weight'], molec_dict['stoich']
        for d in dstoich.keys():
           condlist.append(dstoich_i[d] < dstoich[d])
    for rhsi, rhsti in zip(rhs, rhst):
        compound_id = rhsi.string()
        if rhsti.name == db.CompoundOrFlask.COMPOUND.name:
            compound = db.Compound(db.ID(compound_id), compounds)
        else:
            compound = db.Flask(db.ID(compound_id), flasks)
        structure = compound.get_centroid()
        molec_dict = _calculate_weight(db.Structure(structure), structures,
                                       dstoich)
        weight, dstoich_j = molec_dict['weight'], molec_dict['stoich']
        for d in dstoich.keys():
           condlist.append(dstoich_j[d] < dstoich[d])
    if all(condlist):
        return True
    else:
        return False


def get_energy_and_barriers(
        energy_type,
        es_id,
        elementary_steps,
        model1,
        structures,
        properties,
        es_from_graph):
    """
    Wrapper function Gets the elementary step ID with the lowest energy of the
    corresponding transition state of a reaction.

    Input:
      - energy_type (str): name of the energy property such as
      'electronic_energy' or 'gibbs_free_energy'
      - es_id (str): id of the elementary_step
      - elementary_steps (db.Collection): the elementary step collection
      - model1 (dict): dictionary with the method_family, method, basis_set
      and program keys.
      - structures (db.Collection): the structures collection
      - properties (db.Collection): the properties step collection
      - es_from_graph (db.ElementaryStep): db object for the given es_id

    Returns:
      - energy (float): energy state in the reaction.
      - barriers (tuple): forward and backward energy barriers
      - not_none (bool): returns True if no None was found in the barriers
      tuple
    """
    energy = get_energy_change(
        db.ElementaryStep(
            es_id,
            elementary_steps),
        energy_type,
        model1,
        structures,
        properties)
    barriers = get_barriers_for_elementary_step_by_type(
        es_from_graph, energy_type, model1, structures, properties)

    if None in barriers:
        not_none = False
    else:
        not_none = True

    return energy, barriers, not_none


def get_reactions_and_compounds(manager, pathfinder, dmethod, 
                                calcsmiles, rdkitprop, databases, 
                                debugiter=False, verbose=True):
    """
    Extract the chemical reactions, compounds and transition states from the
    Mongo-DB where the exploration with Chemoton was run.

    Input:
      - db_name (str): name of the database
      - ip (str): internet protocol to the database
      - port (int): port number to the database
      - read_pathfinder (bool, str, optional): either False if no file to
      read, or string with the path to the file
      - write_pathfinder (bool, str, optional): either False if no file to
      write, or string with the path to the new file
      - verbose (bool, optional): if True, enable verbose output. Default
      is True.

    Returns:
      - html_reactions (list): a list of tuples with the indexes of the
      reactant, product, and TS.
      - html_compounds (dict): a dictionary for each compound containing
      relevant information (charge, spin, xyz ...)
    """
    # Get the SCINE collections
    model1 = db.Model(
        dmethod["method_family"],
        dmethod["method"],
        dmethod["basis_set"])
    model1.program = dmethod["program"]
    structures = manager.get_collection("structures")
    reactions = manager.get_collection("reactions")
    flasks = manager.get_collection("flasks")
    compounds = manager.get_collection("compounds")
    properties = manager.get_collection('properties')
    elementary_steps = manager.get_collection('elementary_steps')

    # List of compounds and reactions
    lhs_rxn_list = [
        node for node in pathfinder.graph_handler.graph.nodes if ";0;" in node]
    cmp_idx, numreac = 1, len(lhs_rxn_list)
    cmp_dict, html_reactions, html_compounds = {}, [], {}
    if verbose:
        tmpstr = "## Iterating through the {x} reactions in the network"
        print(tmpstr.format(x=str(numreac)))
    cmp_idx, rxn_idx = 1, 0
    for rxn_id in lhs_rxn_list:
        # Iterate through the reations of the network
        #print(rxn_id)
        tmpstr = '### Iteration {a} out of {b}'
        rxn_idx += 1
        if verbose: print(tmpstr.format(b=str(numreac), a=str(rxn_idx)))
        
        if debugiter is not False:
            # Useful for testing the whole VizChemoton workflow for large CRNs
            if rxn_idx > debugiter:
                print("### WARNING! Debug continue activated in scine_module")
                continue

        rxn = db.Reaction(db.ID(rxn_id[:-3]), reactions)
        reactants = rxn.get_reactants(db.Side.BOTH)
        reactants_type = rxn.get_reactant_types(db.Side.BOTH)
        lhs, rhs = reactants
        s_lhs, s_rhs = len(lhs), len(rhs)
        if dmethod["vfilter"]:
            vfilter = check_natoms(reactants, reactants_type, compounds,
                                    flasks, structures, dmethod["vfilter"])
        else:
            vfilter = True 

        if s_lhs < 3 and s_rhs < 3 and vfilter:
            # Get reactant indexes
            cmp_dict_keys = cmp_dict.keys()
            if len(reactants[0]) == 1:
                node_x = reactants[0][0].string()
                if node_x not in cmp_dict_keys:
                    cmp_dict[node_x] = cmp_idx
                    cmp_idx = cmp_idx + 1
            elif len(reactants[0]) == 2:
                for node_i in [o.string() for o in reactants[0]]:
                    if node_i not in cmp_dict_keys:
                        cmp_dict[node_i] = cmp_idx
                        cmp_idx = cmp_idx + 1
                # flasks (i.e., adducts) are depicted with //
                node_x = "//".join(sorted([o.string() for o in reactants[0]]))
                if node_x not in cmp_dict_keys:
                    cmp_dict[node_x] = cmp_idx
                    cmp_idx = cmp_idx + 1

            # Get product indexes
            if len(reactants[1]) == 1:
                node_y = reactants[1][0].string()
                if node_y not in cmp_dict_keys:
                    cmp_dict[node_y] = cmp_idx
                    cmp_idx = cmp_idx + 1
            elif len(reactants[1]) == 2:
                for node_i in [o.string() for o in reactants[1]]:
                    if node_i not in cmp_dict_keys:
                        cmp_dict[node_i] = cmp_idx
                        cmp_idx = cmp_idx + 1
                # flasks (i.e., adducts) are depicted with //
                node_y = "//".join(sorted([o.string() for o in reactants[1]]))
                if node_y not in cmp_dict_keys:
                    cmp_dict[node_y] = cmp_idx
                    cmp_idx = cmp_idx + 1

            # Get elementary steps and energies
            if "elementary_step_id" in pathfinder.graph_handler.graph.nodes(
                                       data=True)[rxn_id]:
                es_id = db.ID(pathfinder.graph_handler.graph.nodes(
                        data=True)[rxn_id]["elementary_step_id"])
                es_from_graph = db.ElementaryStep(es_id, elementary_steps)
                _energy, _, not_none = get_energy_and_barriers(
                    'electronic_energy', es_id, elementary_steps, model1,
                    structures, properties, es_from_graph)

                step_type = es_from_graph.get_type()
                is_barrierless = step_type == db.ElementaryStepType.BARRIERLESS
                if is_barrierless and not_none and _energy is not None:
                    html_reactions.append([cmp_dict[node_x],
                                           cmp_dict[node_y], None])
                elif not_none:
                    #node_ts = es_from_graph.get_transition_state().string()+";"
                    str_id = es_from_graph.get_transition_state().string()+";"
                    node_ts = str_id + "_" + rxn_id
                    if node_ts not in cmp_dict.keys():
                        cmp_dict[node_ts] = cmp_idx
                        cmp_idx = cmp_idx + 1
                    html_reactions.append(
                        [cmp_dict[node_x], cmp_dict[node_y],
                         cmp_dict[node_ts]])

        elif s_lhs == 3 or s_rhs == 3 and vfilter:
            # Get reactant indexes
            cmp_dict_keys = cmp_dict.keys()
            if len(reactants[0]) == 1:
                node_x = reactants[0][0].string()
                if node_x not in cmp_dict_keys:
                    cmp_dict[node_x] = cmp_idx
                    cmp_idx = cmp_idx + 1
            elif len(reactants[0]) == 3:
                for node_i in [o.string() for o in reactants[0]]:
                    if node_i not in cmp_dict_keys:
                        cmp_dict[node_i] = cmp_idx
                        cmp_idx = cmp_idx + 1
                # flasks (i.e., adducts) are depicted with //
                node_x = "//".join(sorted([o.string() for o in reactants[0]]))
                if node_x not in cmp_dict_keys:
                    cmp_dict[node_x] = cmp_idx
                    cmp_idx = cmp_idx + 1

            # Get product indexes
            if len(reactants[1]) == 1:
                node_y = reactants[1][0].string()
                if node_y not in cmp_dict_keys:
                    cmp_dict[node_y] = cmp_idx
                    cmp_idx = cmp_idx + 1
            elif len(reactants[1]) == 3:
                for node_i in [o.string() for o in reactants[1]]:
                    if node_i not in cmp_dict_keys:
                        cmp_dict[node_i] = cmp_idx
                        cmp_idx = cmp_idx + 1
                # flasks (i.e., adducts) are depicted with //
                node_y = "//".join(sorted([o.string() for o in reactants[1]]))
                if node_y not in cmp_dict_keys:
                    cmp_dict[node_y] = cmp_idx
                    cmp_idx = cmp_idx + 1

            # Get elementary steps and energies
            if "elementary_step_id" in pathfinder.graph_handler.graph.nodes(
                                       data=True)[rxn_id]:
                es_id = db.ID(pathfinder.graph_handler.graph.nodes(
                        data=True)[rxn_id]["elementary_step_id"])
                es_from_graph = db.ElementaryStep(es_id, elementary_steps)
                _energy, _, not_none = get_energy_and_barriers(
                    'electronic_energy', es_id, elementary_steps, model1,
                    structures, properties, es_from_graph)

                step_type = es_from_graph.get_type()
                is_barrierless = step_type == db.ElementaryStepType.BARRIERLESS
                if is_barrierless and not_none and _energy is not None:
                    html_reactions.append([cmp_dict[node_x],
                                           cmp_dict[node_y], None])
                elif not_none:
                    #node_ts = es_from_graph.get_transition_state().string()+";"
                    str_id = es_from_graph.get_transition_state().string()+";"
                    node_ts = str_id + "_" + rxn_id
                    if node_ts not in cmp_dict.keys():
                        cmp_dict[node_ts] = cmp_idx
                        cmp_idx = cmp_idx + 1
                    html_reactions.append(
                        [cmp_dict[node_x], cmp_dict[node_y],
                         cmp_dict[node_ts]])

    # Create a dictionary for the compounds and their properties
    html_compounds = _get_html_compound_dict(pathfinder, model1, cmp_dict, structures, 
                                             compounds, flasks, properties, calcsmiles,
                                             rdkitprop, databases, verbose)
    return html_reactions, html_compounds

def _init_list_fields(rdkitprop):
    """
    Initialize list-based fields for flask compound data storage.

    This internal helper function creates and returns a dictionary in which
    each predefined field name is mapped to an empty list. The resulting
    dictionary is used to accumulate compound-related data (e.g., identifiers,
    structural information, computed properties, and external references)
    for flask entries. Additional RDKit-derived property fields are appended
    dynamically.

    Parameters
    ----------
    rdkitprop : list of str
        List of RDKit property field names to include in the dictionary.
        Each entry will be initialized with an empty list.

    Returns
    -------
    dict
        Dictionary mapping field names to empty lists, ready to be populated
        with compound-specific data.
    """
    return {k: [] for k in [
        "crn_id", "mongodb_id", "xyz", "charge", "multiplicity", "can_id",
        "energy", "method", "basis_set", "program", "solvent", "solvation", 
        "smiles", "inchikey", "xyzdes", "pubchem", "chembl", "chebi"] + rdkitprop}

def _extract_structure_data(structure_obj, model, structures, properties, calcsmiles, rdkitprop, databases, timestmp):
    """
    Extracts detailed structural, energetic, chemical, and database information from a SCINE structure object.

    This function processes a molecular structure object to collect:
      - Cartesian coordinates (xyz)
      - Molecular charge and spin multiplicity
      - Electronic energy in kJ/mol
      - Quantum chemical model details (method, basis set, program, solvent)
      - SMILES representation and canonical compound ID
      - Cartesian-based molecular descriptors
      - Requested RDKit-derived properties
      - Identifiers from public chemical databases (PubChem, ChEMBL, ChEBI)

    Parameters
    ----------
    structure_obj : scine.structure.Structure
        SCINE structure object representing the molecular system.
    model : object
        Quantum chemistry model containing attributes `method`, `basis_set`, `program`, `version`, `solvent`, and `solvation`.
    structures : list
        List of structures used for reference in energy computations.
    properties : dict
        Dictionary of properties for computation and database mapping.
    calcsmiles : tuple(bool, bool)
        Tuple indicating whether to calculate SMILES (`bolsmiles`) and which type (`typsmiles`) to generate.
    rdkitprop : list of str
        List of RDKit property names to compute for the molecule.
    databases : dict
        Dictionary specifying which public chemical databases to query. Example: {"pubchem": True, "chembl": False, ...}.
    timestmp : datetime.datetime
        Timestamp to use for logging or record keeping.

    Returns
    -------
    dict
        A dictionary containing the following keys:
            - "xyz" : list of tuples, each containing element symbol and Cartesian coordinates
            - "charge" : int, molecular charge
            - "multiplicity" : int, spin multiplicity
            - "energy" : float, electronic energy in kJ/mol
            - "method" : str, quantum chemical method
            - "basis_set" : str, basis set used
            - "program" : str, program name and version
            - "solvent" : str, solvent used (if any)
            - "solvation" : str, solvation model (if any)
            - "smiles" : str or False, generated SMILES string
            - "can_id" : str, canonical compound identifier
            - "xyzdes" : dict, Cartesian-based molecular descriptors
            - "pubchem", "chembl", "chebi" : str or False, IDs from public databases
            - additional RDKit properties as requested in `rdkitprop`

    Notes
    -----
    - Energy is converted from Hartree to kJ/mol.
    - RDKit properties are computed only if SMILES generation is successful.
    - Database queries are performed only if the corresponding flag in `databases` is True.
    - The function handles both static properties (like xyz, energy) and dynamic properties requested at runtime.
    """
    dprop = {k:None for k in rdkitprop}
    dpublidbs = {"pubchem": False, "chembl": False, "chebi": False, "chemspi": False}
    xyz = [(str(o.element), tuple(o.position))
           for o in structure_obj.get_atoms()]
    z, s = structure_obj.get_charge(), structure_obj.multiplicity
    bolsmiles, typsmiles = calcsmiles
    dsmiles = convert_struct_to_smiles(
            structure_obj, properties, timestmp, s, typsmiles) if bolsmiles else {'smiles': False, "inchikey": False}
    cancompid = get_canolized_compid(structure_obj, properties, timestmp)
    # smiles calculation
    if bolsmiles and dsmiles['smiles'] != None:
        dprop = get_rdkit_properties(dsmiles['smiles'], rdkitprop)
        # query public databases
        for name in ["pubchem", "chembl", "chebi"]:  # hardcoded
            if databases[name]:
                db_id = get_public_database_id(name, dsmiles['smiles'])["id"]
                dpublidbs[name] = db_id
    e = get_energy_for_structure(
        structure_obj,
        'electronic_energy',
        model,
        structures,
        properties)
    if isinstance(e, (int, float)):
        e_kj = e * utils.KJPERMOL_PER_HARTREE
    else:
        e_kj = 0
    xyzdes = get_cartesian_descriptors(xyz)   
    # create static dictionary
    tmpd = {
        "xyz": xyz,
        "charge": z,
        "multiplicity": s,
        "energy": e_kj,
        "method": model.method,
        "basis_set": model.basis_set,
        "program": f"{model.program} {model.version}",
        "solvent": model.solvent,
        "solvation": model.solvation,
        "smiles": dsmiles['smiles'],
        "inchikey": dsmiles['inchikey'],
        "can_id": cancompid,
        "xyzdes": xyzdes,
        "pubchem": dpublidbs["pubchem"],
        "chembl": dpublidbs["chembl"],
        "chebi": dpublidbs["chebi"],
    }
    # add dynamic requested rdkit properties
    for k in rdkitprop:
        tmpd[k] = dprop[k]
    return tmpd


def _get_compound_and_crnid(pathfinder, cmp_dict, mongoid, compounds, flasks):
    """
    Retrieve the compound or flask object corresponding to a given MongoDB ID 
    and generate its CRN (Chemical Reaction Network) identifier.

    This function checks the type of the node in the graph associated with 
    `pathfinder`. Depending on whether the node represents a Compound or a Flask, 
    it initializes the appropriate object and constructs a CRN ID prefixed with 
    'c' for compounds and 'f' for flasks.

    Parameters
    ----------
    pathfinder : object
        Pathfinder instance containing the graph_handler with the chemical network graph.
    cmp_dict : dict
        Dictionary mapping MongoDB IDs to their numeric identifiers within the CRN.
    mongoid : str or int
        The MongoDB ID of the node (compound or flask) to retrieve.
    compounds : dict or collection
        Collection of Compound objects or data used to initialize a `db.Compound`.
    flasks : dict or collection
        Collection of Flask objects or data used to initialize a `db.Flask`.

    Returns
    -------
    tuple
        - compound : db.Compound or db.Flask
            The instantiated compound or flask object corresponding to `mongoid`.
        - crn_id : str
            The CRN identifier for the object, formatted as:
                - "c<number>" for compounds
                - "f<number>" for flasks

    Notes
    -----
    - The function relies on the node's "type" attribute in the graph to distinguish 
      between compounds and flasks.
    - CRN IDs are derived from `cmp_dict` which maps MongoDB IDs to numeric indices.
    """
    type_object = pathfinder.graph_handler.graph.nodes(data=True)[mongoid]["type"]
    if type_object == db.CompoundOrFlask.COMPOUND.name:
        compound = db.Compound(db.ID(mongoid), compounds)
        crn_id = "c" + str(cmp_dict[mongoid])
    else:
        compound = db.Flask(db.ID(mongoid), flasks)
        crn_id = "f" + str(cmp_dict[mongoid])

    return compound, crn_id

def _get_html_compound_dict(pathfinder, model1, cmp_dict, structures, compounds, flasks, properties, calcsmiles, rdkitprop, databases, verbose=False):    
    """
    Generate a dictionary containing chemical compounds and their computed properties for a chemical reaction network (CRN).

    This function iterates over all compounds in `cmp_dict` and constructs a dictionary 
    mapping each compound's key to its chemical and structural data. It handles three 
    types of compounds differently:
        1. Adducts of multiple aggregates (compound IDs containing "//")
        2. Transition state structures (compound IDs containing ";")
        3. Standard unimolecular compounds

    For each compound, it extracts:
        - Cartesian coordinates (xyz) and descriptors (xyzdes)
        - Charge, multiplicity, and energy
        - Quantum chemical model details (method, basis_set, program, solvent, solvation)
        - SMILES representation (if applicable) and canonical compound ID
        - Public database identifiers (PubChem, ChEMBL, ChEBI)
        - Additional RDKit properties as requested

    Parameters
    ----------
    pathfinder : object
        Pathfinder instance containing the CRN graph and chemical network data.
    model1 : object
        Quantum chemistry model containing method, basis set, program, version, solvent, and solvation attributes.
    cmp_dict : dict
        Dictionary mapping MongoDB IDs to numeric CRN identifiers.
    structures : list
        Collection of structures for energy and property computation.
    compounds : dict or collection
        Collection of Compound objects for CRN initialization.
    flasks : dict or collection
        Collection of Flask objects for CRN initialization.
    properties : dict
        Dictionary of properties used for structure analysis and database queries.
    calcsmiles : tuple(bool, str)
        Tuple indicating whether to compute SMILES (boolean) and the type of SMILES to generate.
    rdkitprop : list of str
        List of RDKit properties to compute for each molecule.
    databases : dict
        Flags indicating which public chemical databases to query for each compound.
    verbose : bool, optional
        If True, prints progress messages during dictionary construction (default is False).

    Returns
    -------
    dict
        A dictionary mapping each compound's CRN numeric key to a dictionary containing:
            - "xyz", "charge", "multiplicity", "energy"
            - "method", "basis_set", "program", "solvent", "solvation"
            - "smiles" and canonical compound ID ("can_id")
            - "xyzdes" descriptors
            - "pubchem", "chembl", "chebi" identifiers
            - Additional RDKit properties as specified
            - "crn_id" : CRN identifier string
            - "mongodb_id" : original MongoDB ID
        For adducts (multi-aggregate compounds), values are merged appropriately.

    Notes
    -----
    - Adducts with IDs containing "//" are combined, their SMILES and CRN IDs concatenated.
    - Transition states with ";" in their IDs do not generate SMILES.
    - Energy is converted to kJ/mol in `_extract_structure_data`.
    - The function uses `_get_compound_and_crnid` and `_extract_structure_data` for data extraction.
    - Returned dictionaries have keys sorted alphabetically via `_sort_dict_keys`.
    """
    if verbose:
        print("## Creating compounds and reaction objects")
    html_compounds = {}
    timestmp = datetime.now().isoformat()
    for compound_id in cmp_dict:
        compound_key = cmp_dict[compound_id]
        if "//" in compound_id:  # adducts of two aggregates 
            # if the user is interested in uploading the data in ioChem-BD,
            # this conditional block should be disregarded by deactivating
            # the following line:
            # continue
            ids = compound_id.split("//")
            html_compounds[int(compound_key)] = _init_list_fields(rdkitprop)
            for _ids in ids:
                compound, crn_id = _get_compound_and_crnid(pathfinder, cmp_dict, _ids, compounds, flasks)
                structure = compound.get_centroid()
                structure_obj = db.Structure(structure, structures)
                struct_data = _extract_structure_data(structure_obj, model1, structures, properties, calcsmiles, rdkitprop, databases, timestmp)
                html_compounds[compound_key]["crn_id"].append(crn_id)
                html_compounds[compound_key]["mongodb_id"].append(_ids)
                for k, v in struct_data.items():
                    html_compounds[compound_key][k].append(v)
            for s, k in [("+", "crn_id"), ("//", "mongodb_id"), ("//", "smiles"), ("//", "inchikey")]:
                copy = html_compounds[compound_key][k].copy()
                tmpstr = s.join([str(o) for o in copy])
                html_compounds[compound_key][k] = tmpstr
            #_sima, _simb = html_compounds[compound_key]['xyzdes']
            #tmpchemsim = [np.mean(s) for s in zip(_sima, _simb)]
            transposed = list(map(list, zip(*html_compounds[compound_key]['xyzdes'])))
            tmpchemsim = [np.mean(s) for s in transposed]
            assert len(tmpchemsim) == 4
            html_compounds[compound_key]['xyzdes'] = tmpchemsim

        elif ";" in compound_id:  # transition state structure
            # get structure and rxn ids
            compound_id, rxn_id = compound_id.split("_")
            html_compounds[compound_key] = {}
            structure = compound_id[0:-1]
            structure_obj = db.Structure(db.ID(structure), structures)
            _calcsmiles = (False, 'placeholder') # TSs do not need SMILES
            struct_data = _extract_structure_data(structure_obj, model1, structures, properties, _calcsmiles, rdkitprop, databases, timestmp)
            crn_id = "ts" + str(compound_key)
            # for ts the mongodb_id correspond to the rxn id
            html_compounds[compound_key] = {
            **struct_data,
            "crn_id": crn_id,
            "mongodb_id": rxn_id[:-3],}

        else:  # unimolecular reaction side
            compound, crn_id = _get_compound_and_crnid(pathfinder, cmp_dict, compound_id, compounds, flasks)
            structure = compound.get_centroid()
            structure_obj = db.Structure(structure, structures)
            struct_data = _extract_structure_data(structure_obj, model1, structures, properties, calcsmiles, rdkitprop, databases, timestmp)
            html_compounds[compound_key] = {
            **struct_data,
            "crn_id": crn_id,
            "mongodb_id": compound_id,
            }
        # Sort keys in alphabetical order
        tmpdict = html_compounds[compound_key].copy()
        html_compounds[compound_key] = _sort_dict_keys(tmpdict)
    return html_compounds

def _sort_dict_keys(d):
    """
    Return a new dictionary with keys sorted alphabetically.

    This function takes an input dictionary and produces a new dictionary 
    where the keys are in ascending alphabetical order. The original dictionary 
    is not modified.

    Parameters
    ----------
    d : dict
        The dictionary whose keys are to be sorted.

    Returns
    -------
    dict
        A new dictionary with the same key-value pairs as `d` but with keys sorted alphabetically.

    Notes
    -----
    - Only the top-level keys are sorted; nested dictionaries are not affected.
    - The ordering is determined by Python's default string comparison.
    """
    return {k: d[k] for k in sorted(d)}

