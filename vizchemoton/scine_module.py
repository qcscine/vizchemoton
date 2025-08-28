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
from .cheminfo_module import (get_cartesian_descriptors, convert_struct_to_smile, 
                              get_rdkit_properties, get_public_database_id)


def get_crn_as_pathfinder(
        ip,
        port,
        db_name,
        dmethod,
        read_pathfinder=False,
        write_pathfinder=False,
        verbose=False):
    """
    Helper function to read and write the chemical reaction network generated
    with SCINE.
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
    Returns a dictionary with the weight and stoichiometry of a structure.
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
    Applies an atom filter to disregard reactions in the final HTML file. 
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


def get_reactions_and_compounds(manager, pathfinder, dmethod, apikey,
                                calcsmiles, rdkitprop, databases, verbose=False):
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
    cmp_idx = 1
    for rxn_id in lhs_rxn_list:
        # Iterate through the reations of the network
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
                    node_ts = es_from_graph.get_transition_state().string()+";"
                    if node_ts not in cmp_dict.keys():
                        cmp_dict[node_ts] = cmp_idx
                        cmp_idx = cmp_idx + 1
                    html_reactions.append(
                        [cmp_dict[node_x], cmp_dict[node_y],
                         cmp_dict[node_ts]])

    # Create a dictionary for the compounds and their properties
    html_compounds = _get_html_compound_dict(pathfinder, model1, cmp_dict, structures, 
                                             compounds, flasks, properties, apikey, calcsmiles,
                                             rdkitprop, databases, verbose)
    return html_reactions, html_compounds

def _init_list_fields(rdkitprop):
    """Initialize all list-based fields for flask compounds."""
    return {k: [] for k in [
        "crn_id", "mongodb_id", "xyz", "charge", "multiplicity",
        "energy", "method", "basis_set", "program", "solvent", "solvation", 
        "smiles", "xyzdes", "pubchem", "chembl", "chebi", "chemspider"] + rdkitprop}

def _extract_structure_data(structure_obj, model, structures, properties, apikey, calcsmiles, rdkitprop, databases):
    """Extracts xyz, charge, multiplicity, energy, and model details from a structure object."""
    dprop = {k:None for k in rdkitprop}
    dpublidbs = {"pubchem": None, "chembl": None, "chebi": None, "chemspi": None}
    xyz = [(str(o.element), tuple(o.position))
           for o in structure_obj.get_atoms()]
    z, s = structure_obj.get_charge(), structure_obj.multiplicity
    dsmiles = convert_struct_to_smile(
        structure_obj) if calcsmiles else {'smiles': None}
    # smiles calculation
    if calcsmiles and dsmiles['smiles'] != None:
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
        "xyzdes": xyzdes,
        "pubchem": dpublidbs["pubchem"],
        "chembl": dpublidbs["chembl"],
        "chebi": dpublidbs["chebi"],
        "chemspider": dpublidbs["chemspi"],
    }
    # add dynamic requested rdkit properties
    for k in rdkitprop:
        tmpd[k] = dprop[k]
    return tmpd


def _get_compound_and_crnid(pathfinder, cmp_dict, mongoid, compounds, flasks):
    """
    Return the compound object and its crn_id
    """
    type_object = pathfinder.graph_handler.graph.nodes(data=True)[mongoid]["type"]
    if type_object == db.CompoundOrFlask.COMPOUND.name:
        compound = db.Compound(db.ID(mongoid), compounds)
        crn_id = "c" + str(cmp_dict[mongoid])
    else:
        compound = db.Flask(db.ID(mongoid), flasks)
        crn_id = "f" + str(cmp_dict[mongoid])

    return compound, crn_id

def _get_html_compound_dict(pathfinder, model1, cmp_dict, structures, compounds, flasks, properties, apikey, calcsmiles, rdkitprop, databases, verbose=False):    
    """
    Create a dictionary with the compounds and their chemical properties (xyz, charge, etc) for the chemical reaction    network.
    """
    if verbose:
        print("## Creating compounds and reaction objects")
    html_compounds = {}
    for compound_id in cmp_dict:
        compound_key = str(cmp_dict[compound_id])
        if "//" in compound_id:  # adducts of two aggregates 
            # if the user is interested in uploading the data in ioChem-BD,
            # this conditional block should be disregarded by deactivating
            # the following line:
            # continue
            ids = compound_id.split("//")
            html_compounds[compound_key] = _init_list_fields(rdkitprop)
            for _ids in ids:
                compound, crn_id = _get_compound_and_crnid(pathfinder, cmp_dict, _ids, compounds, flasks)
                structure = compound.get_centroid()
                structure_obj = db.Structure(structure, structures)
                struct_data = _extract_structure_data(structure_obj, model1, structures, properties, apikey, calcsmiles, rdkitprop, databases)
                html_compounds[compound_key]["crn_id"].append(crn_id)
                html_compounds[compound_key]["mongodb_id"].append(_ids)
                for k, v in struct_data.items():
                    html_compounds[compound_key][k].append(v)
            for s, k in [("+", "crn_id"), ("//", "mongodb_id"), ("//", "smiles")]:
                copy = html_compounds[compound_key][k].copy()
                tmpstr = s.join([str(o) for o in copy])
                html_compounds[compound_key][k] = tmpstr
            _sima, _simb = html_compounds[compound_key]['xyzdes']
            tmpchemsim = [np.mean(s) for s in zip(_sima, _simb)]
            html_compounds[compound_key]['xyzdes'] = tmpchemsim

        elif ";" in compound_id:  # transition state structure
            html_compounds[compound_key] = {}
            structure = compound_id[0:-1]
            structure_obj = db.Structure(db.ID(structure), structures)
            struct_data = _extract_structure_data(structure_obj, model1, structures, properties, apikey, calcsmiles, rdkitprop, databases)
            crn_id = "ts" + compound_key
            html_compounds[compound_key] = {
            **struct_data,
            "crn_id": crn_id,
            "mongodb_id": compound_id,}

        else:  # unimolecular reaction side
            compound, crn_id = _get_compound_and_crnid(pathfinder, cmp_dict, compound_id, compounds, flasks)
            structure = compound.get_centroid()
            structure_obj = db.Structure(structure, structures)
            struct_data = _extract_structure_data(structure_obj, model1, structures, properties, apikey, calcsmiles, rdkitprop, databases)
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
    Custom function to sort alphabetically a dictionary
    """
    return {k: d[k] for k in sorted(d)}


