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
from sklearn.cluster import KMeans

# Project-Specific SCINE imports
import scine_utilities as utils
import scine_database as db
from scine_chemoton.gears.pathfinder import Pathfinder as pf
from scine_database.energy_query_functions import (
    get_energy_change,
    get_barriers_for_elementary_step_by_type,
    get_energy_for_structure)


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


def _convert_xyz_to_smiles(elements, coordinates, charge):
    """
    Convert xyz file to smiles using the external xyz2mol library.
    """
    # deprecated - now implemented in rdkit
    data = {'flag': False}
    try:  
        molformat = xyz2mol(elements, coordinates, charge, use_huckel=False,
                            embed_chiral=False, allow_charged_fragments=True)
        if len(molformat) != 0:
            smiles = MolToSmiles(molformat[0])
            m = MolFromSmiles(smiles)
            smiles = MolToSmiles(m)
            data = {'flag': True, 'smiles': smiles}
            return data
        else:
            return data
    except:  # filter out cases where a SMILES is not feasible
        print("WARNING! Aggregate could not be converted to SMILES format")
        return data


def _convert_struct_to_smile(centroid):
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


def get_reactions_and_compounds(manager, pathfinder, dmethod,
                                calcsmiles=False, verbose=False):
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

    html_compounds = _get_html_compound_dict(pathfinder, model1, cmp_dict, structures, compounds, flasks, properties, calcsmiles, verbose)
    return html_reactions, html_compounds

def _get_html_compound_dict(pathfinder, model1, cmp_dict, structures, compounds, flasks, properties, calcsmiles=False, verbose=False):    
    if verbose:
        print("## Creating compounds and reaction objects")
    html_compounds = {}
    for compound_id in cmp_dict:
        html_compounds[cmp_dict[compound_id]] = {}
        if "//" in compound_id:  # checking the flasks
            # if the user is interested in uploading the data in ioChem-BD,
            # this conditional block should be disregarded by deactivating
            # the following line:
            # continue
            ids = compound_id.split("//")
            html_compounds[cmp_dict[compound_id]]['crn_id'] = []
            html_compounds[cmp_dict[compound_id]]['_mongodb_id'] = []
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = []
            html_compounds[cmp_dict[compound_id]]['xyz'] = []
            html_compounds[cmp_dict[compound_id]]['charge'] = []
            html_compounds[cmp_dict[compound_id]]['multiplicity'] = []
            html_compounds[cmp_dict[compound_id]]['energy'] = []
            html_compounds[cmp_dict[compound_id]]['method'] = []
            html_compounds[cmp_dict[compound_id]]['basis_set'] = []
            html_compounds[cmp_dict[compound_id]]['program'] = []
            html_compounds[cmp_dict[compound_id]]['solvent'] = []
            html_compounds[cmp_dict[compound_id]]['solvation'] = []
            html_compounds[cmp_dict[compound_id]]['_smiles'] = []
            html_compounds[cmp_dict[compound_id]]['smiles'] = []
            html_compounds[cmp_dict[compound_id]]['_chemsim'] = []
            html_compounds[cmp_dict[compound_id]]['chemsim'] = []
            for _ids in ids:
                type_object = pathfinder.graph_handler.graph.nodes(data=True)[
                    _ids]["type"]
                if type_object == db.CompoundOrFlask.COMPOUND.name:
                    compound = db.Compound(db.ID(_ids), compounds)
                    html_compounds[cmp_dict[compound_id]]['crn_id'].append(
                        "c" + str(cmp_dict[_ids]))
                else:
                    compound = db.Flask(db.ID(_ids), flasks)
                    html_compounds[cmp_dict[compound_id]]['crn_id'].append(
                        "f" + str(cmp_dict[_ids]))
                structure = compound.get_centroid()
                structure_obj = db.Structure(structure, structures)
                xyz = [(str(o.element), tuple(o.position))
                       for o in structure_obj.get_atoms()]
                z, s = structure_obj.get_charge(), structure_obj.multiplicity
                smiles = _convert_struct_to_smile(
                    structure_obj) if calcsmiles else {'flag': False}
                e = get_energy_for_structure(
                    structure_obj,
                    'electronic_energy',
                    model1,
                    structures,
                    properties)
                if isinstance(e, (int, float)):
                    e_kj = e * utils.KJPERMOL_PER_HARTREE
                else:
                    e_kj = 0
                html_compounds[cmp_dict[compound_id]
                               ]['_mongodb_id'].append(_ids)
                html_compounds[cmp_dict[compound_id]]['xyz'].append(xyz)
                html_compounds[cmp_dict[compound_id]]['charge'].append(z)
                html_compounds[cmp_dict[compound_id]]['multiplicity'].append(s)
                html_compounds[cmp_dict[compound_id]]['energy'].append(e_kj)
                html_compounds[cmp_dict[compound_id]
                               ]['method'].append(model1.method)
                html_compounds[cmp_dict[compound_id]
                               ]['basis_set'].append(model1.basis_set)
                html_compounds[cmp_dict[compound_id]]['program'].append(
                    model1.program + " " + model1.version)
                html_compounds[cmp_dict[compound_id]
                               ]['solvent'].append(model1.solvent)
                html_compounds[cmp_dict[compound_id]
                               ]['solvation'].append(model1.solvation)
                if smiles['flag']:
                    html_compounds[cmp_dict[compound_id]
                                   ]['_smiles'].append(smiles['smiles'])
                else:
                    html_compounds[cmp_dict[compound_id]
                                   ]['_smiles'].append('None')
                html_compounds[cmp_dict[compound_id]]['_chemsim'].append(get_cartesian_descriptors(xyz))
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = "//".join(
                html_compounds[cmp_dict[compound_id]]['_mongodb_id'])
            html_compounds[cmp_dict[compound_id]]['smiles'] = "//".join(
                html_compounds[cmp_dict[compound_id]]['_smiles'])
            _sima, _simb = html_compounds[cmp_dict[compound_id]]['_chemsim']  # assuming always two
            html_compounds[cmp_dict[compound_id]]['chemsim'] = [np.mean(s) for s in zip(_sima, _simb)]

        elif ";" in compound_id:  # checking the transitions states
            structure = compound_id[0:-1]
            structure_obj = db.Structure(db.ID(structure), structures)
            xyz = [(str(o.element), tuple(o.position))
                   for o in structure_obj.get_atoms()]
            z, s = structure_obj.get_charge(), structure_obj.multiplicity
            e = get_energy_for_structure(
                structure_obj,
                'electronic_energy',
                model1,
                structures,
                properties)
            if isinstance(e, (int, float)):
                e_kj = e * utils.KJPERMOL_PER_HARTREE
            else:
                e_kj = 0
            html_compounds[cmp_dict[compound_id]
                           ]['crn_id'] = "ts" + str(cmp_dict[compound_id])
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = compound_id
            html_compounds[cmp_dict[compound_id]]['xyz'] = xyz
            html_compounds[cmp_dict[compound_id]]['charge'] = z
            html_compounds[cmp_dict[compound_id]]['multiplicity'] = s
            html_compounds[cmp_dict[compound_id]]['energy'] = e_kj
            html_compounds[cmp_dict[compound_id]]['method'] = model1.method
            html_compounds[cmp_dict[compound_id]
                           ]['basis_set'] = model1.basis_set
            html_compounds[cmp_dict[compound_id]
                           ]['program'] = model1.program + " 5.0.3"
            html_compounds[cmp_dict[compound_id]]['solvent'] = model1.solvent
            html_compounds[cmp_dict[compound_id]
                           ]['solvation'] = model1.solvation
            html_compounds[cmp_dict[compound_id]]['smiles'] = 'None'
            html_compounds[cmp_dict[compound_id]]['chemsim'] = 'None' 
        else:  # checking compounds
            type_object = pathfinder.graph_handler.graph.nodes(data=True)[
                compound_id]["type"]
            if type_object == db.CompoundOrFlask.COMPOUND.name:
                compound = db.Compound(db.ID(compound_id), compounds)
                html_compounds[cmp_dict[compound_id]
                               ]['crn_id'] = "c" + str(cmp_dict[compound_id])
            else:
                html_compounds[cmp_dict[compound_id]
                               ]['crn_id'] = "f" + str(cmp_dict[compound_id])
                compound = db.Flask(db.ID(compound_id), flasks)
            structure = compound.get_centroid()
            structure_obj = db.Structure(structure, structures)
            xyz = [(str(o.element), tuple(o.position))
                   for o in structure_obj.get_atoms()]
            z, s = structure_obj.get_charge(), structure_obj.multiplicity
            e = get_energy_for_structure(
                structure_obj,
                'electronic_energy',
                model1,
                structures,
                properties)
            smiles = _convert_struct_to_smile(
                structure_obj) if calcsmiles else {'flag': False}
            if isinstance(e, (int, float)):
                e_kj = e * utils.KJPERMOL_PER_HARTREE
            else:
                e_kj = 0
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = compound_id
            html_compounds[cmp_dict[compound_id]]['xyz'] = xyz
            html_compounds[cmp_dict[compound_id]]['charge'] = z
            html_compounds[cmp_dict[compound_id]]['multiplicity'] = s
            html_compounds[cmp_dict[compound_id]]['energy'] = e_kj
            html_compounds[cmp_dict[compound_id]
                           ]['method'] = model1.method  # model_obj.method
            # model_obj.basis_set
            html_compounds[cmp_dict[compound_id]
                           ]['basis_set'] = model1.basis_set
            # model_obj.program+" "+model_obj.version
            html_compounds[cmp_dict[compound_id]
                           ]['program'] = model1.program + " 5.0.3"
            html_compounds[cmp_dict[compound_id]
                           ]['solvent'] = model1.solvent  # model_obj.solvent
            # model_obj.solvation
            html_compounds[cmp_dict[compound_id]
                           ]['solvation'] = model1.solvation
            if smiles['flag']:
                html_compounds[cmp_dict[compound_id]]['smiles'] = smiles['smiles']
                #html_compounds[cmp_dict[compound_id]]['chemsim'] = smiles['chemsim']
            else:
                html_compounds[cmp_dict[compound_id]]['smiles'] = 'None'
                #html_compounds[cmp_dict[compound_id]]['chemsim'] = [0, 0, 0]
            html_compounds[cmp_dict[compound_id]]['chemsim'] = get_cartesian_descriptors(xyz)
    return html_compounds


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


