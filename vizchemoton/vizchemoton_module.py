'''
Enric Petrus, December 2024. Added SCINE helper functiosn to link with the
amk-tools generation of html files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link
amk-tools and grrm-tools, generating interactive
HTML dashboards to visualize GRRM-generated reaction networks.
'''

# Standard Library Imports
from collections import Counter
import json
import copy

# Third-Party Library Imports
import yaml
import numpy as np
import bokeh.plotting
import bokeh.models as bkm
import RXVisualizer as arxviz
import networkx as nx
from xyz2mol import xyz2mol
from rdkit.Chem import MolToSmiles, MolFromSmiles

# Project-Specific SCINE imports
import scine_utilities as utils
import scine_database as db
from scine_chemoton.gears.pathfinder import Pathfinder as pf
from scine_database.energy_query_functions import (
    get_energy_change,
    get_barriers_for_elementary_step_by_type,
    get_energy_for_structure)


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


def get_crn_as_pathfinder(
        ip,
        port,
        db_name,
        dmethod,
        read_pathfinder=False,
        write_pathfinder=False,
        verbose=False):

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

    #########################################################################

    # # # Load Pathfinder and assign NetworkX Digraph
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

    # # # List of compounds and reactions
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

    if verbose:
        print("## Creating compounds and reaction objects")
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
                smiles = convert_struct_to_smile(
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
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = "//".join(
                html_compounds[cmp_dict[compound_id]]['_mongodb_id'])
            html_compounds[cmp_dict[compound_id]]['smiles'] = "//".join(
                html_compounds[cmp_dict[compound_id]]['_smiles'])

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
            smiles = convert_struct_to_smile(
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
            else:
                html_compounds[cmp_dict[compound_id]]['smiles'] = 'None'
    return html_reactions, html_compounds


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


def build_dashboard(G,title,outfile,size=(1400,800), layout_function=nx.kamada_kawai_layout,  map_field="energy", verbose=True):
    """
    Wrapper function to generate HTML visualizations for a given network.

    Input:
    - G (nx.Graph): object as generated from RXReader. For profile support, it should contain a graph["pathList"] property.
    - title (str): title for the visualization.
    - outfile (str): name of the output HTML file.
    - size (tuple): tuple of integers, size of the final visualization in pixels.
    - layout_function (nx.object, optional): Function to generate graph layout.
    - map_field (str): name of the field used for node coloring.

    Output:
    - lay (bokey.obj): Bokeh layout as generated by full_view_layout()
    """
    if verbose: print("## Writing {f1} output file".format(f1=outfile))

    ### Define sizing
    w1 = int(size[0]*4/7)
    w2 = int(size[0]*3/7)
    wu = int(size[0]/7)
    h = int(size[1]*6/8)

    sizing_dict = {'w1':w1,'w2':w2,'wu':wu,'h':h}

    ### Define custom classes

    style_template = """
    {% block postamble %}
	<script type="text/javascript" src="https://cdn.jsdelivr.net/gh/dgarayr/jsmol_to_bokeh/jsmol_to_bokeh.min.js"></script>
    <style>
    .bk-root .bk-btn-default {
        font-size: 1.2vh;
    }
    .bk-root .bk-input {
        font-size: 1.2vh;
        padding-bottom: 5px;
        padding-top: 5px;
    }
    .bk-root .bk {
        font-size: 1.2vh;
    }
    .bk-root .bk-clearfix{
        padding-bottom: 0.8vh;
    }
    </style>
    {% endblock %}
    """
    posx = layout_function(G)
    # Add model field to all nodes and edges & also vibrations
    arxviz.add_models(G)

    # Bokeh-powered visualization via RXVisualizer
    bk_fig,bk_graph = arxviz.bokeh_network_view(G,positions=posx,graph_title=title,width=w1,height=h,
                                                map_field=map_field,hide_energy=True)

    # bk_graph.selection_policy = bkm.NodesAndLinkedEdges()
    bk_graph.selection_policy = bkm.EdgesAndLinkedNodes()

    ### Modify the hovering tools here to add additional fields, removing the previous ones first
    valid_tools = [tool for tool in bk_fig.tools if tool.description]
    old_hovers = [tool for tool in valid_tools if "hover" in tool.description]
    for tool in old_hovers:
        bk_fig.tools.remove(tool)


    # custom edge hovering to reduce noise
    #
    hover_edgeJS = '''
    var erend = graph.edge_renderer.data_source
    var label1 = String.fromCharCode(916).concat("E1")
    var label2 = String.fromCharCode(916).concat("E2")
    if (cb_data.index.indices.length > 0) {
        var ndx = cb_data.index.indices[0]
        var tsname = erend.data["name"][ndx]
        if (tsname.includes('TSb')){
            hover.tooltips = [["tag","@name"]]
        } else {
            hover.tooltips = [["tag","@name"],["charge","@charge"],
                                ["multiplicity","@multiplicity"],["formula","@formula"],
                                [label1,"@deltaE1"],[label2,"@deltaE2"]]
        }
    }
    '''
        
    hide_barrlessJS = '''
        var erend = graph.edge_renderer.data_source
        var nrend = graph.node_renderer.data_source
        var edgenames = erend.data["name"]
        var nodenames = nrend.data["name"]
        var numEdges = edgenames.length
        var numNodes = nodenames.length
        var statusCounter = counter[0]
        var connectedNodes = []
        var labsNodes = figure.center[2].source.data
        
    
        if (statusCounter == 0){
        // remove and set 1
            for (let j = 0; j < numEdges; j++){
                var is_tsb = edgenames[j].includes("TSb")
                if (is_tsb) {
                    erend.data["start"][j] = null            
                    erend.data["end"][j] = null            
                }
                else {
                    connectedNodes.push(erend.data["start"][j])
                    connectedNodes.push(erend.data["end"][j])
                }
            }
            for (let i = 0; i < numNodes; i++){
                var nname = nodenames[i]
                if (!connectedNodes.includes(nname)) {
                    nrend.data["index"][i] = null
                    labsNodes["nnames"][i] = " "
                }
            }
            statusCounter = 1
        } else {
        // restore and set counter back to zero
             for (let j = 0; j < numEdges; j++){
                var is_tsb = edgenames[j].includes("TSb")
                if (is_tsb){
                    erend.data["start"][j] = backupEdgeRoutes["start"][j]            
                    erend.data["end"][j] = backupEdgeRoutes["end"][j]
                }
            }
            for (let i = 0; i < numNodes; i++) {
                if (nrend.data["index"][i] == null){
                    nrend.data["index"][i] = backupNodes["index"][i]
                    labsNodes["nnames"][i] = backupNodes["name"][i]
                }
            }
            statusCounter = 0
        }
        counter[0] = statusCounter
        nrend.change.emit()
        erend.change.emit()
        '''

    # Custom locateMolecule function to support search by SMILES
    locateMolecule = """
		// source - source object for JSMol
		// pass graph and fetch node and edge renderers
		// from fig, we modify x_range and y_range. Default plot starts from -1.2 to 1.2,
		var nrend = graph.node_renderer.data_source
		var erend = graph.edge_renderer.data_source
		var layout = graph.layout_provider.graph_layout
		// fetch the query in the data sources, choosing the appropiate renderer depending on the query
		var mol_query = text_input.value
		if (mol_query.includes("TS") || mol_query.includes("ts")) {
			var renderer = erend
			var other_renderer = nrend
		} else {
			var renderer = nrend
			var other_renderer = erend
		}
		var pool_names = renderer.data["name"]
		var pool_smiles = renderer.data["smiles"]

		// split species joined by + sign
		var pool_species = pool_names.reduce((acc,name) =>
					{acc.push(name.split("+"));
					return acc},[])
        // smiles are separated by // instead
		var pool_smiles_split = pool_smiles.reduce((acc,smiles) =>
					{acc.push(smiles.split("//"));
					return acc},[])

		// function to match results in the array
		var getSubstringIndices = function(arr,query){
			return arr.reduce(
					function(matches,tgt,i){
							if (tgt.includes(query))
						{matches.push(i)};
							return matches;
					},
					[]);
		}

		if (!mol_query.includes("+")) {
			var ndx1 = getSubstringIndices(pool_species,mol_query)
			var ndx2 = getSubstringIndices(pool_smiles_split,mol_query)
		} else {
			var ndx_u1 = pool_names.indexOf(mol_query)
			var ndx_u2 = pool_smiles.indexOf(mol_query)

			if (ndx_u1 < 0) {var ndx1 = []} 
            else {var ndx1 = [ndx_u1]}

            if (ndx_u2 < 0) {var ndx2 = []} 
            else {var ndx2 = [ndx_u2]}
		}
        
        // check both -> only choose the ones having matches, if both do, prefer SMILES
        if ((ndx1.length == 0) && (ndx2.length == 0)){
            var ndx = []
        } else if ((ndx1.length > 0) && (ndx2.length == 0)){
            var ndx = ndx1 
        } else if (ndx2.length > 0) {
            var ndx = ndx2
        }

		// locate positions of the node or of the nodes defining an edge
		if (mol_query.includes("TS") || mol_query.includes("ts")) {
			var n1 = renderer.data["start"][ndx]
			var n2 = renderer.data["end"][ndx]
			var pos1 = layout[n1]
			var pos2 = layout[n2]
			var positions = new Array(2)
			positions[0] = 0.5*(pos1[0]+pos2[0])
			positions[1] = 0.5*(pos1[1]+pos2[1])
		} else {
			var positions = layout[pool_names[ndx[0]]]
		}
		if (ndx.length > 0) {
			// clearing other sel. avoids problems for model loading sometimes
			other_renderer.selected.indices = []
			renderer.selected.indices = ndx
			fig.x_range.start = positions[0] - 0.5
			fig.x_range.end = positions[0] + 0.5
			fig.y_range.start = positions[1] - 0.5
			fig.y_range.end = positions[1] + 0.5
		}
		"""
    hover_node = bkm.HoverTool(description="Node hover",renderers=[bk_graph.node_renderer],
                               tooltips=[("tag","@name"),("charge","@charge"),("multiplicity","@multiplicity"),
                                         ("formula","@formula"),("smiles","@smiles")],
                               formatters={"@energy":"printf"})
    bk_fig.add_tools(hover_node)
    hover_edge = bkm.HoverTool(description="Edge hover",renderers=[bk_graph.edge_renderer],
                               formatters={"@energy":"printf"},line_policy="interp")
    hover_edge.callback = bkm.CustomJS(args={"hover":hover_edge,"graph":bk_graph},code=hover_edgeJS)
    bk_fig.add_tools(hover_edge)

    highl_callback = bkm.CustomJS(args={"graph":bk_graph}, code=arxviz.js_callback_dict["highlightNeighbors"])

    # We need edge backups
    edgesource = bk_graph.edge_renderer.data_source
    backup_edges = {"start":copy.deepcopy(edgesource.data["start"]),
					"end":copy.deepcopy(edgesource.data["end"])}
    backup_nodes = {"index":bk_graph.node_renderer.data_source.data["index"],
                    "name":bk_graph.node_renderer.data_source.data["name"]}

    hide_barrless_callback = bkm.CustomJS(args={"graph":bk_graph,"figure":bk_fig,"counter":[0],
                                                "backupNodes":backup_nodes,
                                                "backupEdgeRoutes":backup_edges}, code=hide_barrlessJS)

    lay = arxviz.full_view_layout(bk_fig,bk_graph,sizing_dict=sizing_dict)

    # add a button to the layout
    b_highlight = bkm.Button(label="Highlight neighbors",max_width=int(w1/6),align="center")
    b_highlight.js_on_click(highl_callback)
    b_hidebarrless = bkm.Button(label="Hide barrierless",max_width=int(w1/6),align="center")
    b_hidebarrless.js_on_click(hide_barrless_callback)


    sel_row = lay.children[0][0].children[2]
    sel_row.children[1].max_width = int(w1/6)
    sel_row.children = sel_row.children[0:2] + [b_highlight,b_hidebarrless] + [sel_row.children[-1]]

    # Modify the callback of the locate molecule button
    text_input = sel_row.children[0]
    js_mol_locator_nw = bkm.CustomJS(args = {"graph":bk_graph,"fig":bk_fig,"text_input":text_input},
								     code = locateMolecule)
    sel_button = sel_row.children[1]
    sel_button.js_event_callbacks['button_click'] = [js_mol_locator_nw]
    sel_button.js_on_click(js_mol_locator_nw)

    bokeh.plotting.output_file(outfile,title=title,mode="cdn")
    bokeh.plotting.save(lay,template=style_template)

    return lay,bk_fig,bk_graph

def scale_xyz_list(xyz, displ_vector=np.zeros(3)):
    """
    Bohr-to-angstrom scaling of a list of XYZ coordinates of the form
    [atom, [x, y, z]].

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with
    atom being a string and x,y,z floats.
    - displ_vector (np.Array, optional): for translating the geometry.

    Output:
    - xyz_nw (list): scaled XYZ coordinates in the same format as the input.
    """

    bohr_to_ang = 0.529
    xyz_arr = np.array([item[1] for item in xyz]) * bohr_to_ang + displ_vector
    xyz_nw = [[item[0], list(xyz_arr[ii])] for ii, item in enumerate(xyz)]
    return xyz_nw


def xyz_list_to_xyz_block(xyz):
    """
    Transform a list of xyz coordinates [atom, [x, y, z]] into a string block.

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with
    atom being a string and x,y,z floats.

    Output:
    - xyz_block (str): newline-joined block of the form
    a1,x1,y1,z2\na2,x2,y2,z2...
    """

    xyz_block = "\n".join(["%s %.6f %.6f %.6f" %
                           (item[0], *item[1]) for item in xyz])
    return xyz_block


def formula_from_xyz_block(xyz):
    """
    Generates the molecular formula for a given XYZ geometry.

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with
    atom being a string and x,y,z floats.

    Output:
    - formula (str): molecular formula from the input geometry.
    """
    labels = [item[0] for item in xyz]
    counter_list = sorted(Counter(labels).items())
    formula = ""
    for atom, ct in counter_list:
        if ct == 1:
            formula += atom
        else:
            formula += "%s%d" % (atom, ct)
    return formula


def sort_edge_names(edge_tuple):
    """
    Helper function to sort edge tuples lexicographically.

    Input:
    - edge_tuple (tuple): edge specification as a pair of node names.

    Output:
    - lexico_tuple (tuple): lexicographically sorted tuple.
    """
    n1, n2 = [int(nd) for nd in edge_tuple]
    srt_pair = sorted([n1, n2])
    lexico_tuple = tuple((str(nd) for nd in srt_pair))
    return lexico_tuple


def preprocess_compounds(compounds):
    """
    Helper function to process compounds properties.
    """
    tgt_vars = ["energy", "charge", "multiplicity"]
    for comp in compounds.values():
        for vv in tgt_vars:
            if not isinstance(comp[vv], list):
                comp[vv] = [comp[vv]]


def build_graph_edges(reaction_list):
    """
    Build graph edges.
    """
    return [(item[0], item[1], {"tsidx": item[2]}) for item in reaction_list]


def get_node_name_and_geometry(comp, dist_adduct, bohr_to_ang):
    """
    Helper function to retrieve the node name and geometry.
    """
    if isinstance(comp["crn_id"], list):
        node_name = "+".join(comp["crn_id"])
        xyz_list = comp["xyz"]
        xyz0_arr = np.array([item[1] for item in xyz_list[0]]) * bohr_to_ang
        cntr = xyz0_arr.mean(axis=0)
        xyz0 = [[item[0], list(xyz0_arr[ii])]
                for ii, item in enumerate(xyz_list[0])]
        xyz_full = xyz0

        for ii, xyz in enumerate(xyz_list[1:]):
            displ_vec = cntr + (ii + 1) * dist_adduct
            xyz_arr = np.array([it[1] for it in xyz]) * bohr_to_ang + displ_vec
            xyz_nw = [[item[0], list(xyz_arr[ii])]
                      for ii, item in enumerate(xyz)]
            xyz_full += xyz_nw
    else:
        node_name = comp["crn_id"]
        xyz_list = [comp["xyz"]]
        xyz_arr = np.array([item[1] for item in xyz_list[0]]) * bohr_to_ang
        xyz_full = [[item[0], list(xyz_arr[ii])]
                    for ii, item in enumerate(xyz_list[0])]

    return node_name, xyz_full, xyz_list


def add_node_attributes(graph, compounds, node_renaming, dist_adduct,
                        bohr_to_ang):
    """
    Add nodes attributes to the graph for building the HTML file.
    """
    for nd in graph.nodes(data=True):
        comp = compounds[nd[0]]
        tmp = get_node_name_and_geometry(comp, dist_adduct, bohr_to_ang)
        node_name, xyz_full, xyz_list = tmp
        node_renaming[nd[0]] = node_name
        strtmp = "%s %.6f %.6f %.6f"
        _xyz = "\n".join([strtmp % (item[0], *item[1]) for item in xyz_full])
        nd[1]["geometry"] = _xyz
        nd[1]["energy"] = sum(comp["energy"])
        nd[1]["ZPVE"] = 0.0
        nd[1]["name"] = node_name
        nd[1]["degree"] = graph.degree(nd[0])
        nd[1]["charge"] = "//".join([str(item) for item in comp["charge"]])
        tmpstr = [str(item) for item in comp["multiplicity"]]
        nd[1]["multiplicity"] = "//".join(tmpstr)
        tmpstr = [formula_from_xyz_block(xyz) for xyz in xyz_list]
        nd[1]["formula"] = "//".join(tmpstr)
        nd[1]["neighbors"] = list(graph.neighbors(nd[0]))
        nd[1]["smiles"] = comp.get("smiles", "None")


def add_edge_attributes(graph, compounds):
    """
    Add edges attributes to the graph for building the HTML file.
    """
    compounds_renamed = {}
    for _, comp in compounds.items():
        if isinstance(comp["crn_id"], list):
            new_key = "+".join(comp["crn_id"])
            compounds_renamed[new_key] = comp
        else:
            compounds_renamed[comp["crn_id"]] = comp
    for ii, ed in enumerate(graph.edges(data=True)):
        e1, e2 = [sum(compounds_renamed[nd]["energy"]) for nd in ed[0:2]]
        if ed[2]["tsidx"] == "None":
            e_ts = max(e1, e2)
            ed[2]["name"] = "TSb_%04d" % ii
            ed[2]["geometry"] = None
            ed[2]["energy"] = 0.0
            ed[2]["ZPVE"] = 0.0
            delta_e1 = (e_ts - e1, ed[0])
            delta_e2 = (e_ts - e2, ed[1])
            ed[2]["deltaE1"] = "%.2f (%s)" % delta_e1
            ed[2]["deltaE2"] = "%.2f (%s)" % delta_e2
            continue
        ts_compound = compounds[ed[2]["tsidx"]]
        xyz_list = [ts_compound["xyz"]]
        geom = scale_xyz_list(xyz_list[0])
        ed[2]["geometry"] = xyz_list_to_xyz_block(geom)
        ed[2]["name"] = ts_compound["crn_id"]
        e_ts = sum(ts_compound["energy"])
        delta_e1 = (e_ts - e1, ed[0])
        delta_e2 = (e_ts - e2, ed[1])
        ed[2]["deltaE1"] = "%.2f (%s)" % delta_e1
        ed[2]["deltaE2"] = "%.2f (%s)" % delta_e2
        ed[2]["energy"] = e_ts
        ed[2]["ZPVE"] = 0.0
        tmpstr = [str(item) for item in ts_compound["charge"]]
        ed[2]["charge"] = "//".join(tmpstr)
        tmpstr = [str(item) for item in ts_compound["multiplicity"]]
        ed[2]["multiplicity"] = "//".join(tmpstr)
        tmpstr = [formula_from_xyz_block(xyz) for xyz in xyz_list]
        ed[2]["formula"] = "//".join(tmpstr)


def process_graph(reaction_list, compounds, dist_adduct=3.0):
    """
    Wrapper function to generate a nx.Graph from a list of reactions and a
    dictionary of compounds, including XYZ-formatted geometries where
    individual geometries of the species forming adducts are joined.
    """
    bohr_to_ang = 0.529177
    graph = nx.Graph()
    edge_list = build_graph_edges(reaction_list)
    graph.add_edges_from(edge_list)
    node_renaming = {}
    preprocess_compounds(compounds)
    add_node_attributes(graph, compounds, node_renaming, dist_adduct,
                        bohr_to_ang)
    nx.relabel_nodes(graph, node_renaming, copy=False)

    # update neighbors after renaming
    for nd in graph.nodes(data=True):
        nd[1]["neighbors"] = list(graph.neighbors(nd[0]))

    add_edge_attributes(graph, compounds)

    return graph
