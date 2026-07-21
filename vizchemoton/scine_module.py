"""
Enric Petrus, December 2024. Added SCINE helper function to link with the
amk-tools generation of HTML files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link
amk-tools and grrm-tools, generating interactive
HTML dashboards to visualize GRRM-generated reaction networks.
"""

# Standard library imports
from datetime import datetime

# Third-party library imports
import numpy as np
from typing import Any, Dict

# Project-specific SCINE imports
import scine_utilities as utils
import scine_database as db
from scine_chemoton.gears.pathfinder import Pathfinder as pf
from scine_database.energy_query_functions import (
    get_energy_change,
    get_barriers_for_elementary_step_by_type,
    get_energy_for_structure,
)
from .cheminfo_module import (
    get_cartesian_descriptors,
    convert_struct_to_smiles,
    get_rdkit_properties,
    get_public_database_id,
)


def get_crn_as_pathfinder(
    ip,
    port,
    db_name,
    dmethod,
    pf_graph_new,
    pf_costs_new,
    verbose=False,
):
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
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!TO-DO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        print("## Connecting to the MongoDB")
    manager.connect()
    model = db.Model(
        dmethod["method_family"], dmethod["method"], dmethod["basis_set"]
    )
    model.program = dmethod["program"]
    if dmethod["solvent"] is not False:
        model.solvent = dmethod["solvent"]
    if dmethod["solvation"] is not False:
        model.solvation = dmethod["solvation"]

    # Load Pathfinder and assign NetworkX Digraph
    pathfinder = pf(manager)
    pf_graph_mode, pf_graph_file = pf_graph_new
    pf_costs_mode, pf_costs_file, pf_costs_init = pf_costs_new
    
    ## Load Pathfinder Graph
    if pf_graph_mode == "read":
         if verbose:
            print("## Reading pathfinder object with name " + pf_graph_file)
         pathfinder.load_graph(pf_graph_file)
    elif pf_graph_mode == "experimental.expand":
        # Experimental option which relies for the time being in a side branch of 
        # pathfinder containing the load_and_expand_graph() method. Needs to be 
        # tested and merged to the main Chemoton version. 
        if verbose:
            print("## (exp!) Expanding pathfinder object with name " + pf_graph_file)
        pathfinder.options.model = model
        pathfinder.options.use_structure_model = True
        pathfinder.options.structure_model = model
        pathfinder.load_and_expand_graph(pf_graph_file)
        print(pf_graph_file, pf_graph_file[:-5])
        pathfinder.export_graph(pf_graph_file[:-5]+"_expanded.json")
    elif pf_graph_mode == "write":
        if verbose:
            print("## Writing pathfinder object with name " + pf_graph_file)
        pathfinder.options.graph_handler = "barrier"
        pathfinder.options.model = model
        pathfinder.options.use_structure_model = True
        pathfinder.options.structure_model = model
        pathfinder.build_graph()
        pathfinder.export_graph(pf_graph_file)
    
    ## Load Pathfinder Costs
    if pf_costs_mode == "read":  # import previous compound costs
        if verbose:
            tmpstr = pf_graph_file + " " + pf_costs_file
            print("## Reading pathfinder objects with names " + tmpstr)
        pathfinder.load_graph(pf_graph_file, pf_costs_file)
    elif pf_costs_mode == "write":  # calculate compounds costs
        if verbose:
            print("## Writing pathfinder object with name " + pf_costs_file)
        pathfinder.set_start_conditions(pf_costs_init)
        pathfinder.calculate_compound_costs()
        pathfinder.update_graph_compound_costs()
        pathfinder.export_compound_costs()
        pathfinder.export_graph(pf_costs_file)
    elif pf_costs_mode == "ignore":  # ignore costs - add dummy values
        if verbose:
            print("## Ignoring the calculation of compound costs")
        for ni in pathfinder.graph_handler.graph.nodes:
            if ";" not in ni:  # not a rxn node
                pathfinder.compound_costs[ni] = 1

    return manager, pathfinder, model


def _calculate_weight(
    structure: db.Structure, structures: db.Collection, dstoich, verbose=False
):
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
    molec_dict: Dict[str, Any] = {}
    structure.link(structures)
    atoms = structure.get_atoms()
    weight, tmp = 0.0, list()
    for e in atoms.elements:
        weight += utils.ElementInfo.mass(e)
        tmp.append(str(e))
    molec_dict["weight"] = weight
    molec_dict["stoich"] = {d: tmp.count(d) for d in dstoich}
    if verbose:
        print(molec_dict["stoich"])
    return molec_dict


def get_energy_and_barriers(
    energy_type,
    es_id,
    elementary_steps,
    model,
    structures,
    properties,
    es_from_graph,
):
    """
    Wrapper function Gets the elementary step ID with the lowest energy of the
    corresponding transition state of a reaction.

    Input:
      - energy_type (str): name of the energy property such as
      'electronic_energy' or 'gibbs_free_energy'
      - es_id (str): id of the elementary_step
      - elementary_steps (db.Collection): the elementary step collection
      - model (dict): dictionary with the method_family, method, basis_set
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
        db.ElementaryStep(es_id, elementary_steps),
        energy_type,
        model,
        structures,
        properties,
    )
    barriers = get_barriers_for_elementary_step_by_type(
        es_from_graph, energy_type, model, structures, properties
    )

    if None in barriers:
        not_none = False
    else:
        not_none = True

    return energy, barriers, not_none


def get_reactions_and_compounds(
    manager,
    pathfinder,
    model,
    calcsmiles,
    rdkitprop,
    databases,
    debugiter=False,
    verbose=True,
):
    """
    Extract the chemical reactions, compounds and transition states from the
    Mongo-DB where the exploration with Chemoton was run.

    Input:
      - (TO-DO)
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
    structures = manager.get_collection("structures")
    reactions = manager.get_collection("reactions")
    flasks = manager.get_collection("flasks")
    compounds = manager.get_collection("compounds")
    properties = manager.get_collection("properties")
    elementary_steps = manager.get_collection("elementary_steps")

    # List of compounds and reactions
    lhs_rxn_list = [
        node for node in pathfinder.graph_handler.graph.nodes if ";0;" in node
    ]
    cmp_idx, numreac = 1, len(lhs_rxn_list)
    cmp_dict, html_reactions, html_compounds = {}, [], {}
    if verbose:
        tmpstr = "## Iterating through the {x} reactions in the network"
        print(tmpstr.format(x=str(numreac)))
    cmp_idx, rxn_idx = 1, 0
    for rxn_id in lhs_rxn_list:
        # Iterate through the reations of the network
        tmpstr = "### Iteration {a} out of {b}"
        rxn_idx += 1
        if verbose:
            print(tmpstr.format(b=str(numreac), a=str(rxn_idx)))

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
        # Get reactant indexes
        cmp_dict_keys = cmp_dict.keys()
        if len(reactants[0]) == 1:
            node_x = reactants[0][0].string()
            if node_x not in cmp_dict_keys:
                cmp_dict[node_x] = cmp_idx
                cmp_idx = cmp_idx + 1
        elif len(reactants[0]) >= 2:
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
        elif len(reactants[1]) >= 2:
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
        if (
            "elementary_step_id"
            in pathfinder.graph_handler.graph.nodes(data=True)[rxn_id]
        ):
            es_id = db.ID(
                pathfinder.graph_handler.graph.nodes(data=True)[rxn_id][
                    "elementary_step_id"
                ]
            )
            es_from_graph = db.ElementaryStep(es_id, elementary_steps)
            _energy, _, not_none = get_energy_and_barriers(
                "electronic_energy",
                es_id,
                elementary_steps,
                model,
                structures,
                properties,
                es_from_graph,
            )

            step_type = es_from_graph.get_type()
            is_barrierless = step_type == db.ElementaryStepType.BARRIERLESS
            if is_barrierless and not_none and _energy is not None:
                html_reactions.append(
                    [
                        cmp_dict[node_x],
                        cmp_dict[node_y],
                        None,
                        rxn_id[:-3],
                        es_id.string(),
                    ]
                )
            elif not_none:
                str_id = (
                    es_from_graph.get_transition_state().string() + ";"
                )
                node_ts = str_id + "_" + rxn_id
                if node_ts not in cmp_dict.keys():
                    cmp_dict[node_ts] = cmp_idx
                    cmp_idx = cmp_idx + 1
                html_reactions.append(
                    [
                        cmp_dict[node_x],
                        cmp_dict[node_y],
                        cmp_dict[node_ts],
                        rxn_id[:-3],
                        es_id.string(),
                    ]
                )

    # Create a dictionary for the compounds and their properties
    html_compounds = _get_html_compound_dict(
        pathfinder,
        model,
        cmp_dict,
        structures,
        compounds,
        flasks,
        properties,
        calcsmiles,
        rdkitprop,
        databases,
        verbose,
    )
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
    return {
        k: []
        for k in [
            "crn_id",
            "mongodb_id",
            "xyz",
            "charge",
            "multiplicity",
            "energy",
            "method",
            "basis_set",
            "program",
            "solvent",
            "solvation",
            "pfcost",
            "smiles",
            "inchikey",
            "xyzdes",
            "pubchem",
            "chembl",
            "chebi",
        ]
        + rdkitprop
    }


def _extract_structure_data(
    structure_obj,
    model,
    structures,
    properties,
    calcsmiles,
    rdkitprop,
    databases,
    timestmp,
):
    """
    Extracts detailed structural, energetic, chemical, and database information
    from a SCINE structure object.

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
        Quantum chemistry model containing attributes `method`, `basis_set`,
        `program`, `version`, `solvent`, and `solvation`.
    structures : list
        List of structures used for reference in energy computations.
    properties : dict
        Dictionary of properties for computation and database mapping.
    calcsmiles : tuple(bool, bool)
        Tuple indicating whether to calculate SMILES (`bolsmiles`) and which
        type (`typsmiles`) to generate.
    rdkitprop : list of str
        List of RDKit property names to compute for the molecule.
    databases : dict
        Dictionary specifying which public chemical databases to query.
        Example: {"pubchem": True, "chembl": False, ...}.
    timestmp : datetime.datetime
        Timestamp to use for logging or record keeping.

    Returns
    -------
    dict
        A dictionary containing the following keys:
            - "xyz" : list of tuples, each containing element symbol and
               Cartesian coordinates
            - "charge" : int, molecular charge
            - "multiplicity" : int, spin multiplicity
            - "energy" : float, electronic energy in kJ/mol
            - "method" : str, quantum chemical method
            - "basis_set" : str, basis set used
            - "program" : str, program name and version
            - "solvent" : str, solvent used (if any)
            - "solvation" : str, solvation model (if any)
            - "smiles" : str or False, generated SMILES string
            - "xyzdes" : dict, Cartesian-based molecular descriptors
            - "pubchem", "chembl", "chebi" : str or False, IDs from public
               databases
            - additional RDKit properties as requested in `rdkitprop`

    Notes
    -----
    - Energy is converted from Hartree to kJ/mol.
    - RDKit properties are computed only if SMILES generation is successful.
    - Database queries are performed only if the corresponding flag in
      `databases` is True.
    - The function handles both static properties (like xyz, energy) and
      dynamic properties requested at runtime.
    """
    dprop = {k: None for k in rdkitprop}
    dpublidbs = {
        "pubchem": False,
        "chembl": False,
        "chebi": False,
    }
    xyz = [
        (str(o.element), tuple(o.position)) for o in structure_obj.get_atoms()
    ]
    z, s = structure_obj.get_charge(), structure_obj.multiplicity
    bolsmiles, typsmiles = calcsmiles
    dsmiles = (
        convert_struct_to_smiles(
            structure_obj, properties, timestmp, s, typsmiles
        )
        if bolsmiles
        else {"smiles": False, "inchikey": False}
    )
    # smiles calculation
    if bolsmiles and dsmiles["smiles"] is not None:
        dprop = get_rdkit_properties(dsmiles["smiles"], rdkitprop)
        # query public databases
        for name in ["pubchem", "chembl", "chebi"]:  # hardcoded
            if databases[name]:
                db_id = get_public_database_id(name, dsmiles["smiles"])["id"]
                dpublidbs[name] = db_id
    e = get_energy_for_structure(
        structure_obj, "electronic_energy", model, structures, properties
    )
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
        "smiles": dsmiles["smiles"],
        "inchikey": dsmiles["inchikey"],
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
    `pathfinder`. Depending on if the node represents a Compound or a Flask,
    it initializes the appropriate object and constructs a CRN ID prefixed with
    'c' for compounds and 'f' for flasks.

    Parameters
    ----------
    pathfinder : object
        Pathfinder instance containing the graph_handler with the etwork graph.
    cmp_dict : dict
        Dictionary mapping MongoDB IDs to their numeric identifiers in the CRN.
    mongoid : str or int
        The MongoDB ID of the node (compound or flask) to retrieve.
    compounds : dict or collection
        Collection of Compound objects used to initialize a `db.Compound`.
    flasks : dict or collection
        Collection of Flask objects or data used to initialize a `db.Flask`.

    Returns
    -------
    tuple
        - compound : db.Compound or db.Flask
            The instantiated compound or flask corresponding to `mongoid`.
        - crn_id : str
            The CRN identifier for the object, formatted as:
                - "c<number>" for compounds
                - "f<number>" for flasks

    Notes
    -----
    - The function relies on the node's attribute in the graph to distinguish
      between compounds and flasks.
    - CRN IDs come from `cmp_dict` mapping MongoDB IDs to numeric indices.
    """
    type_object = pathfinder.graph_handler.graph.nodes(data=True)[mongoid][
        "type"
    ]
    if type_object == db.CompoundOrFlask.COMPOUND.name:
        compound = db.Compound(db.ID(mongoid), compounds)
        crn_id = "c" + str(cmp_dict[mongoid])
    else:
        compound = db.Flask(db.ID(mongoid), flasks)
        crn_id = "f" + str(cmp_dict[mongoid])

    return compound, crn_id


def _get_html_compound_dict(
    pathfinder,
    model,
    cmp_dict,
    structures,
    compounds,
    flasks,
    properties,
    calcsmiles,
    rdkitprop,
    databases,
    verbose=False,
):
    """
    Generate a dictionary containing chemical compounds and their computed
    properties for a chemical reaction network (CRN).

    This function iterates over all compounds in `cmp_dict` and constructs a
    dictionary mapping each compound's key to its chemical and structural data
    It handles three types of compounds differently:
        1. Adducts of multiple aggregates (compound IDs containing "//")
        2. Transition state structures (compound IDs containing ";")
        3. Standard unimolecular compounds

    For each compound, it extracts:
        - Cartesian coordinates (xyz) and descriptors (xyzdes)
        - Charge, multiplicity, and energy
        - Quantum model (method, basis_set, program, solvent, solvation)
        - SMILES representation (if applicable)
        - Public database identifiers (PubChem, ChEMBL, ChEBI)
        - Additional RDKit properties as requested

    Parameters
    ----------
    pathfinder : object
        Pathfinder instance containing the CRN graph and chemical network data.
    model : object
        Quantum model containing method, basis set, program, version, solvent,
        and solvation attributes.
    cmp_dict : dict
        Dictionary mapping MongoDB IDs to numeric CRN identifiers.
    structures : list
        Collection of structures for energy and property computation.
    compounds : dict or collection
        Collection of Compound objects for CRN initialization.
    flasks : dict or collection
        Collection of Flask objects for CRN initialization.
    properties : dict
        Dictionary of properties used for structure analysis and database
        queries.
    calcsmiles : tuple(bool, str)
        Tuple indicating whether to compute SMILES (boolean) and the type of
        SMILES to generate.
    rdkitprop : list of str
        List of RDKit properties to compute for each molecule.
    databases : dict
        Flags indicating which public chemical databases to query for each
        compound.
    verbose : bool, optional
        If True, prints progress messages during dictionary construction
        (default is False).

    Returns
    -------
    dict
        A dictionary mapping each compound's CRN numeric key to a dictionary
        containing:
            - "xyz", "charge", "multiplicity", "energy"
            - "method", "basis_set", "program", "solvent", "solvation"
            - "xyzdes" descriptors
            - "pubchem", "chembl", "chebi" identifiers
            - Additional RDKit properties as specified
            - "crn_id" : CRN identifier string
            - "mongodb_id" : original MongoDB ID
        For adducts, values are merged appropriately.

    Notes
    -----
    - Adducts with IDs containing "//" are combined, their SMILES and CRN IDs
      concatenated.
    - Transition states with ";" in their IDs do not generate SMILES.
    - Energy is converted to kJ/mol in `_extract_structure_data`.
    - The function uses `_get_compound_and_crnid` and `_extract_structure_data`
      for data extraction.
    - Returned dictionaries keys sorted alphabetically via `_sort_dict_keys`.
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
                compound, crn_id = _get_compound_and_crnid(
                    pathfinder, cmp_dict, _ids, compounds, flasks
                )
                tmpcost = _get_rescaled_pfcost(pathfinder, _ids)
                structure = compound.get_centroid()
                structure_obj = db.Structure(structure, structures)
                struct_data = _extract_structure_data(
                    structure_obj,
                    model,
                    structures,
                    properties,
                    calcsmiles,
                    rdkitprop,
                    databases,
                    timestmp,
                )
                html_compounds[compound_key]["crn_id"].append(crn_id)
                html_compounds[compound_key]["mongodb_id"].append(_ids)
                html_compounds[compound_key]["pfcost"].append(tmpcost)
                for k, v in struct_data.items():
                    html_compounds[compound_key][k].append(v)
            for s, k in [
                ("+", "crn_id"),
                ("//", "mongodb_id"),
                ("//", "smiles"),
                ("//", "inchikey"),
            ]:
                copy = html_compounds[compound_key][k].copy()
                tmpstr = s.join([str(o) for o in copy])
                html_compounds[compound_key][k] = tmpstr
            # _sima, _simb = html_compounds[compound_key]['xyzdes']
            # tmpchemsim = [np.mean(s) for s in zip(_sima, _simb)]
            transposed = list(
                map(list, zip(*html_compounds[compound_key]["xyzdes"]))
            )
            tmpchemsim = [np.mean(s) for s in transposed]
            html_compounds[compound_key]["xyzdes"] = tmpchemsim
            tmpcost = np.mean(html_compounds[compound_key]["pfcost"])
            html_compounds[compound_key]["pfcost"] = round(tmpcost, 0)


        elif ";" in compound_id:  # transition state structure
            # get structure and rxn ids
            compound_id, rxn_id = compound_id.split("_")
            html_compounds[compound_key] = {}
            structure = compound_id[0:-1]
            structure_obj = db.Structure(db.ID(structure), structures)
            _calcsmiles = (False, "placeholder")  # TSs do not need SMILES
            struct_data = _extract_structure_data(
                structure_obj,
                model,
                structures,
                properties,
                _calcsmiles,
                rdkitprop,
                databases,
                timestmp,
            )
            crn_id = "ts" + str(compound_key)
            # for ts the mongodb_id correspond to the rxn id
            html_compounds[compound_key] = {
                **struct_data,
                "crn_id": crn_id,
                "mongodb_id": rxn_id[:-3],
                "pfcost": 0,
            }

        else:  # unimolecular reaction side
            compound, crn_id = _get_compound_and_crnid(
                pathfinder, cmp_dict, compound_id, compounds, flasks
            )
            structure = compound.get_centroid()
            structure_obj = db.Structure(structure, structures)
            struct_data = _extract_structure_data(
                structure_obj,
                model,
                structures,
                properties,
                calcsmiles,
                rdkitprop,
                databases,
                timestmp,
            )
            tmpcost = _get_rescaled_pfcost(pathfinder, compound_id)
            html_compounds[compound_key] = {
                **struct_data,
                "crn_id": crn_id,
                "mongodb_id": compound_id,
                "pfcost": tmpcost,
            }
        # Sort keys in alphabetical order
        tmpdict = html_compounds[compound_key].copy()
        html_compounds[compound_key] = _sort_dict_keys(tmpdict)
    return html_compounds

def _get_rescaled_pfcost(pathfinder, compound_id, threshold=5):
    """
    Return the pathfinder cost value normalized. Because Pathfinder uses a 
    sentinel value of 1e+13. (TO-DO)
    """
    #tmpcost = round(np.log10(pathfinder.compound_costs[compound_id]), 0)
    tmpcost = round(np.log10(pathfinder.compound_costs[compound_id]), 4)
    if tmpcost > threshold:
        tmpcost = threshold
    return tmpcost

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
        A new dictionary with the same key-value pairs as `d` but with keys
        sorted alphabetically.

    Notes
    -----
    - Only the top-level keys are sorted; nested dictionaries are not affected.
    - The ordering is determined by Python's default string comparison.
    """
    return {k: d[k] for k in sorted(d)}


def get_reaction_mechanism_from_A_to_B(source, target, model, pathfinder, 
        manager, npaths=15):
    """
    Determines the most likely reaction mechanisms between a selected source and 
    target compounds. Selection of the most likely mechanisms is done by minimizing
    the overall compound cost.[1]

    [1] Paul L. Tuertscher and Markus Reiher, J. Chem. Inf. Model. 2023, 63, 1, 147-160
    """
    reactions = manager.get_collection("reactions")
    elemsteps = manager.get_collection("elementary_steps")
    structures = manager.get_collection("structures")
    properties = manager.get_collection("properties")
    unique_paths = pathfinder.find_unique_paths(source, target, npaths)
    reaction_dict = {}
    for count, path in enumerate(unique_paths):
        labels, energies, mongodb_ids = [], [], []
        count_ts, count_int = 1, 1
        current_energy = 0
        labels.append("Reactants")
        energies.append(0)
        el_path_str = pathfinder.get_elementary_step_sequence(path[0])
        overall_rxn = pathfinder.get_overall_reaction_equation(path[0])
        cost = round(path[1],1)
        for node in [node for node in path[0]]:
            if ";" in node:  # reaction node
                reaction = db.Reaction(db.ID(node[:-3]))  # the final part '0;0' is deleted
                reaction.link(reactions)
                tmppf = pathfinder.graph_handler.graph.nodes(data=True)
                es_id = db.ID(tmppf[node]["elementary_step_id"])
                side = int(node[-2])
                es_from_graph = db.ElementaryStep(es_id, elemsteps)
                tmpobj = get_energy_and_barriers('electronic_energy', es_id, elemsteps, 
                                                 model, structures, properties, es_from_graph)
                _energy, barrier, not_None = tmpobj
                if not_None:
                    energy = _energy * utils.KJPERMOL_PER_HARTREE
                    elemstepi = db.ElementaryStep(es_id, elemsteps)
                    if elemstepi.get_type() != db.ElementaryStepType.BARRIERLESS:
                        labels.append("TS" + str(count_ts))
                        energies.append(current_energy + barrier[side])
                        count_ts += 1
                        mongodb_ids.append(reaction.get_id().string())  # add rxn id as ts id
                    current_energy += energy
                    labels.append("Int" + str(count_int))
                    energies.append(current_energy)
                    count_int += 1
            else:  # compound node
                mongodb_ids.append(node)
        assert len(labels) == len(energies) == len(mongodb_ids)
        reaction_dict[count] = {"labels":labels, "energies":energies, 
                                "mongodb_ids": mongodb_ids, "cost": cost}
        return reaction_dict

# Static data for Diego. Here I paste the ''reaction_dict'' because otherwise you would need to connect to the MongoDB. 

deleteme = {0: {'labels': ['Reactants', 'TS1', 'Int1', 'TS2', 'Int2', 'Int3'], 'energies': [0, 108.20642677229661, -272.3224795944654, 44.248300059354676, -428.8909335908028, -570.1682927482391], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c99486b7e549636e443a', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 215.2}, 1: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'Int5'], 'energies': [0, -20.353868617363986, 50.05478238807537, -90.49531596626933, -45.3660015758912, -267.31091605768137, 49.2598635961387, -423.87937005401875, -565.1567292114551], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4444', '6649c98d86b7e5495e6e4494', '66543e3686b7e558597e3a0a', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 233.9}, 2: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'Int4'], 'energies': [0, -20.353868617363986, 109.75097210882427, -270.781709285327, 45.78907036849307, -427.3501632816643, -568.6275224391006], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e447c', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 238.4}, 3: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'Int5'], 'energies': [0, -20.353868617363986, 49.98921739251415, -90.45694769686204, -11.96522105479562, 71.17160701759713, 387.7423866714172, -85.39684697874023, -226.6742061361766], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e4484', '6649c98d86b7e5495e6e44b6', '664cc08f86b7e53d14486d42', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 246.3}, 4: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 44.44305332124469, -158.24593193712715, -36.29150621617681, -233.89235168332135, -188.7630372929432, -410.7079517747334, -94.1371721209133, -567.2764057710707, -708.553764928507], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4430', '6649c98c86b7e5495e6e442a', '6649c99586b7e549636e4482', '6649c98d86b7e5495e6e4494', '66543e3686b7e558597e3a0a', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 279.7}, 5: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'Int5'], 'energies': [0, -20.353868617363986, 44.44305332124469, -158.24593193712715, 28.116043981082953, -56.25961189939282, 260.31116775442723, -212.82806589573016, -354.1054250531665], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4430', '6649c98c86b7e5495e6e442a', '6649c99586b7e549636e4474', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 284.1}, 6: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 49.98921739251415, -90.45694769686204, 46.82430176202928, -90.40638952944806, -45.27707513906993, -267.2219896208601, 49.34879003295998, -423.7904436171974, -565.0678027746337], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e4484', '6649c98d86b7e5495e6e44b6', '668074c686b7e53f8c72cbce', '6649c98d86b7e5495e6e4494', '66543e3686b7e558597e3a0a', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 286.9}, 7: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'Int5'], 'energies': [0, -20.353868617363986, 44.44305332124469, -158.24593193712715, 171.62168342773307, -107.07828112029341, 67.52087918142226, 0.05550731322706781, -141.2218518442093], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4430', '6649c98c86b7e5495e6e442a', '66546d8486b7e558597e3a3e', '6649c98c86b7e5495e6e4436', '6649c99486b7e549636e444c', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 289.9}, 8: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 44.44305332124469, -158.24593193712715, -37.360769506831204, -96.69200601392438, -18.20027937185796, 64.93654870053479, 381.50732835435485, -91.63190529580257, -232.90926445323893], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4430', '6649c98c86b7e5495e6e442a', '665442b586b7e558597e3a3c', '6649c98d86b7e5495e6e44b6', '664cc08f86b7e53d14486d42', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 290.7}, 9: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 50.05478238807537, -90.49531596626933, 46.836491660035975, -90.44475779885535, -11.953031156788924, 71.18379691560382, 387.7545765694239, -85.38465708073353, -226.66201623816988], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4444', '6649c98d86b7e5495e6e4494', '668074c686b7e53f8c72cbce', '6649c98d86b7e5495e6e44b6', '664cc08f86b7e53d14486d42', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 298.9}, 10: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 50.05478238807537, -90.49531596626933, -44.18730999151321, -166.14173571246354, 20.220240205746563, -64.15541567472921, 252.41536397909084, -220.72386967106655, -362.0012288285029], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4444', '6649c98d86b7e5495e6e4494', '6649c99586b7e549636e4482', '6649c98c86b7e5495e6e442a', '6649c99586b7e549636e4474', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 307.0}, 11: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'TS5', 'Int6', 'Int7'], 'energies': [0, -20.353868617363986, 49.98921739251415, -90.45694769686204, -31.125711189768865, -28.903021773659276, 93.05140394729106, -104.54944151985349, -59.420127129475354, -281.36504161126555, 35.205738042554515, -437.93349560760294, -579.2108547650394], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e4484', '6649c98d86b7e5495e6e44b6', '665442b586b7e558597e3a3c', '6649c98c86b7e5495e6e442a', '6649c99586b7e549636e4482', '6649c98d86b7e5495e6e4494', '66543e3686b7e558597e3a0a', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 307.3}, 12: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 49.98921739251415, -90.45694769686204, -31.125711189768865, -28.903021773659276, 157.45895414455083, 73.08329826407505, 389.6540779178951, -83.4851557322623, -224.76251488969865], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e4484', '6649c98d86b7e5495e6e44b6', '665442b586b7e558597e3a3c', '6649c98c86b7e5495e6e442a', '6649c99586b7e549636e4474', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 311.7}, 13: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 50.05478238807537, -90.49531596626933, -44.18730999151321, -166.14173571246354, 163.72587965239669, -114.9740848956298, 59.62507540608587, -7.840296462109322, -149.1176556195457], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4444', '6649c98d86b7e5495e6e4494', '6649c99586b7e549636e4482', '6649c98c86b7e5495e6e442a', '66546d8486b7e558597e3a3e', '6649c98c86b7e5495e6e4436', '6649c99486b7e549636e444c', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 312.7}, 14: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'TS5', 'Int6', 'Int7'], 'energies': [0, -20.353868617363986, 50.05478238807537, -90.49531596626933, -44.18730999151321, -166.14173571246354, -45.256573282167594, -104.58780978926077, -26.09608314719435, 57.040744925198396, 373.61152457901846, -99.52770907113896, -240.80506822857532], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4444', '6649c98d86b7e5495e6e4494', '6649c99586b7e549636e4482', '6649c98c86b7e5495e6e442a', '665442b586b7e558597e3a3c', '6649c98d86b7e5495e6e44b6', '664cc08f86b7e53d14486d42', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 313.5}, 15: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 49.98921739251415, -90.45694769686204, -31.125711189768865, -28.903021773659276, 300.9645935912009, 22.264629043174473, 196.86378934489014, 129.39841747669496, -11.878941680741406], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e4484', '6649c98d86b7e5495e6e44b6', '665442b586b7e558597e3a3c', '6649c98c86b7e5495e6e442a', '66546d8486b7e558597e3a3e', '6649c98c86b7e5495e6e4436', '6649c99486b7e549636e444c', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 317.4}, 16: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 50.05478238807537, -90.49531596626933, 123.40182593526104, 62.5434455669844, 279.8843904074332, 24.064478284809276, 198.66363858652494, 131.19826671832976, -10.079092439106603], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4444', '6649c98d86b7e5495e6e4494', '6649c99586b7e549636e4478', '6649c98d86b7e5495e6e44bc', '664cc08e86b7e53d14486d22', '6649c98c86b7e5495e6e4436', '6649c99486b7e549636e444c', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 328.6}, 17: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'TS5', 'Int6', 'Int7'], 'energies': [0, -20.353868617363986, 44.44305332124469, -158.24593193712715, -37.360769506831204, -96.69200601392438, 40.58924344496694, -96.6414478465104, -51.51213345613227, -273.4570479379224, 43.11373171589764, -430.02550193425975, -571.3028610916961], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4430', '6649c98c86b7e5495e6e442a', '665442b586b7e558597e3a3c', '6649c98d86b7e5495e6e44b6', '668074c686b7e53f8c72cbce', '6649c98d86b7e5495e6e4494', '66543e3686b7e558597e3a0a', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 331.3}, 18: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 49.98921739251415, -90.45694769686204, 35.587012582995484, 49.60748049904416, 215.30562167755218, 6.127405317512334, 322.6981849713324, -150.441048678825, -291.71840783626135], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e4484', '6649c98d86b7e5495e6e44b6', '6649c99586b7e549636e4466', '6649c98c86b7e5495e6e4442', '6649c99586b7e549636e4460', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 334.4}, 19: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'TS5', 'Int6', 'Int7'], 'energies': [0, -20.353868617363986, 44.44305332124469, -158.24593193712715, -36.29150621617681, -233.89235168332135, -96.56054405701605, -233.84179351590737, -155.35006687384094, -72.2132388014482, 244.35754085237187, -228.78169279778555, -370.0590519552219], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4430', '6649c98c86b7e5495e6e442a', '6649c99586b7e549636e4482', '6649c98d86b7e5495e6e4494', '668074c686b7e53f8c72cbce', '6649c98d86b7e5495e6e44b6', '664cc08f86b7e53d14486d42', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 344.8}, 20: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'Int6'], 'energies': [0, -20.353868617363986, 44.44305332124469, -158.24593193712715, 25.432606032515025, -220.76151669445713, -55.063375515949105, -264.24159187598895, 52.32918777783112, -420.8100458723263, -562.0874050297627], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4430', '6649c98c86b7e5495e6e442a', '664cc08e86b7e53d14486d20', '6649c98c86b7e5495e6e4442', '6649c99586b7e549636e4460', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 355.3}, 21: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'TS5', 'Int6', 'Int7'], 'energies': [0, -20.353868617363986, 49.98921739251415, -90.45694769686204, 46.82430176202928, -90.40638952944806, -44.09838355469194, -166.0528092756423, 20.309166642567817, -64.06648923790796, 252.50429041591212, -220.63494323424533, -361.9123023916817], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e4484', '6649c98d86b7e5495e6e44b6', '668074c686b7e53f8c72cbce', '6649c98d86b7e5495e6e4494', '6649c99586b7e549636e4482', '6649c98c86b7e5495e6e442a', '6649c99586b7e549636e4474', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 360.0}, 22: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'TS5', 'Int6', 'Int7'], 'energies': [0, -20.353868617363986, 50.05478238807537, -90.49531596626933, 46.836491660035975, -90.44475779885535, -31.11352129176217, -28.89083187565258, 157.4711440425575, 73.09548816208175, 389.66626781590185, -83.4729658342556, -224.75032499169197], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4444', '6649c98d86b7e5495e6e4494', '668074c686b7e53f8c72cbce', '6649c98d86b7e5495e6e44b6', '665442b586b7e558597e3a3c', '6649c98c86b7e5495e6e442a', '6649c99586b7e549636e4474', '6649c98d86b7e5495e6e4486', '6654a36d86b7e558597e3a42', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 364.3}, 23: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'TS5', 'Int6', 'Int7'], 'energies': [0, -20.353868617363986, 49.98921739251415, -90.45694769686204, 46.82430176202928, -90.40638952944806, -44.09838355469194, -166.0528092756423, 163.81480608921794, -114.88515845880855, 59.71400184290712, -7.751370025288068, -149.02872918272442], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99586b7e549636e4484', '6649c98d86b7e5495e6e44b6', '668074c686b7e53f8c72cbce', '6649c98d86b7e5495e6e4494', '6649c99586b7e549636e4482', '6649c98c86b7e5495e6e442a', '66546d8486b7e558597e3a3e', '6649c98c86b7e5495e6e4436', '6649c99486b7e549636e444c', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 365.8}, 24: {'labels': ['Reactants', 'Int1', 'TS1', 'Int2', 'TS2', 'Int3', 'TS3', 'Int4', 'TS4', 'Int5', 'TS5', 'Int6', 'Int7'], 'energies': [0, -20.353868617363986, 50.05478238807537, -90.49531596626933, 46.836491660035975, -90.44475779885535, -31.11352129176217, -28.89083187565258, 300.97678348920766, 22.27681894118117, 196.87597924289685, 129.41060737470164, -11.866751782734724], 'mongodb_ids': ['660d49de86b7e52a0c3bc06c', '6649c98c86b7e5495e6e4424', '6649c99486b7e549636e4444', '6649c98d86b7e5495e6e4494', '668074c686b7e53f8c72cbce', '6649c98d86b7e5495e6e44b6', '665442b586b7e558597e3a3c', '6649c98c86b7e5495e6e442a', '66546d8486b7e558597e3a3e', '6649c98c86b7e5495e6e4436', '6649c99486b7e549636e444c', '6649c98c86b7e5495e6e4438', '6649c98d86b7e5495e6e449a'], 'cost': 370.1}}

