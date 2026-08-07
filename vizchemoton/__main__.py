"""
Enric Petrus, December 2024. Added SCINE helper functions to link with the
amk-tools generation of html files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link
amk-tools and grrm-tools, generating interactive HTML dashboards to
visualize GRRM-generated reaction networks.
"""

# Standard library imports
import sys

# Third-party library imports
import networkx as nx

# Project-specific SCINE imports
from .text_module import (
    vizchemoton_header,
    write_compound_reactions_files,
    read_compound_reactions_files,
    load_config,
    review_compound_file,
    upgrade_compound_file,
    read_filter_file,
)
from .scine_module import (get_crn_as_pathfinder, 
                           get_reactions_and_compounds,
                           get_reaction_mechanism_from_A_to_B)
from .html_module import (process_graph, 
                          build_dashboard, 
                          aggregate_property, 
                          save_graph,
                          generate_paths,
                          get_energy_ref)
from .cheminfo_module import (db_node_check, 
                              compute_cheminf_props)


def fill_keys_config(config_dict,test_keys,arg={}):
    """Checks if subsections are defined in a dictionary (from config file): if not,
    adds a blank element (dict) to them"""
    for k in test_keys:
        if k not in config_dict.keys():
            config_dict[k] = arg 
    return None


def main():
    # Read Command Line Interface arguments 
    if len(sys.argv) == 1:
        # Take default config.yaml
        configfile = "config.yaml"
    elif len(sys.argv) == 2:
        # Take custom config.yaml
        configfile = sys.argv[1]
    else:
        print("Pass a valid YAML configuration file when running VizChemoton")
        print("python -m vizchemoton [config.yaml]")
        raise ValueError
    
    # Load configuration
    config = load_config(configfile)
    # Mongo-DB settings and quantum chemistry model
    scine_conf = config.get("scine",{})
    db_active = scine_conf.get("active",False)
    db_name = scine_conf.get("name","default")
    ip = scine_conf.get("ip","localhost")
    port = scine_conf.get("port","27017")
    method_conf = scine_conf.get("method", {})
    method_family = method_conf.get("method_family", "gfn2")
    method = method_conf.get("method", "gfn2")
    basis_set = method_conf.get("basis_set", "")
    program = method_conf.get("program", "xtb")
    solvent = method_conf.get("solvent", "")
    solvation = method_conf.get("solvation", "")
    vfilter = method_conf.get("vfilter", "")
    dict_method = {
        "method_family": method_family,
        "method": method,
        "basis_set": basis_set,
        "program": program,
        "solvent": solvent,
        "solvation": solvation,
        "vfilter": vfilter,
    }
    verbose = scine_conf.get("verbose",False)
    pathfinder_conf = scine_conf.get("pathfinder",{})
    
    # Pathfinder object properties
    pf_graph_file = pathfinder_conf.get("path_graph",None)
    pf_graph_mode = pathfinder_conf.get("mode_graph","read")
    pf_costs_file = pathfinder_conf.get("path_costs",None)
    pf_costs_mode = pathfinder_conf.get("mode_costs","read")
    pf_costs_init = pathfinder_conf.get("init_costs",{})
    pf_costs_recu = pathfinder_conf.get("recu_costs",True)

    # Cheminformatics properties
    cheminfo = config.get("cheminfo",{})
    fill_keys_config(cheminfo,["rdkit","database"])
    
    rdkitobj = cheminfo["rdkit"].get("active",False)
    smiles_method = cheminfo["rdkit"].get("method","hybrid")
    rdkitprop = cheminfo["rdkit"].get("props",None)
    pubchem = cheminfo["database"].get("pubchem",False)
    chembl = cheminfo["database"].get("chembl",False)
    chebi = cheminfo["database"].get("chebi",False)
    databases = {
        "pubchem": pubchem,
        "chembl": chembl,
        "chebi": chebi,
    }

    # Output text files
    reactions_file = config["files"]["path_reactions"]
    reactions_mode = config["files"]["mode_reactions"]
    compounds_file = config["files"]["path_compounds"]
    compounds_mode = config["files"]["mode_compounds"]
    output_graph_file = config["files"].get("path_graphml",None)

    # Graphical user interface (HTML) generation
    html_info = config.get("html",{})
    dist_adduct = html_info.get("dist_adduct",3.0)
    size = tuple(html_info.get("size",[1400,800]))
    layout_function = html_info.get("layout","random")
    map_field = html_info.get("map_field","energy")
    node_size = float(html_info.get("node_size",25))
    title_html = html_info.get("title","VizChemoton graph")
    output_file = html_info.get("path_network","network_html")
    filter_file = html_info.get("filter_file", None)
    add_editor = html_info.get("addEditor",True)
    ## qualitative or quantitative palette selection
    if "Rank" in map_field:
        palette = ["#d01414", "#f0ce0e", "#12ba14"]
        qual_map = dict(zip([0, 1, 2], palette))
    else:
        palette = html_info.get("color_palette","Viridis256")
        qual_map = {}


    ### Exploration of paths
    path_search = html_info.get("pathsearch",{})
    if path_search:
        source = path_search["source"]
        target = path_search["target"]
        max_length = path_search.get("max_length",6)
        Npaths = path_search.get("Npaths",1)
        use_costs = path_search.get("use_costs",False)


    # Extract reaction network data
    vizchemoton_header()
    if db_active:  # the Mongo-DB is reachable
        # read the pathfinder object (to speed-up the process)
        reactions, compounds = [], {}
        manager, pathfinder, model = get_crn_as_pathfinder(
            ip,
            int(port),
            db_name,
            dict_method,
            pf_graph_new=(pf_graph_mode, pf_graph_file),
            pf_costs_new=(pf_costs_mode, pf_costs_file, pf_costs_init),
            recursive_cost=pf_costs_recu,
            verbose=verbose,
        )
<<<<<<< HEAD
    #    reactions, compounds = get_reactions_and_compounds(
    #        manager,
    #        pathfinder,
    #        model,
    #        (rdkitobj, smiles_method),
    #        rdkitprop,
    #        databases=databases,
    #        verbose=verbose,
    #    )

    # create reactions and compounds
    if reactions_mode == "write" and compounds_mode == "write" and db_active:
=======
    if reactions_mode == "write" and compounds_mode == "write":
>>>>>>> 35d8407dbda28332eaec5dc9229ebcc6ebf301a2
        reactions, compounds = get_reactions_and_compounds(
            manager,
            pathfinder,
            model,
            (rdkitobj, smiles_method),
            rdkitprop,
            databases=databases,
            verbose=verbose,
        )
        write_compound_reactions_files(
<<<<<<< HEAD
            reactions,
            compounds,
            reactions_file,
            compounds_file,
            verbose=verbose,
        )
    elif compounds_mode == "read" and reactions_mode == "read":
        reactions, compounds = read_compound_reactions_files(
            reactions_file, compounds_file, verbose=verbose
        )
=======
               reactions,
               compounds,
               reactions_file,
               compounds_file,
               verbose=verbose,
           )
>>>>>>> 35d8407dbda28332eaec5dc9229ebcc6ebf301a2
    elif compounds_mode == "review":
        reactions, compounds = read_compound_reactions_files(
            reactions_file, compounds_file, verbose=verbose
        )
        compounds = review_compound_file(compounds_file)
<<<<<<< HEAD
=======
    elif compounds_mode == "read" and reactions_mode == "read":
        reactions, compounds = read_compound_reactions_files(
            reactions_file, compounds_file, verbose=verbose
        )
>>>>>>> 35d8407dbda28332eaec5dc9229ebcc6ebf301a2
    elif compounds_mode == "upgrade":
        reactions, compounds = read_compound_reactions_files(
            reactions_file, compounds_file, verbose=verbose
        )
        compounds = upgrade_compound_file(compounds_file, rdkitprop, databases)
    
    # Preparing HTML GUI
    graph = process_graph(reactions, compounds, dist_adduct)
    if output_graph_file and output_graph_file != "None":
        save_graph(graph,output_graph_file)

    kwargs_dash = {
        "custom_hovers": [],
        "palette": palette,
        "qual_mapping": qual_map,
    }
    if "Rank" in map_field:
        db_name = map_field.replace("Rank", "")
        db_node_check(graph, compounds, db_name)
        kwargs_dash["custom_hovers"] += [
            (f"{db_name}Ids", f"@{db_name}InfoStr")
        ]

    # Adapting collision of modifications
    if rdkitprop:
        compute_cheminf_props(graph, rdkitprop)
        property_hovers = [(prop, f"@{prop}Str") for prop in rdkitprop]
        kwargs_dash["custom_hovers"] += property_hovers

    # For color-mapping cheminf properties, we need some preprocessing
    if (map_field not in ["energy", "degree"]) and ("Rank" not in map_field):
        fun = "mean"
        mapping_values, mapping_flags = aggregate_property(graph, map_field)
        map_field_name = map_field + "_" + fun
        field_to_nodes = {
            nd: {map_field_name: val}
            for nd, val in zip(graph.nodes, mapping_values)
        }
        nx.set_node_attributes(graph, field_to_nodes)
    else:
        map_field_name = map_field

    # Enable post-filtering of exported selections
    if filter_file:
        flag, filter_mapping = read_filter_file(filter_file)
        if flag:
            graph_work = graph.copy()
            out_nodes = []
            for nd in graph.nodes:
                if nd not in filter_mapping["nodes"]:
                    out_nodes.append(nd)
            graph_work.remove_nodes_from(out_nodes)
            graph = graph_work

    # Management of path info 
    if path_search:
        path_list = generate_paths(graph,source,target,max_length=max_length,
                                   Npaths_filt=Npaths,check_costs=use_costs)
        #path_list = get_reaction_mechanism_from_A_to_B(source, target, model, pathfinder,
        #manager, compounds, npaths=Npaths)
        graph.graph["pathList"] = path_list
        nd_ref,e_ref = get_energy_ref(graph,path_list)
        kwargs_dash["alt_ref_energy"] = e_ref

    build_dashboard(
        graph,
        compounds,
        title_html,
        output_file,
        size=size,
        layout_function=layout_function,
        map_field=map_field_name,
        node_size=node_size,
        add_editor=add_editor,
        **kwargs_dash)

if __name__ == "__main__":
    main()
