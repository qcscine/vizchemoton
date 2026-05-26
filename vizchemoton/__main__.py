"""
Enric Petrus, December 2024. Added SCINE helper functions to link with the
amk-tools generation of html files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link
amk-tools and grrm-tools, generating interactive HTML dashboards to
visualize GRRM-generated reaction networks.
"""

import sys
import networkx as nx
from .text_module import (
    vizchemoton_header,
    write_compound_reactions_files,
    read_compound_reactions_files,
    load_config,
    review_compound_file,
    upgrade_compound_file,
    read_filter_file,
)
from .scine_module import get_crn_as_pathfinder, get_reactions_and_compounds
from .html_module import process_graph, build_dashboard, aggregate_property
from .cheminfo_module import db_node_check, compute_cheminf_props


def main():

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

    # Parameters from config
    db_active = config["scine"]["active"]
    db_name = config["scine"]["name"]
    ip = config["scine"]["ip"]
    port = config["scine"]["port"]
    dict_method = config["scine"]["method"]
    verbose = config["scine"]["verbose"]

    pathfinder_file = config["scine"]["pathfinder"]["path"]
    pathfinder_mode = config["scine"]["pathfinder"]["mode"]
    reactions_file = config["files"]["reactions"]["path"]
    reactions_mode = config["files"]["reactions"]["mode"]
    compounds_file = config["files"]["compounds"]["path"]
    compounds_mode = config["files"]["compounds"]["mode"]

    rdkitobj = config["cheminfo"]["rdkit"]["active"]
    smiles_method = config["cheminfo"]["rdkit"]["method"]
    rdkitprop = config["cheminfo"]["rdkit"]["props"]
    pubchem = config["cheminfo"]["database"]["pubchem"]
    chembl = config["cheminfo"]["database"]["chembl"]
    chebi = config["cheminfo"]["database"]["chebi"]
    databases = {
        "pubchem": pubchem,
        "chembl": chembl,
        "chebi": chebi,
    }

    dist_adduct = config["html"]["dist_adduct"]
    size = tuple(config["html"]["size"])
    layout_function = config["html"]["layout"]
    map_field = config["html"]["map_field"]
    node_size = float(config["html"]["node_size"])
    title_html = config["html"]["title"]
    output_file = config["html"]["path"]
    filter_file = config["html"].get("filter_file", None)

    # qualitative or quantitative palette selection
    if "Rank" in map_field:
        palette = ["#d01414", "#f0ce0e", "#12ba14"]
        qual_map = dict(zip([0, 1, 2], palette))
    else:
        palette = "Viridis256"
        qual_map = {}

    # start of Vizchemoton
    vizchemoton_header()
    if db_active:  # the Mongo-DB is reachable
        # read the pathfinder object (to speed-up the process)
        reactions, compounds = [], {}
        if pathfinder_mode == "read":
            manager, pathfinder = get_crn_as_pathfinder(
                ip,
                int(port),
                db_name,
                dict_method,
                write_pathfinder=False,
                read_pathfinder=pathfinder_file,
                verbose=verbose,
            )
            reactions, compounds = get_reactions_and_compounds(
                manager,
                pathfinder,
                dict_method,
                (rdkitobj, smiles_method),
                rdkitprop,
                databases=databases,
                verbose=verbose,
            )

        elif pathfinder_mode == "write":  # write the pathfinder object
            manager, pathfinder = get_crn_as_pathfinder(
                ip,
                int(port),
                db_name,
                dict_method,
                write_pathfinder=pathfinder_file,
                read_pathfinder=False,
                verbose=verbose,
            )
            reactions, compounds = get_reactions_and_compounds(
                manager,
                pathfinder,
                dict_method,
                (rdkitobj, smiles_method),
                rdkitprop,
                databases,
                verbose=verbose,
            )

        # write the reactions and compounds
        if reactions_mode == "write" and compounds_mode == "write":
            write_compound_reactions_files(
                reactions,
                compounds,
                reactions_file,
                compounds_file,
                verbose=verbose,
            )

    if compounds_mode == "review":
        reactions, compounds = read_compound_reactions_files(
            reactions_file, compounds_file, verbose=verbose
        )
        compounds = review_compound_file(compounds_file)
    elif compounds_mode == "read":
        reactions, compounds = read_compound_reactions_files(
            reactions_file, compounds_file, verbose=verbose
        )
    elif compounds_mode == "upgrade":
        reactions, compounds = read_compound_reactions_files(
            reactions_file, compounds_file, verbose=verbose
        )
        compounds = upgrade_compound_file(compounds_file, rdkitprop, databases)

    graph = process_graph(reactions, compounds, dist_adduct)
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

    # adapting collision of modifications
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

    build_dashboard(
        graph,
        compounds,
        title_html,
        output_file,
        size=size,
        layout_function=layout_function,
        map_field=map_field_name,
        node_size=node_size,
        **kwargs_dash,
    )


if __name__ == "__main__":
    main()
