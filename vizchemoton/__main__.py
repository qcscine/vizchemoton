'''
Enric Petrus, December 2024. Added SCINE helper functions to link with the
amk-tools generation of html files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link
amk-tools and grrm-tools, generating interactive HTML dashboards to
visualize GRRM-generated reaction networks.
'''

import networkx as nx
#from .vizchemoton_module import (vizchemoton_header, get_crn_as_pathfinder,
#                                 get_reactions_and_compounds,
#                                 write_compound_reactions_files,
#                                 read_compound_reactions_files, process_graph,
#                                 build_dashboard, load_config)

from .text_module import (vizchemoton_header, write_compound_reactions_files,
                          read_compound_reactions_files, load_config)
from .scine_module import (get_crn_as_pathfinder, get_reactions_and_compounds)
from .html_module import (process_graph, build_dashboard)
from .cheminfo_module import (pubchem_node_check)

def main():
    # Load configuration
    config = load_config()

    # Parameters from config
    db_active = config["db"]["active"]
    db_name = config["db"]["name"]
    ip = config["db"]["ip"]
    port = config["db"]["port"]
    dict_method = config["method"]

    pathfinder_file = config["files"]["pathfinder"]["path"]
    pathfinder_mode = config["files"]["pathfinder"]["mode"]

    reactions_file = config["files"]["reactions"]["path"]
    reactions_mode = config["files"]["reactions"]["mode"]

    compounds_file = config["files"]["compounds"]["path"]
    compounds_mode = config["files"]["compounds"]["mode"]

    output_file = config["output"]["file"]
    title_html = config["output"]["title"]
    smiles = config["output"]["smiles"]
    verbose = config["output"]["verbose"]

    dist_adduct = config["graph"]["dist_adduct"]
    size = tuple(config["graph"]["size"])
    layout_function = config['graph']['layout']
    #layout_function = getattr(nx, f"{config['graph']['layout']}_layout")
    map_field = config["graph"]["map_field"]
    node_size = float(config["graph"]["node_size"])

    # Start of Vizchemoton
    vizchemoton_header()
    if db_active:  # the Mongo-DB is reachable
        # read the pathfinder object (to speed-up the process)
        reactions, compounds = [], {}
        if pathfinder_mode == 'read':
            manager, pathfinder = get_crn_as_pathfinder(ip, int(
                port), db_name, dict_method, write_pathfinder=False,
                read_pathfinder=pathfinder_file, verbose=verbose)
            reactions, compounds = get_reactions_and_compounds(
                manager, pathfinder, dict_method, calcsmiles=smiles,
                verbose=verbose)

        elif pathfinder_mode == 'write':  # write the pathfinder object
            manager, pathfinder = get_crn_as_pathfinder(ip, int(
                port), db_name, dict_method, write_pathfinder=pathfinder_file,
                read_pathfinder=False, verbose=verbose)
            reactions, compounds = get_reactions_and_compounds(
                manager, pathfinder, dict_method, calcsmiles=smiles,
                verbose=verbose)

        # write the reactions and compounds
        if reactions_mode == 'write' and compounds_mode == 'write':
            write_compound_reactions_files(
                reactions,
                compounds,
                reactions_file,
                compounds_file,
                verbose=verbose)

    # else: # the Mongo-DB is not reachable, or not necessary as reactions and
    # compounds are stored in separate files
    reactions, compounds = read_compound_reactions_files(
        reactions_file, compounds_file, verbose=verbose)
    graph = process_graph(reactions, compounds, dist_adduct)

    kwargs_dash =  {}
    if map_field == "pubchemRank":
        pubchem_node_check(graph,compounds)
        kwargs_dash["custom_hovers"] = [("pubchemIds","@pubchemInfo")]
    
    build_dashboard(
        graph,
        compounds,
        title_html,
        output_file,
        size=size,
        layout_function=layout_function,
        map_field=map_field,
        node_size=node_size,
        **kwargs_dash)


if __name__ == '__main__':
    main()
