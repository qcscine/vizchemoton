#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher
Group. See LICENSE.txt for details.
"""

# Standard library imports
import os
import unittest

# Third party imports
from scine_chemoton.gears import HoldsCollections
from scine_chemoton.gears.pathfinder import Pathfinder as pf
import scine_database as db
from scine_database import test_database_setup as db_setup
from bokeh.plotting import Figure

# Local imports
from vizchemoton.tests.resources import resources_root_path
from vizchemoton.vizchemoton_module import (get_reactions_and_compounds,
                                            convert_struct_to_smile,
                                            read_compound_reactions_files,
                                            process_graph,
                                            build_dashboard)


class VizChemotonTests(unittest.TestCase, HoldsCollections):
    """
    Tests the main -and thus most critical- functions of VizChemoton.
    """

    def custom_setup(self, manager: db.Manager) -> None:
        """
        Initializes a custom database.
        """
        self._required_collections = [
           "manager",
           "elementary_steps",
           "structures",
           "reactions",
           "compounds",
           "flasks",
           "properties"]
        self.initialize_collections(manager)

    def test_convert_struct_to_smiles(self):
        """
        Test that convert_struct_to_smiles() behaves properly, converting
        simple cartesian files into SMILES.
        """
        test_molec = ["test_carbondioxide.xyz",
                      "test_h2o2.xyz",
                      "test_hoocohcl.xyz",
                      "test_ozonide.xyz",
                      "test_ozone.xyz",
                      "test_flask.xyz",
                      "test_flask_3h2o.xyz",
                      "test_ts.xyz",
                      "test_h6.xyz"]
        # connect to test DB
        manager = db_setup.get_clean_db("chemoton_test_compound_creation")
        self.custom_setup(manager)
        # add structure data
        rr = resources_root_path()
        manager.init()
        lcentroids = []
        for ipath in test_molec:
            structure = db.Structure()
            structure.link(self._structures)
            structure.create(os.path.join(rr, ipath), 0, 1)
            lcentroids.append(structure)
        dsmiles = {}
        for ipath, icentr in zip(test_molec, lcentroids):
            dsmiles[ipath] = convert_struct_to_smile(icentr)
        # check five typical ozonation products
        print(dsmiles)
        assert dsmiles["test_carbondioxide.xyz"]['smiles'] == 'O=C=O'
        assert dsmiles["test_h2o2.xyz"]['smiles'] == 'OO'
        assert dsmiles["test_hoocohcl.xyz"]['smiles'] == 'OOC(O)Cl'
        assert dsmiles["test_ozonide.xyz"]['smiles'] == 'C1COOO1'
        assert dsmiles["test_ozone.xyz"]['smiles'] == 'O=[O+][O-]'
       # test two flask examples - tricky for smiles generation
        assert dsmiles["test_flask.xyz"]['smiles'] == 'Cl.O.[C-]#[O+]'
        assert dsmiles["test_flask_3h2o.xyz"]['smiles'] == 'O.O.O'
        assert dsmiles["test_ts.xyz"]['smiles'] == 'Cl.[O-][O+]=CO'
        assert dsmiles["test_h6.xyz"]["flag"] == False

    def test_get_reactions_and_compounds(self):
        """
        Tests that get_reactions_and_compounds() correctly reads the reaction
        network data from a pathfinder object.
        """
        # prepare settings for creating a generic crn
        params = {
            "n_compounds": 7,
            "n_flasks": 1,
            "n_reactions": 8,
            "max_r_per_c": 7,
            "name": "test_pathfinder_build_graph",
            "max_n_products_per_r": 3,
            "max_n_educts_per_r": 3,
            "max_s_per_c": 1,
            "barrier_limits": (10, 80),
            "n_inserts": 2,
        }
        # create a generic crn
        manager = db_setup.get_random_db(**params)
        self.custom_setup(manager)
        # define arbitrary parameters
        dmethod = {
            "method_family": "FAKE",
            "method": "FAKE",
            "basis_set": "F-AKE",
            "program": "FA-KE",
            "vfilter": None}
        # pylint: disable=no-member
        model1 = db.Model(
            dmethod["method_family"],
            dmethod["method"],
            dmethod["basis_set"])
        # construct arbitrary pathfinder object
        pathfinder = pf(manager)
        pathfinder.options.model = model1
        pathfinder.options.graph_handler = "barrier"
        pathfinder.options.use_structure_model = True
        pathfinder.options.structure_model = model1
        pathfinder.build_graph()
        # test the get_reactions_and_compounds()
        reactions, compounds = get_reactions_and_compounds(
            manager, pathfinder, dmethod, verbose=True)
        assert len(reactions) != 0
        assert isinstance(reactions, list)
        assert len(compounds.keys()) != 0
        assert isinstance(compounds, dict)

    def test_process_graph_and_build_dashboard(self):
        """
        Tests that the conversion to a NetworkX object is successfuly
        -and consistently- done.
        """
        rr = resources_root_path()
        compounds_file = "test_compounds.json"
        reaction_file = "test_reactions.csv"
        rfile, cfile = os.path.join(
            rr, reaction_file), os.path.join(
            rr, compounds_file)
        reactions, compounds = read_compound_reactions_files(
            rfile, cfile, verbose=True)
        graph = process_graph(reactions, compounds, dist_adduct=3.0)
        assert len(graph.edges) == 24
        assert len(graph.nodes) == 25
        outfile, title = os.path.join(rr, "test_network.html"), 'test_network'
        bokehobj = build_dashboard(graph, title, outfile)
        assert any(isinstance(x, Figure)
                   for x in bokehobj), "No Figure in build_dashboard output"
