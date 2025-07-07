#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

# Standard library imports
from json import dumps
from typing import List, Tuple
import os
import inspect
import unittest
import pytest

# Local application tests imports
from scine_chemoton.gears import HoldsCollections
from scine_chemoton.engine import Engine
from scine_chemoton.gears.reaction import BasicReactionHousekeeping

# Third party imports
import scine_database as db
import scine_utilities as utils
from scine_database import test_database_setup as db_setup
from scine_chemoton.gears.pathfinder import Pathfinder as pf

# Local imports
from ..vizchemoton_module import get_reactions_and_compounds

class VizChemotonTests(unittest.TestCase, HoldsCollections):

    def custom_setup(self, manager: db.Manager) -> None:
        self._required_collections = ["manager", "elementary_steps", "structures", "reactions", "compounds", "flasks",
                                      "properties"]
        self.initialize_collections(manager)

    def tearDown(self) -> None:
        self._manager.wipe()

    # Capture std out end err
    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_pathfinder_build_graph_and_find_paths(self):
        # prepare settings for creating a generic crn
        n_compounds = 7
        n_reactions = 8
        max_r_per_c = 7
        max_n_products_per_r = 3
        max_n_educts_per_r = 3
        max_s_per_c = 1
        max_el_steps_per_r = 1
        barrier_limits = (10, 80)
        n_inserts = 2
        n_flasks = 1
        # create a generic crn
        manager = db_setup.get_random_db(
            n_compounds,
            n_flasks,
            n_reactions,
            max_r_per_c,
            "test_pathfinder_build_graph",
            max_n_products_per_r,
            max_n_educts_per_r,
            max_s_per_c,
            max_el_steps_per_r,
            barrier_limits,
            n_inserts,
        )
        self.custom_setup(manager)
        # define arbitrary parameters
        dmethod = {"method_family": "FAKE", "method": "FAKE", "basis_set": "F-AKE", "program": "FA-KE"}
        model1 = db.Model(dmethod["method_family"], dmethod["method"], dmethod["basis_set"])
        # construct arbitrary pathfinder object
        pathfinder = pf(manager)
        pathfinder.options.model = model1 
        pathfinder.options.graph_handler = "barrier"
        pathfinder.options.use_structure_model = True
        pathfinder.options.structure_model = model1
        pathfinder.build_graph()
        # test the get_reactions_and_compounds() 
        reactions, compounds = get_reactions_and_compounds(manager, pathfinder, dmethod)
        assert len(reactions) != 0
        assert isinstance(reactions, list)
        assert len(compounds.keys()) != 0
        assert isinstance(compounds, dict)











