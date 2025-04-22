
API
===

Reaction Templates
------------------

.. attribute:: scine_art.reaction_template.MoFrAtIndices: TypeAlias = Tuple[int, int, int]

      A triplet of integers pointing at a specific atom within a template.
      Templates group atoms into fragments and fragments into molecule, this triplet
      encodes this grouping ``(<mol_idx>, <frag_idx>, <atom_idx>)``.
      The first integer indicating the index of the molecule the atom is in,
      the second index pointing to the fragment within the molecule, the third
      pointing at an atom within this fragment.

.. attribute:: scine_art.reaction_template.MoFrAtPair: TypeAlias = Tuple[MoFrAtIndices, MoFrAtIndices]

      Matches two atoms. Both atoms are encoded based on the molecule-fragment-atom
      grouping.

.. attribute:: scine_art.reaction_template.BondChangesContainer: TypeAlias = SideContainer[AssosDissos[MoFrAtPair]]

      A nested dictionary encoding all bond associations (``'assos'``) and
      dissociations (``'dissos'``) that are part of a reaction template.
      Bond modifications are given per side of the reaction (outer dictionary key).

      Example

         >>> bond_changes: BondChangesContainer = ...
         >>> (mol_idx1, frag_idx1, atom_idx1), (mol_idx2, frag_idx2, atom_idx2) = \
         >>>      bond_changes['lhs']['assos']


.. attribute:: scine_art.reaction_template.SideConversionTree: TypeAlias = List[List[List[MoFrAtIndices]]]

      Maps one atom encoded by a molecule-fragment-atom index set onto another that is
      encoded the same way.

      Example

         >>> mapping: SideConversionTree = ...
         >>> rhs_mol_idx, rhs_frag_idx, rhs_atom_idx = \
         >>>      mapping[lhs_mol_idx][lhs_frag_idx][lhs_atom_idx]

.. attribute:: scine_art.reaction_template.ShapesByMoFrAt: TypeAlias = List[List[List[Optional[masm.shapes.Shape]]]]

      Lists atom shapes for atoms encoded by a molecule-fragment-atom triple

      Example

         >>> shapes: ShapesByMoFrAt = ...
         >>> atom_shape = shapes[mol_idx][frag_idx][atom_idx]


.. attribute:: scine_art.reaction_template.Match: TypeAlias = Dict[str, List[Tuple[Tuple[int, int], Tuple[int, int]]]]

      Indices references the given lists of molecules in the generating call.

      Example

         >>> molecules = [masm_mol1, masm_mol2]
         >>> match: Match = matching_function(molecules, ...)
         >>> for association in match['assos']:
         >>>     mol_idx1, atom_idx1 = association[0]
         >>>     mol_idx2, atom_idx2 = association[1]
         >>>     # connect the specified atoms to build product structures
         >>>     #   (repeat with disconnections for dissociations)

.. automodule:: scine_art.reaction_template
   :members:
   :exclude-members: Match, MoFrAtPair, MoFrAtIndices, ShapesByMoFrAt, SideConversionTree, BondChangesContainer

Reaction Template Database
--------------------------

.. automodule:: scine_art.database
   :members:


Molecule Comparisons
--------------------

.. automodule:: scine_art.molecules
   :members:


IO
--

.. automodule:: scine_art.io
   :members:


Experimental: Atom Matching
---------------------------

.. automodule:: scine_art.experimental
   :members:
