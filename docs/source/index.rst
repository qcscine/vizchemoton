.. vizchemoton master file

VizChemoton
=================================

VizChemoton allows the **visualization of chemical reaction networks (CRNs)** constructed by `Chemoton <https://github.com/qcscine/chemoton>`_
through the generation of standalone HTML files with `amk-tools <https://github.com/dgarayr/amk_tools>`_, allowing the user to easily interact
with the network **just via browser**. Furthermore, the reaction network data is interoperable as two plain text files: a JSON file containing
all compounds and their cheminformatics properties, and a CSV file containing all the elementary steps, reactions and transition states in the
network.

.. note::

   If you use VizChemoton, please cite the original article and the corresponding `Zenodo version <https://zenodo.org/records/14803803>`_:
   
   Enric Petrus, Diego Garay-Ruiz, Thomas Weymuth, Markus Reiher, Thomas B. Hofstetter. VizChemoton 2.0.0, ChemRxiv (2026)  
   DOI: `(not available yet) <https://doi.org>`_

.. image:: ../example_crn_html.png
   :alt: Example Image

*The image above shows which type of visualization (network on the left, structure on
the right) the HTML file permits*

Below we depict the table of contents of this documentation. 

.. toctree::
   :maxdepth: 3

   readme
   inst
   inputfile
   examples
   changelog
