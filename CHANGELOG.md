Changelog
=========


Release 1.0.0
-------------

**New Features**
* Add SMILES search in the HTML file. This is done by adding a new key in the compounds.json dictionary where the canonical SMILES with rdkit are computed. The calculation of SMILES can be switched on (or off) in the config file.
* Add static HTML documentation folder.
* Add atom filter when iterating the reaction collection in order to make the HTML more interpretable for large CRNs. 

Release 0.1.0
-------------

First implementation of VizChemoton. The package is a light-weight alternative to Heron to access and visualize the compounds, reactions, and transition states of a reaction network, without requiring any hard-to-installl dependencies nor any particular operating system.

**Technical Details**

VizChemoton is implemented in Python and JavaScript, and it can be installed via PIP. To generate the HTML file of the reaction network, it requires the following SCINE modules to extract the chemical data from the exploration:

* SCINE Chemoton, which is necessary to create a Pathfinder object for each exploration,
the SCINE database wrapper, which is necessary to query the MongoDB to obtain reaction and energy data,
* SCINE Utilities, which a library of common functionality used across all SCINE modules.
Apart from the SCINE dependencies, VizChemoton also requires the amk-tools package to convert the data from the exploration into a single HTML file.
These components do not need to be installed manually; rather, they are automatically installed when setting up VizChemoton.

**Current Features**

* Generate the HTML file from a local (or remote) MongoDB where the exploration data is collected. The user can define the method family (e.g. "cc", "dft"), the specific method, the basis set, and the program.
* Generate the HTML file from a JSON file written by Pathfinder containing a set of elementary steps and reactions. The Pathfinder JSON file can be either read of written using VizChemoton.
* Generate the HTML file from a pair of custom data files: compounds.json and reactions.csv. This allows running VizChemoton without an active MongoDB, as well as storing the reaction data in plain text files.
* Define the name, size and layout of the dashboard displaying the reaction network in the HTML file.
* Search any given compound in the reaction network by its index name or SMILES. The former can be obtained searching the structure in ioChem-BD at the server hosted in the Barcelona Supercomputer Center (more information in the documentation).
