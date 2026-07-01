Changelog
=========


Release 2.0.0
-------------

Second release of VizChemoton. This version expands and improves the functionalities of the three main output formats (JSON, CSV, and HTML), enabling the transition from prototype development to routine applications.

**Technical Details**

As a result of the expanded functionality, a considerable refactor of the codebase has been carried out to improve modularity and maintainability.

* Refactored `vizchemoton_module` into four dedicated modules: `scine_module`, `cheminfo_module`, `html_module` and `text_module`.
* Added RDKit as a required dependency for handling cheminformatics functionality.
* Approved compatibility with Python 3.8 and 3.10.
* Added unitary tests covering the main functionalities of all four modules.
* Extended CSV outputs to include MongoDB IDs for elementary steps and reactions, facilitating traceability between VizChemoton outputs (HTML, JSON, and CSV) and the original MongoDB database.

**New Features**

* Added SMILES search functionality to the HTML dashboard. Canonical SMILES are computed with RDKit and stored in the `compounds.json` dictionary. SMILES generation can be enabled or disabled through the configuration file.
* Added a static HTML documentation folder.
* Added atom filtering during reaction collection iteration to improve interpretability for large CRNs.
* Added a toggle button in the HTML dashboard to hide barrierless reactions, simplifying the analysis of large and complex CRNs.
* Added a custom JSON dump function to improve the readability of `compounds.json`.
* Added a custom graph layout function for HTML visualization. This replaces the `kamada_kawai` layout from NetworkX, which becomes memory-intensive for large CRNs. The new implementation uses Cartesian geometric descriptors combined with K-Means clustering to determine node positions before graph rendering.
* Added the `node_size` keyword argument to `build_dashboard()` to control node size in the HTML visualization.
* Added node-coloring functionality in the HTML dashboard based on quantitative properties (e.g., energies) or qualitative annotations (e.g., PubChem matches).
* Added queries to the Python APIs of three public chemical structure databases: PubChem, ChEMBL and ChEBI.
* Added InChIKeys to the JSON and HTML outputs, enabling searches for specific target compounds.
* Added an `Export current nodes` button to generate subsets of the reaction network and improve interpretability.

Release 1.0.0
-------------

First implementation of VizChemoton. The package is a light-weight alternative to Heron to access and visualize the compounds, reactions, and transition states of a reaction network, without requiring any hard-to-installl dependencies nor any particular operating system.

**Technical Details**

VizChemoton is implemented in Python and JavaScript, and it can be installed via PIP. To generate the HTML file of the reaction network, it requires the following SCINE modules to extract the chemical data from the exploration:

* SCINE Chemoton, which is necessary to create a Pathfinder object for each exploration, the SCINE database wrapper, which is necessary to query the MongoDB to obtain reaction and energy data.
* SCINE Utilities, which a library of common functionality used across all SCINE modules.

Apart from the SCINE dependencies, VizChemoton also requires the amk-tools package to convert the datafrom the exploration into a single HTML file. These components do not need to be installed manually; rather, they are automatically installed when setting up VizChemoton.

**Current Features**

* Generate the HTML file from a local (or remote) MongoDB where the exploration data is collected. The user can define the method family (e.g. "cc", "dft"), the specific method, the basis set, and the program.
* Generate the HTML file from a JSON file written by Pathfinder containing a set of elementary steps and reactions. The Pathfinder JSON file can be either read of written using VizChemoton.
* Generate the HTML file from a pair of custom data files: compounds.json and reactions.csv. This allows running VizChemoton without an active MongoDB, as well as storing the reaction data in plain text files.
* Define the name, size and layout of the dashboard displaying the reaction network in the HTML file.
* Search any given compound in the reaction network by its index name or SMILES. The former can be obtained searching the structure in ioChem-BD at the server hosted in the Barcelona Supercomputer Center (more information in the documentation).
