Changelog
=========


Release 2.0.0
-------------

Second release of VizChemoton. We expanded and improved the functionalities of the three main output files (JSON, CSV and HTML) created by VizChemoton to move from prototype to routine applications.

**Technical Details**

As a result of the increase in VizChemoton's functionalities, a considerable refactor has been done to accomodate the new functionalities. 

* Refactored "vizchemoton\_module" into four different modules: "scine\_module", "cheminfo\_module", "html\_module" and "text\_module".
* RDKit is now a required dependency to manage cheminformatic results.
* Approved compatibility with Python 3.8
* Added unitary test to check the main functionalities of the four modules
* CSV file contains the MongoDB IDs of the elementary steps and reactions. This facilitates tracing back reaction results from the VizChemoton files (HTML, JSON and CSV) to the original MongoDB.  

**New Features**
* Added SMILES search in the HTML file. This is done by adding a new key in the compounds.json dictionary where the canonical SMILES with rdkit are computed. The calculation of SMILES can be switched on (or off) in the config file.
* Added static HTML documentation folder.
* Added atom filter when iterating the reaction collection in order to make the HTML more interpretable for large CRNs. 
* Added button in the HTML dashboard for hiding barrierless reactions, thus easing the interpretation of large/complex CRNs. 
* Added a custom JSON dump function to improve the readability of the compounds.json file. 
* Added a custom layout function for plotting the graph in HTML format. This bypasses the use of the "kamada\_kawai" option in NetworkX because it consumes too much memory for large CRNs. The custom function uses Cartesian geometric descriptors for each structure in order to then perform a K-Means clustering. The (x,y) position of each structure according to the clustering is then input in the NetworkX plotting function.
* Added the kwargs "node\_size" in the build\_dashboard() to control the size of the nodes in the HTML. 
* Added plotting functionality to color the nodes of the HTML based on quantitative (e.g., energy) or qualitative (e.g., PubChem matches).
* Added query to the Python APIs of three public chemical structural databases: PubChem, ChEMBL and ChEBI. 
* Added InChiKeys in the JSON and HTML so that it is possible to search for specific target compounds.
* Added "Export current nodes" button to create a subset of the reaction network to ease interpretability.

Release 1.0.0
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
