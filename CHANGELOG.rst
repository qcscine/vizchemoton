Changelog
=========

Release 1.0.0
-------------

Initial Features
^^^^^^^^^^^^^^^^

- HTML dashboard to display the results of reaction exploration performed with Chemoton:
    - 3D interactive visualization of compounds, flasks, and transition states using JSMol
    - Node information: charge, multiplicity, formula, and tag
    - Edge information: charge, multiplicity, formula, tag, and forward/backward energy barriers
    - Zoom in and out of the reaction network
    - Highlight neighboring nodes and edges

- Integration with online repository: ioChem-BD:
    - Chemical data in the HTML file can be accessed via the ioChem-BD platform
    - Structure-based search supported via the Find module of ioChem-BD

- Generation of plain-text reaction network data from MongoDB:
    - Create `compounds.json` containing chemical information for all nodes
    - Create `reaction.csv` containing all reactions in the network
    - These files allow usage without an active connection to MongoDB

