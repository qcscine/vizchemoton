Input file
=================================

This manual explains the configuration parameters of the *config.yaml* file. A key decision is
whether to connect to the MongoDB database, where Chemoton's exploration data is stored. By default,
VizChemoton renders the final `.html` file without requiring a MongoDB connection 
(`scine.active`: `False`), as the reactions and compounds data are already stored in the *resources* folder.
However, if you wish to adapt VizChemoton to your own system, you need to configure the MongoDB
connection (`scine`) and specify the computational methodology (`method`) that was used. 


.. code-block:: yaml

    scine:
      active: False
      name: "ozone_tme_lcpbe"
      ip: "localhost"
      port: "27017"
      verbose: True
      method:
        method_family: "dft"
        method: "lc-pbe"
        basis_set: "def2-svp"
        solvation: False
        solvent: False
        program: "orca"
        vfilter: False
      pathfinder:
        mode_graph: "write"
        path_graph: "./vizchemoton/resources/pathfinder_graph_tme_dft.json"
        mode_costs: "write"
        path_costs: "./vizchemoton/resources/pathfinder_costs_tme_dft.json"
        init_costs: {"660d49de86b7e52a0c3bc06c": 1, # ozone
                     "660d4a1086b7e52a0c3bc070": 1,   # tme
                     "6649c98d86b7e5495e6e44ac": 1} # water   
        recu_costs: True
    
    cheminfo:
      verbose: True
      rdkit:
          active: False
          method: "hybrid"
          props: ["MolLogP", "CalcTPSA", "MolWt"]
      database:
          pubchem: False
          chembl: False
          chebi: False
    
    files:
      mode_reactions: "write"
      path_reactions: "./vizchemoton/resources/reactions_tme_dft.csv"
      mode_compounds: "write"
      path_compounds: "./vizchemoton/resources/compounds_tme_dft.json"
      path_graphml: "./vizchemoton/resources/network_tme_dft.graphml"

    html:
      pathsearch:
        use_costs: True
        Npaths: 5
        max_length: 6
        source: "c1"
        target: "c48"
      dist_adduct: 3.0
      size:  [1400, 800]
      layout: "kamada_kawai"
      map_field: "pfcost"
      node_size: 25
      title: "Chemoton graph"
      add_editor: True
      path: "network_tme_dft.html"

.. warning::

   Word of caution: the use of database APIs can lead to issues if the APIs cannot be reached, or the URLs have changed.

----

Here we break down all the parameters that can be defined in the input file:

1. scine 
--------

- **active** (``bool``): Enables (``True``) or disables (``False``) the use of the MongoDB. If it is
  disabled, then this section is omitted.  
- **name** (``str``): Name of the MongoDB.
- **ip** (``str``): IP address of the MongoDB server.
- **port** (``str``): Port number for MongoDB communication.
- **method.method_family** (``str``): Specifies the family of the computational method (e.g., ``dft``, ``gfn2``).
- **method.method** (``str``): Name of the computational method (e.g., ``lc-pbe``, ``gfn2``).
- **method.basis_set** (``str``): Basis set used in the calculation (e.g., ``def2-svp``).
- **method.solvent** (``str``): Solvent used in the calculation (e.g., ``water``).
- **method.solvation** (``str``): Solvation model used in the calculation (e.g., ``CPCM``).  
- **method.program** (``str``): Quantum chemistry program used (e.g., ``orca``, ``xtb``).
- **method.vfilter** (``str``): Filter out compounds with an upper atom threshold (e.g., ``{"C": 10, "O": 6}``).
- **pathfinder.mode_graph** (``str``): Either read a preexisting pathfinder graph file (``read``), or write a new one (``write``). Even if
  it is set to ``read``, it will be necessary to have an active connection to the MongoDB.
- **pathfinder.path_graph** (``str``): Path to the pathfinder graph file.
- **pathfinder.mode_costs** (``str``): Either read a preexisting pathfinder costs file (``read``), or write a new one (``write``). Even if
  it is set to ``read``, it will be necessary to have an active connection to the MongoDB. 
- **pathfinder.path_costs** (``str``): Path to the pathfinder costs file.
- **pathfinder.init_costs** (``dict``): Dictionary containing the compound or flasks IDs as keys (``str``) and the custom initial concentrations
as values (``float``).
- **pathfinder.recu_costs** (``bool``): Updates recursively the compound costs (``True``) or it only updates the compound costs once (``False``). Beware that for large reaction networks, updating recursively the compound costs might become a computational bottleneck.
- **verbose** (``bool``): Print additional logs during running time.

2. cheminfo
-----------
- **rdkit.active** (``bool``): Enables (``True``) or disables (``False``) the calculation of RDKit mol objects and
  their cheminformatic properties.
- **rdkit.active** (``str``): Method of choice to calculate RDKit mol objects (``xyz2mol``, ``scinebos`` and ``hybrid``).
- **rdkit.props** (``list``): Cheminformatic properties to be calculated. Names **must** match the native function names in RDKit.   
- **database.pubchem** (``bool``): Query PubChem API to determine which compounds are reported in this database.
- **database.chembl** (``bool``): Query ChEMBL API to determine which compounds are reported in this database.
- **database.chebi** (``bool``): Query ChEBI API to determine which compounds are reported in this database.
- **verbose** (``bool``): Print additional logs during running time.

3. files
--------
- **mode_reactions** (``str``): Either read a preexisting file (``read``) or write a new file (``write``). If one
  sets it to ``read``, because there is a preexisting file, it is not necessary to have an active connection
  to the MongoDB.
- **path_reactions** (``str``): Path to the CSV file containing reaction data.
- **mode_compounds** (``str``): Either read a preexisting file (``read``) or write a new file (``write``). If one
  sets it to ``read``, because there is a preexisting file, it is not necessary to have an active connection
  to the MongoDB.
- **path_compounds** (``str``): Path to the JSON file containing compound data.
- **path_graphml** (``str``): Path to create a GRAPHML file containing the reaction network data (optional).

4. html
-------
- **pathsearch.use_costs** (``bool``): Whether to filter out the most likely reaction pathways based on the previously calculated compound costs (``True``) or not (``False``).
- **pathsearch.Npaths** (``int``)`: Number of paths to be calculated between the source and target compounds.
- **pathsearch.max_length** (``int``)`: Maximum length of the reaction pathways between the source and target compounds.   
- **pathsearch.source** (``str``)`: Selected source compound based on its crn_id.
- **pathsearch.target** (``str``)`: Selected target compound based on its crn_id.
- **dist_adduct** (``float``): Distance threshold for adduct detection.
- **size** (``list[int, int]``): Graph size in pixels (``[width, height]``).
- **layout** (``str``): Graph layout algorithm (e.g., ``kamada_kawai``, ``random``).
- **map_field** (``str``): Property used for node mapping (e.g., ``energy``, ``pfcost``).
- **path** (``str``): Output file name for the generated network visualization.
- **title** (``str``): Title of the network visualization.
- **node_size** (``int``): Size of the nodes in the HTML file.
- **add_Editor** (``bool``): Add a predefined Kekule.js and RDKit.js canvas to the HTML. This allows the user to obtain the 
  SMILES and InChIKey from a custom 2D structure. 

