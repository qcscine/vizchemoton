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
      method:
        method_family: "dft"
        method: "lc-pbe"
        basis_set: "def2-svp"
        solvation: False
        solvent: False
        program: "orca"
        vfilter: False
      pathfinder:
        path: "crn_tme_lcpbe.json"
        mode: "read"
      verbose: True
    
    cheminfo:
      rdkit:
          active: False
          method: "hybrid"
          props: ["MolLogP", "CalcTPSA", "MolWt"]
      database:
          pubchem: False
          chembl: False
          chebi: False
      verbose: False
    
    files:
      reactions:
        path: "reactions_tme_dft.csv"
        mode: "read"
      compounds:
        path: "compounds_tme_dft.json"
        mode: "read"
    
    html:
      dist_adduct: 3.0
      size:  [1400, 800]
      layout: "kamada_kawai"
      map_field: "energy"
      node_size: 25
      title: "Chemoton graph"
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
- **method.method_family** (``str``): Specifies the family of the computational method (e.g., ``dft``).
- **method.method** (``str``): Name of the computational method (e.g., ``lc-pbe``).
- **method.basis_set** (``str``): Basis set used in the calculation (e.g., ``def2-svp``).
- **method.solvent** (``str``): Solvent used in the calculation (e.g., ``water``).
- **method.solvation** (``str``): Solvation model used in the calculation (e.g., ``CPCM``).  
- **method.program** (``str``): Quantum chemistry program used (e.g., ``orca``).
- **method.vfilter** (``str``): Filter out compounds with an upper atom threshold (e.g., ``{"C": 10, "O": 6}``).
- **pathfinder.path** (``str``): Path to the json file containing CRN data.
- **pathfinder.mode** (``str``): Either read a preexisting file (``read``), or write a new file (``write``). Even if
  it is set to ``read``, it will be necessary to have an active connection to the MongoDB.
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
- **reactions.path** (``str``): Path to the CSV file containing reaction data.
- **reactions.mode** (``str``): Either read a preexisting file (``read``) or write a new file (``write``). If one
  sets it to ``read``, because there is a preexisting file, it is not necessary to have an active connection
  to the MongoDB.
- **compounds.path** (``str``): Path to the JSON file containing compound data.
- **compounds.mode** (``str``): Either read a preexisting file (``read``) or write a new file (``write``). If one
  sets it to ``read``, because there is a preexisting file, it is not necessary to have an active connection
  to the MongoDB.

4. html
-------
- **dist_adduct** (``float``): Distance threshold for adduct detection.
- **size** (``list[int, int]``): Graph size in pixels (``[width, height]``).
- **layout** (``str``): Graph layout algorithm (e.g., ``kamada_kawai``).
- **map_field** (``str``): Property used for node mapping (e.g., ``degree``).
- **path** (``str``): Output file name for the generated network visualization.
- **title** (``str``): Title of the network visualization.
- **node_size** (``int``): Size of the nodes in the HTML file.

