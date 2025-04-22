Input file
=================================

This manual explains the configuration parameters of the *config.yaml* file. A key decision is
whether to connect to the MongoDB database, where Chemoton's exploration data is stored. By default,
VizChemoton renders the final `.html` file without requiring a MongoDB connection 
(`db.active`: `False`), as the reactions and compounds data are stored in the *resources* folder.
However, if you wish to adapt VizChemoton to your own system, you need to configure the MongoDB
connection (`db`) and specify the computational methodology (`method`) that was used.  


.. code-block:: yaml

    db:
      active: False
      name: "ozone_tme_lcpbe"
      ip: "localhost"
      port: "8889"

    method:
      method_family: "dft"
      method: "lc-pbe"
      basis_set: "def2-svp"
      program: "orca"

    files:
      pathfinder:
        path: "./vizchemoton/resources/crn_pathfinder.json"
        mode: "read"
      reactions:
        path: "./vizchemoton/resources/reactions.csv"
        mode: "read"
      compounds:
        path: "./vizchemoton/resources/compounds.json"
        mode: "read"

    output:
      file: "network.html"
      title: "Chemoton graph"
      verbose: True

    graph:
      dist_adduct: 3.0
      size: [1400, 800]
      layout: "kamada_kawai"
      map_field: "degree"


Here we break down all the parameters that can be defined in the input file:

----

1. Database 
--------------------

- **active** (``bool``): Enables (``True``) or disables (``False``) the use of the MongoDB. If it is
  disabled, then this section is omitted.  
- **name** (``str``): Name of the MongoDB.
- **ip** (``str``): IP address of the MongoDB server.
- **port** (``str``): Port number for MongoDB communication.

2. Method
------------------------------------

- **method_family** (``str``): Specifies the family of the computational method (e.g., ``dft``).
- **method** (``str``): Name of the computational method (e.g., ``lc-pbe``).
- **basis_set** (``str``): Basis set used in the calculation (e.g., ``def2-svp``).
- **program** (``str``): Quantum chemistry program used (e.g., ``orca``).

3. Files
-------------------------

pathfinder
^^^^^^^^^^

- **path** (``str``): Path to the json file containing CRN data.
- **mode** (``str``): Either read a preexisting file (``read``), or write a new file (``write``). Even if
  it is set to ``read``, it will be necessary to have an active connection to the MongoDB.

reactions
^^^^^^^^^

- **path** (``str``): Path to the csv file containing reaction data.
- **mode** (``str``): Either read a preexisting file (``read``) or write a new file (``write``). If one
  sets it to ``read``, because there is a preexisting file, it is not necessary to have an active connection
  to the MongoDB.

compounds
^^^^^^^^^

- **path** (``str``): Path to the json file containing compound data.
- **mode** (``str``): Either read a preexisting file (``read``) or write a new file (``write``). If one
  sets it to ``read``, because there is a preexisting file, it is not necessary to have an active connection
  to the MongoDB.

4. Settings 
-----------------------------

- **dist_adduct** (``float``): Distance threshold for adduct detection.
- **size** (``list[int, int]``): Graph size in pixels (``[width, height]``).
- **layout** (``str``): Graph layout algorithm (e.g., ``kamada_kawai``).
- **map_field** (``str``): Property used for node mapping (e.g., ``degree``).

5. Output 
----------------------

- **file** (``str``): Output file name for the generated network visualization.
- **title** (``str``): Title of the network visualization.
- **verbose** (``bool``): Enables (``True``) or disables (``False``) the call to the print statements during
  runtime of the code.


