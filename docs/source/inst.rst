Installation
============

Here we provide a short installation guide for VizChemoton based on using an Anaconda environment (e.g., ``my_env``).

.. code-block:: bash

   # Initialize conda
   source $(conda info --base)/etc/profile.d/conda.sh

   # Create the conda environment with Python 3.6
   conda create --name my_env python=3.6

   # Activate the conda environment
   conda activate my_env

   # Install Chemoton
   git clone https://github.com/qcscine/chemoton.git
   cd chemoton
   git checkout 3.1.0
   # make sure python3 from my_env is called
   python3 -m pip install -r requirements.txt
   python3 -m pip install .
   cd ..

   # Install amk-tools
   git clone https://gitlab.com/dgarayr/amk_tools.git
   cd amk_tools
   python3 -m pip install -e .
   cd ..

   # Install VizChemoton
   git clone https://github.com/petrusen/vizchemoton.git
   cd vizchemoton
   python3 -m pip install -r ./requirements.txt
   python3 -m pip install .

