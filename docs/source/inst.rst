Installation
============

Here we provide a short installation guide for VizChemoton.

.. code-block:: bash

   # Initialize a Python environment
   python3.8 -m venv env4vizchemoton
   cd env4vizchemoton
   source bin/activate
   
   # Install Chemoton
   git clone https://github.com/qcscine/chemoton.git
   cd chemoton
   # depending on the Chemoton version that was used for the exploration
   git checkout 3.1.0
   python3 -m pip install -r requirements.txt
   python3 -m pip install .
   cd ..
   
   # Install VizChemoton
   git clone https://github.com/petrusen/vizchemoton.git
   cd vizchemoton
   python3 -m pip install -r ./requirements.txt
   python3 -m pip install .

