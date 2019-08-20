Installation
============

General Advice
--------------

We strongly recommend that users install 
`conda <https://conda.io/en/latest/miniconda.html>`_ first.
Then, using conda, one can create a new environment and install RiboPy
inside that environment.
For example:

.. code:: bash

      conda create -n ribo python=3.7
      conda activate ribo

Also, we highly encourage using conda for the installation of 
scientific Python packages: numpy, pandas and h5py.   

.. code:: bash
      
      conda install numpy pandas h5py

.. Warning::
   This section needs update!



Requirements
------------

RiboPy is a python package and requires Python 3.6 or a later version.
The python packages, that RiboPy requires, are automatically installed
by pip. The command line interface (CLI) requires terminal access.

If you are going to use bam files, you need bedtools so that
RiboPy can convert bam files to bed entries internally.


Using pip
---------

.. code:: bash

    pip install ribopy

From Github
-----------

Alternatively, you can install the latest version from github.

.. code:: bash

      pip install git+https://github.com/ribosomeprofiling/ribopy.git
