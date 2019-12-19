Installation
============

General Advice
--------------

We strongly recommend that users install
`conda <https://conda.io/en/latest/miniconda.html>`_ first.
Then, using conda, one can install RiboPy in a separate conda environment.

For example, the following command will install RiboPy inside a conda environment named *ribo*.

.. code:: bash

      conda env create -f environment.yaml

The environment file above, `environment.yaml`, can be obtained from the
`RiboFlow repository <https://github.com/ribosomeprofiling/riboflow/blob/master/environment.yaml>`_.


Also, we highly encourage using conda for the installation of
scientific Python packages: numpy, pandas and h5py.

Requirements
------------

RiboPy is a python package and requires Python 3.6 or a later version.
The python packages, that RiboPy requires, are automatically installed
by pip. The command line interface (CLI) requires terminal access.

Using pip
---------

.. code:: bash

    pip install ribopy

From Github
-----------

Alternatively, you can install the latest version from github.

.. code:: bash

      pip install git+https://github.com/ribosomeprofiling/ribopy.git
