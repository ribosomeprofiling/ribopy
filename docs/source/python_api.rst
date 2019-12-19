.. _python-api:

Python API
==========


.. toctree::
   :maxdepth: 1

   api/api_walkthrough.rst


Quick reference
---------------

Ribo Attributes
~~~~~~~~~~~~~~~

Essential Ribo class attributes

.. autosummary::
    ribopy.Ribo
    ribopy.Ribo.experiments
    ribopy.Ribo.minimum_length
    ribopy.Ribo.maximum_length
    ribopy.Ribo.metagene_radius
    ribopy.Ribo.left_span
    ribopy.Ribo.right_span
    ribopy.Ribo.format_version


Getter Functions
~~~~~~~~~~~~~~~~
Methods for reading ribosome profiling data.

.. autosummary::
    ribopy.Ribo.get_metagene
    ribopy.Ribo.get_region_counts
    ribopy.Ribo.get_coverage
    ribopy.Ribo.get_rnaseq

Plot Functions
~~~~~~~~~~~~~~
Some essential plots for ribosome profiling analysis,

.. autosummary::
    ribopy.Ribo.plot_metagene
    ribopy.Ribo.plot_lengthdist
    ribopy.Ribo.plot_region_counts


Ribo
----

.. autoclass:: ribopy.Ribo
    :members:
