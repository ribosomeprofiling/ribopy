.. _python-api:

Python API
==========


.. toctree::
   :maxdepth: 1


Quick reference
---------------

Ribo Attributes
~~~~~~~~~~~~~~~

Essential Ribo class attributes

.. autosummary:: 
    ribopy.ribo
    ribopy.ribo.Ribo
    ribopy.ribo.Ribo.experiments
    ribopy.ribo.Ribo.minimum_length
    ribopy.ribo.Ribo.maximum_length
    ribopy.ribo.Ribo.metagene_radius
    ribopy.ribo.Ribo.left_span
    ribopy.ribo.Ribo.right_span
    ribopy.ribo.Ribo.format_version
    
    
Getter Functions
~~~~~~~~~~~~~~~~
Methods for reading ribosome profiling data.
    
.. autosummary::    
    ribopy.ribo.Ribo.get_metagene
    ribopy.ribo.Ribo.get_region_counts
    ribopy.ribo.Ribo.get_coverage
    ribopy.ribo.Ribo.get_rnaseq
    
    
Plot Functions
~~~~~~~~~~~~~~
Some essential plots for ribosome profiling analysis,
    
.. autosummary::
    ribopy.ribo.Ribo.plot_metagene
    ribopy.ribo.Ribo.plot_lengthdist
    ribopy.ribo.Ribo.plot_region_counts

    
Ribo
----

.. autoclass:: ribopy.ribo.Ribo
    :members:
 
    
