# -*- coding: utf-8 -*-

import numpy as np
from ..settings import *


def find_site_coverage( coverage, radius, annotation, site_type ):
    """
    Finds the coverage of a given site ( start / stop) for a given radius
    In ribosome-profiling jargon, this is  known as "metagene analysis".

    Parameters
    ----------
    coverage: This is an ordered dict where 
              each keyword is a reference name each value is an array
              whose length is the nucleotide length of the refrence
              This array has the coverage of the given reference 
              in nucleotide level.
              So, for example, coverage["MYC"][10] 
              gives the number of nucleotides mapped to the 10th position
              of the MYC gene.

    radius: Number of nucleotides to the left and right of the site of interest.
            
    annotation: An list of triplets where each triplet is of the form
                (  [UTR3_start, UTR3_END], [CDS_start, CDS_end], [UTR5_Start, UTR5_end]) 

    site_type: Valid choices are "start" and "stop"
               Finds the coverage around the start or stop site accordingly

    Returns
    -------
    2-dimensional numpy array of size (#references) X ( (2*radius) + 1 )
    Each row corresponds to one reference and each column corresponds to a
    nucleotide position (relative to the site). The site position will be at the middle column.
    

    Notes
    -----
    Note that by our design, the order of references in the 
    coverage OrderedDict and in annotation must be the same.
   
    """   

    if site_type == "start":
        ref_coordinates = list( map(lambda x: (x[0][0], x[1][0], x[2][1]) , annotation ) )
    elif site_type == "stop":
        ref_coordinates = list( map(lambda x: (x[0][0], x[1][1], x[2][1]) , annotation ) )
    else:
        raise IOError("Site type can be either start or stop."
                      "Provided type is " + site_type)

    array_dimensions = (len(ref_coordinates), 2*radius + 1)
    site_coverage = np.zeros(shape = array_dimensions, dtype = SITE_COVERAGE_DT)

    for ref_index, ref_items in enumerate(list(coverage.items() )) :
        ref_name, ref_coverage = ref_items
        this_lower_coord, this_site_coord, this_upper_coord = \
                   ref_coordinates[ref_index]

        index_range_lower = this_site_coord - radius
        index_range_upper = this_site_coord + radius + 1
        index_range_raw = np.arange(index_range_lower, index_range_upper)

        index_range_logical = np.all( [index_range_raw >= this_lower_coord ,\
                               index_range_raw < this_upper_coord ] ,\
                               axis = 0)
         
        this_coverage = site_coverage[ref_index, :]
        for k in range( 2*radius + 1 ):
            if index_range_logical[k]:
                this_coverage[k] = coverage[ref_name][this_site_coord  - radius + k]

        site_coverage[ref_index, :] = this_coverage
  

    return site_coverage
