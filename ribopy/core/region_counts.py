# -*- coding: utf-8 -*-

import numpy as np
from ..settings import *

def get_extended_boundaries(annotation, left_span, right_span):
    extended_annotation = list()

    """
        Note that a is of the form
         ( (UTR5_orig_start, UTR5_orig_end),
           (CDS_orig_start, CDS_orig_end),
           (UTR3_orig_start, UTR3_orig_end)
          )
    """
    for a in annotation:

        UTR5_start = 0
        #UTR5_junction_start = max(a[1][0] - left_span, 0) # old version
        UTR5_junction_start = max(a[0][1] - left_span, 0)
        UTR5_end = UTR5_junction_start
        UTR5_junction_end   = min(a[1][0] + right_span + 1 , a[1][1] )
        CDS_start = UTR5_junction_end
        CDS_end = max( a[1][1] - left_span, CDS_start )
        UTR3_junction_start = CDS_end
        #UTR3_junction_end = min( a[1][1] + right_span + 1 , a[2][1])# old version
        UTR3_junction_end = min( a[1][1] + right_span + 1 , max( a[1][1] , a[2][1] ) )
        UTR3_start = UTR3_junction_end
        UTR3_end   = a[2][1] # this could have been max( a[2][1], UTR3_end )  

        extended_annotation.append( ( (UTR5_start, UTR5_end),
                                      (UTR5_junction_start, UTR5_junction_end),
                                      (CDS_start, CDS_end),
                                      (UTR3_junction_start, UTR3_junction_end),
                                      (UTR3_start, UTR3_end) )
                                  )

    return extended_annotation

def find_region_counts( coverage, annotation, left_span, right_span ):
    """
    Finds the number of counts for each region.

    Each transcript is partitioned into 5 regions:
    UTR5, UTR5_junction, CDS, UTR3_junction and UTR3

    The boundaries of these regions are as follows:
    Note that the second coordinate component is excluded

    UTR5 = [0, start_site - left_span )
    UTR5_junction: [start_site - left_span, start_site + right_span + 1 )
    CDS = [ start_site + right_span + 1, stop_site - left_span )
    UTR3_junction: [stop_site - left_span, stop_site + right_span + 1 )
    UTR3: [ stop_site + right_span  +1, REF_length)

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

    left_span: Number of nucleotides to the left of the start / stop site.
               This is used to determine the junction regions.

    annotation: An list of triplets where each triplet is of the form
                (  [UTR3_star, UTR3_END],
                   [CDS_start, CDS_end],
                   [UTR5_Start, UTR5_end])

    Returns
    -------
    region_counts: a 2-D numpy array of dimension #references X 5
                   each coulmn corresponds to a region
    """

    extended_region_boundaries = \
            get_extended_boundaries(annotation, left_span, right_span)

    # initialize all region counts to 0
    # Note that 5 is for the 5 regions of interest
    # UTR5, UTR5_junction, CDS, UTR3_junction, UTR3
    region_counts = np.zeros(shape = (len(annotation), 5),
                             dtype = REGION_COUNTS_DT)

    # Since the two lists are ordered the same way
    # we can safely zip them
    # Note that coverage is coming from and ordereddict
    iter_items = zip( extended_region_boundaries,
                      coverage.values(),
                      list(range(len(annotation))) )

    for boundaries, cover, ref_index in iter_items:
        for region_ind, start_stop in enumerate(boundaries):
            region_counts[ref_index, region_ind] = \
                np.add.reduce( cover[ start_stop[0]: start_stop[1] ] )

    return region_counts
