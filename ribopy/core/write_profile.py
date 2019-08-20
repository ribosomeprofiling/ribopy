# -*- coding: utf-8 -*-
from collections import OrderedDict

import h5py
import pandas as pd
import numpy as np


def write_profile( profile_array, lengths_array, hdf5_handle, metadata ):

    start_site_coverage = merge_metagene() 
    stop_site_coverage = merge_metagene()

    # Do this in array - pyothonic way
    #UTR3_counts = merge_region_counts()

     
