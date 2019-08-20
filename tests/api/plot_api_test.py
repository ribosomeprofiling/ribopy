# -*- coding: utf-8 -*-
import unittest
import os
from io import StringIO, BytesIO

import h5py

from ribopy import Ribo
from ribopy import create
from ribopy.merge import merge_ribos
from ribopy.settings import *
from ribopy.core.exceptions import *

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)

from multilength_test_data import *
from api_test_base import ApiTestBase

####################################################################
#  NOTE
# For the plot functions, we are currently only checking that
# They don't raise an exception
####################################################################

class TestCreate(ApiTestBase):
 
    ################################################################
    
    def test_plot_metagene(self):
        try: 
            self.sample_ribo.plot_metagene(
                title       = "Start Site Plot",
                site_type   = "start",
                experiments = "merzifon",
                normalize   = True
            )
            
            self.sample_ribo.plot_metagene(
                title       = "Start Site Plot",
                site_type   = "start",
                range_lower = 3,
                range_upper = 4,
                experiments = "merzifon",
                normalize   = True
            )
            
            self.sample_ribo.plot_metagene(
                site_type   = "stop",
                experiments = ["ankara"],
                output_file = BytesIO() 
            )
        except Exception as e:
            self.fail("plot_metagene raised an exception.\n" + str(e))
        
    def test_plot_lengthdist(self):
        try:
            self.sample_ribo.plot_lengthdist( 
                                 region_type = "CDS", 
                                 experiments = self.sample_ribo.experiments,  
                                 normalize   = False,
                                 colors      = PLOT_COLORS )
                                 
            self.sample_ribo.plot_lengthdist( 
                                 region_type = "UTR5", 
                                 experiments = "merzifon",  
                                 normalize   = True,
                                 output_file = "Length Dist")
                                 
            self.sample_ribo.plot_lengthdist( 
                                 region_type = "UTR3_junction", 
                                 experiments = "merzifon",  
                                 normalize   = True,
                                 output_file = BytesIO())
        except Exception as e:
            self.fail("plot_lengthdist raised an exception.\n" + str(e))
            
    def test_plot_region_counts(self):
        try:
            self.sample_ribo.plot_region_counts(
                               experiments      = self.sample_ribo.experiments,
                               title            = "",
                               range_lower      = 2, 
                               range_upper      = 4, 
                               horizontal       = True)
                               
            self.sample_ribo.plot_region_counts(
                               experiments      = "ankara",
                               title            = "",
                               range_lower      = 2, 
                               horizontal       = True,
                               output_file      = BytesIO())
        except Exception as e:
            self.fail("plot_region_counts raised an exception.\n" + str(e))
            

        
if __name__ == '__main__':
        
    unittest.main()
