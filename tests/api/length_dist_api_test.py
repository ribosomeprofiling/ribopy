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

class TestGetLengthDist(ApiTestBase):
    def test_get_length_dist(self):
        length_dist = self.sample_ribo.get_length_dist(
                               region_name = "CDS",
                               experiments = ["merzifon", "ankara"])
        self.assertTrue(np.allclose(length_dist["merzifon"], [3,6,11,6]))
        
        self.assertTrue(np.allclose(length_dist["ankara"], [0,8,0,0]))
        
if __name__ == '__main__':
        
    unittest.main()
