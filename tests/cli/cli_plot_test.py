# -*- coding: utf-8 -*-
import unittest
from unittest import mock
from unittest.mock import patch
import os
from io import StringIO
import subprocess

import numpy as np
import h5py

from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.settings import *

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)


from multilength_test_data import *

###########################################

NPROCESS = 4

ribo_output_file = "cli_create.ribo"


class TestCLIBase(unittest.TestCase):

    def run_command(self, command):
      """
      Command must be a list.
      See examples below,
      """

      command_str    = tuple( map( str, command ) )
      process        = subprocess.Popen(command_str, 
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE)
      output, error  = process.communicate()
      output_str     = output.decode()
      error_str      = error.decode()

      return (output_str, error_str) 

    def setUp(self):
        self.tmp_files = list()

        files_and_contents = \
            [ (REF_LEN_FILE,     TRANSCRIPT_LENGTHS),
              (ANNOTATION_FILE,  TRANSCRIPT_ANNOTATION),
              (ALIGNMENT_FILE_1, READ_SET_1),
              (METADATA_FILE_1,  METADATA_EXPERIMENT_STR_1),
              (RIBO_META_FILE_1, RIBO_METADATA_STR_1) ]

        for t_file, content_str in files_and_contents:
            with open(t_file, "w") as output_stream:
                print(content_str, file = output_stream)
            self.tmp_files.append(t_file)
            
        self.out_pdf = "plot_cli_out.pdf"
        self.tmp_files.append(self.out_pdf)

        if os.path.isfile(ribo_output_file):
            os.remove(ribo_output_file)

        create_command = ["ribopy", "create", 
                          "--alignmentfile",  ALIGNMENT_FILE_1,
                          "--name",          "merzifon", 
                          "--reference",     "hg38",
                          "--lengths",        REF_LEN_FILE, 
                          "--annotation",     ANNOTATION_FILE,
                          "--metageneradius", METAGENE_RADIUS, 
                          "--leftspan",       LEFT_SPAN, 
                          "--rightspan",      RIGHT_SPAN,
                          "--lengthmin",      2, 
                          "--lengthmax",      5, 
                          "--expmeta",        METADATA_FILE_1,
                          "--ribometa",       RIBO_META_FILE_1,
                          "-n",               4, 
                          ribo_output_file]

        output, error      = self.run_command(create_command)

        #print("Output:", output)
        self.tmp_files.append( ribo_output_file )
        if output:
            print(output)
        if error:
            print(error)


    def tearDown(self):
        for f in self.tmp_files:
            if os.path.isfile(f):
                os.remove(f)

#######################################################################        




class TestCLIPlot(TestCLIBase):

    """
    def tearDown(self):
        TestCLIBase.tearDown(self)
        if os.path.exists(self.out_pdf):
            os.remove(self.out_pdf)
    """


    def test_plot_metagene_start(self):
        command_pieces = ["ribopy", "plot", "metagene",
                          "--lowerlength", "2", 
                          "--upperlength", "3", 
                          "-o",            self.out_pdf,
                          "--site",        "start", 
                          ribo_output_file, "merzifon"]
        output, error = self.run_command(command_pieces)
        self.assertTrue(os.path.exists(self.out_pdf) )


    def test_plot_lengthdist(self):
        command_pieces = ["ribopy", "plot", "lengthdist", 
                          "-o",        self.out_pdf, 
                          "--title",  "title_here",
                          "--region", "CDS", 
                          ribo_output_file, 
                          "merzifon"]

        output, error = self.run_command(command_pieces)

        self.assertTrue(os.path.exists(self.out_pdf) )

    def test_plot_region_counts(self):
        command_pieces = ["ribopy", "plot", "regioncounts",
                          "-o",       self.out_pdf, 
                          "--title",  "title_here",
                          ribo_output_file, 
                          "merzifon"]
        output, error = self.run_command(command_pieces)
        self.assertTrue(os.path.exists(self.out_pdf) )


if __name__ == '__main__':
    unittest.main()
