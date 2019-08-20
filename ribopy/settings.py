# -*- coding: utf-8 -*-
from collections import OrderedDict

import numpy as np

#####################################################
# C O M M O N    N A M I N G   C O N V E N T I O N S
#####################################################
EXPERIMENTS_name = "experiments"
REFERENCE_name = "reference"

UTR5_name = "UTR5"
UTR3_name = "UTR3"
UTR5_JUNCTION_name = "UTR5_junction"
UTR3_JUNCTION_name = "UTR3_junction"
CDS_name = "CDS"

REGION_names          = [UTR5_name , CDS_name, UTR3_name ]
EXTENDED_REGION_names = [UTR5_name , UTR5_JUNCTION_name , 
                         CDS_name  , UTR3_JUNCTION_name , 
                         UTR3_name ]

METAGENE_name        = "metagene"
METAGENE_RADIUS_name = "metagene_radius"

LEFT_SPAN_name     = "left_span"
RIGHT_SPAN_name    = "right_span"

REF_DG_REFERENCE_NAMES   = "reference_names"
REF_DG_REFERENCE_LENGTHS = "reference_lengths"
REF_DG_START_SITE_COV    = "start_site_coverage"
REF_DG_STOP_SITE_COV     = "stop_site_coverage"
REF_DG_REGION_COUNTS     = "region_counts"
REF_DG_COVERAGE          = "coverage"

REF_ANNOTATION_NAME = "annotation"

LENGTH_MIN_name = "length_min"
LENGTH_MAX_name = "length_max"

RNASEQ_name = "rnaseq"


#####################################################
# H D F 5    S E T T I N G S
#####################################################
DEFAULT_COMPRESSION = "gzip"
DEFAULT_FLETCHER32  = True

SITE_COVERAGE_DT       = np.uint32
TRANSCRIPT_COVERAGE_DT = np.uint16
REGION_COUNTS_DT       = np.uint32
DEFAULT_COUNT_DT       = np.uint32
RNASEQ_DT              = np.float32

#####################################################
# B E D   C O L U M N   N A M E S
#####################################################
BED_CHROM  = "chrom"
BED_START  = "start"
BED_END    = "end"
BED_NAME   = "name"
BED_SCORE  = "score"
BED_STRAND = "strand"

BED_COLUMN_NAMES = [BED_CHROM,
                    BED_START, BED_END,
                    BED_NAME, BED_SCORE,
                    BED_STRAND ]

################################################################
# D A T A F R A M E   C O L U M N   N A M E S & S E T T I N G S
################################################################
DF_CSV_SEPARATOR = ","

DF_TRANSCRIPT      = "transcript"
DF_EXPERIMENT_NAME = "experiment"
DF_READLENGTH      = "read_length"

#####################################################
# M E T A D A T A
#####################################################
# Standard metadata for ribo files

ATTRS_REFERENCE       = "reference"
ATTRS_VERSION         = "ribopy_version"
ATTRS_FORMAT_VERSION  = "format_version"
ATTRS_TIME            = "time"
ATTRS_LENGTH_MIN      = "length_min"
ATTRS_LENGTH_MAX      = "length_max"
ATTRS_METAGENE_RADIUS = "metagene_radius"
ATTRS_LEFT_SPAN       = "left_span"
ATTRS_RIGHT_SPAN      = "right_span"
ATTRS_TOTAL_READS     = "total_reads"


# These are the attributes to be attributes to be copied 
# when merging ribo files.
ATTRIBUTES_TO_BE_COPIED_FOR_MERGE = \
    ( ATTRS_REFERENCE, ATTRS_VERSION, ATTRS_FORMAT_VERSION,
      ATTRS_LENGTH_MIN, ATTRS_LENGTH_MAX, ATTRS_METAGENE_RADIUS,
      ATTRS_LEFT_SPAN, ATTRS_RIGHT_SPAN )

ATTRS_ESSENTIAL_COMPATIBILITY = (ATTRS_REFERENCE, 
                                 ATTRS_FORMAT_VERSION,
                                 ATTRS_LENGTH_MIN, 
                                 ATTRS_LENGTH_MAX,
                                 ATTRS_METAGENE_RADIUS,
                                 ATTRS_LEFT_SPAN, 
                                 ATTRS_RIGHT_SPAN)

# This determines how the attrributes are displayed
# when information of rthe ribo file is requested
RIBO_METADATA_FOR_DISPLAY = OrderedDict(
    [ 
    (ATTRS_FORMAT_VERSION, "Ribo File Version"),
    (ATTRS_VERSION,        "RiboPy Version"),
    (ATTRS_TIME,           "Creation Time"),
    (ATTRS_REFERENCE,      "Reference"),
    (ATTRS_LENGTH_MIN,     "Min Read Length"),
    (ATTRS_LENGTH_MAX,     "Max Read Length"),
    (ATTRS_METAGENE_RADIUS,"Metagene Radius"),
    (ATTRS_LEFT_SPAN,      "Left Span"),
    (ATTRS_RIGHT_SPAN,     "Right Span") 
    ]
)

# The following metadata must match
# when merging experiments into one ribo file
RIBO_ESSENTIAL_METADATA =\
{
    ATTRS_METAGENE_RADIUS:\
        "number of nucleotides to the right "
        "and left of the start / stop sites",
    ATTRS_LEFT_SPAN: \
        "For 5/3 prime junction regions: "
        "number of nucleotides to the left of start / stop site",
    ATTRS_RIGHT_SPAN:\
        "For 5/3 prime junction regions: "
        "number of nucleotides to the right of start / stop site"
}

# USER METADATA is stored under this name
# in the attributes of the ribo file or the experiment
USER_METADATA = "metadata"

# Some essential metadata 
# The values correspond to display names
"""
LIBRARY_ESSENTIAL_METADATA = OrderedDict(\
    [ ("cell_line", "Cell Line"),
      ("sra", "SRA Number"),
      ("link", "External Link"),
      ("digestion_enzyme", "Digestion Enzyme"),
      ("digestion_duration", "Digestion Duration"),
      ("digestion_enzyme_volume", "D. Enz. Volume"),
      ("exp_prep_method", "Library Prep. Method"),
      ("chip", "Chip"),
      ("notes", "Notes:")]
)
"""

# All other metadata will be stored in "metadata" attribute of the
# sample in json format


#####################################################
# P L O T 
#####################################################

PLOT_COLORS = ["blue",  "red",    "green", 
               "brown", "orange", "violet", "black"]
