# -*- coding: utf-8 -*-

import numpy as np

###########################################
REF_LEN_FILE     = "ref_lengths.tsv"
TEST_RIBO        = "test_create.ribo"
TEST_RIBO_2      = "test_2.ribo"
ANNOTATION_FILE  = "annotation.bed"
ALIGNMENT_FILE_1 = "alignment_file_1.bed"
ALIGNMENT_FILE_2 = "alignment_file_2.bed"
METADATA_FILE_1  = "exp_metadata.yaml"
RIBO_META_FILE_1 = "ribo_metadata.yaml"
RNASEQ_FILE_1    = "rnaseq_1.tsv"
RNASEQ_FILE_2    = "rnaseq_2.tsv"



#TRANSCRIPT_LENGTHS=\
"""GAPDH    20
VEGFA   22
MYC 17"""

TRANSCRIPT_LENGTHS=\
"GAPDH\t20\nVEGFA\t22\nMYC\t17"


TRANSCRIPT_ANNOTATION=\
"""GAPDH    0   5  UTR5    0   +
GAPDH   5  15    CDS 0   +
GAPDH   15  20    UTR3 0   +
VEGFA    0   4  UTR5    0   +
VEGFA   4  16    CDS 0   +
VEGFA   16  22    UTR3 0   +
MYC    0   3  UTR5    0   +
MYC   3  10    CDS 0   +
MYC   10  17    UTR3 0   +"""


READ_SET_1=\
"""GAPDH 1  3  len_2_UTR5_1 0   +
GAPDH 3  5  len_2_UTR5_junc_1 0   +
GAPDH 4  6  len_2_UTR5_junc_2 0   +
GAPDH 4  6  len_2_UTR5_junc_3 0   +
GAPDH 6  8  len_2_UTR5_junc_4 0   +
GAPDH 7  9  len_2_UTR5_junc_5 0   +
GAPDH 9  11 len_2_CDS_1 0   +
GAPDH 10  12  len_2_CDS_2 0   +
GAPDH 13  15  len_2_UTR3_junc_1 0   +
GAPDH 13  15  len_2_UTR3_junc_2 0   +
GAPDH 14  16  len_2_UTR3_junc_3 0   +
GAPDH 15  17  len_2_UTR3_junc_4 0   +
GAPDH 18  20  len_2_UTR3_1 0   +
GAPDH 19  21  len_2_UTR3_1 0   +
GAPDH 0  3  len_3_UTR5_junc_1 0   +
GAPDH 1  4  len_3_UTR5_junc_2 0   +
GAPDH 3  6  len_3_UTR5_1 0   +
GAPDH 3  6  len_3_UTR5_2 0   +
GAPDH 3  6  len_3_UTR5_3 0   +
GAPDH 4  7  len_3_UTR5_4 0   +
GAPDH 4  7  len_3_UTR5_5 0   +
GAPDH 5  8  len_3_UTR5_6 0   +
GAPDH 6  9  len_3_UTR5_7 0   +
GAPDH 8  11  len_3_CDS_1 0   +
GAPDH 8  11  len_3_CDS_2 0   +
GAPDH 10  13  len_3_CDS_3 0   +
GAPDH 13  16  len_3_UTR3_junc_1 0   +
GAPDH 13  16  len_3_UTR3_junc_2 0   +
GAPDH 13  16  len_3_UTR3_junc_3 0   +
GAPDH 13  16  len_3_UTR3_junc_4 0   +
GAPDH 14  17  len_3_UTR3_junc_5 0   +
GAPDH 14  17  len_3_UTR3_junc_6 0   +
GAPDH 15  18  len_3_UTR3_junc_7 0   +
GAPDH 15  18  len_3_UTR3_junc_8 0   +
GAPDH 16  19  len_3_UTR3_junc_9 0   +
GAPDH 3  7  len_4_UTR5_junc_1 0   +
GAPDH 3  7  len_4_UTR5_junc_2 0   +
GAPDH 4  8  len_4_UTR5_junc_3 0   +
GAPDH 4  8  len_4_UTR5_junc_4 0   +
GAPDH 4  8  len_4_UTR5_junc_5 0   +
GAPDH 5  9  len_4_UTR5_junc_6 0   +
GAPDH 5  9  len_4_UTR5_junc_7 0   +
GAPDH 5  9  len_4_UTR5_junc_8 0   +
GAPDH 5  9  len_4_UTR5_junc_9 0   +
GAPDH 7  11  len_4_UTR5_junc_10 0   +
GAPDH 9  13  len_4_CDS_1 0   +
GAPDH 9  13  len_4_CDS_2 0   +
GAPDH 10  14  len_4_CDS_3 0   +
GAPDH 11  15  len_4_CDS_4 0   +
GAPDH 11  15  len_4_CDS_5 0   +
GAPDH 11  15  len_4_CDS_6 0   +
GAPDH 13  17  len_4_UTR3_junc_1 0   +
GAPDH 13  17  len_4_UTR3_junc_2 0   +
GAPDH 13  17  len_4_UTR3_junc_3 0   +
GAPDH 13  17  len_4_UTR3_junc_4 0   +
GAPDH 14  18  len_4_UTR3_junc_5 0   +
GAPDH 14  18  len_4_UTR3_junc_6 0   +
GAPDH 15  19  len_4_UTR3_junc_7 0   +
GAPDH 15  19  len_4_UTR3_junc_8 0   +
GAPDH 16  20  len_4_UTR3_junc_9 0   +
GAPDH 1  6  len_5_UTR5_1 0   +
GAPDH 3  8  len_5_UTR5_junc_1 0   +
GAPDH 5  10  len_5_UTR5_junc_1 0   +
GAPDH 8  13  len_5_CDS_ 1   +
GAPDH 10  15  len_5_CDS_2 0   +
GAPDH 14  19  len_5_UTR5_junc_1 0   +
VEGFA 3  5  len_2_UTR5_junc_1 0   +
VEGFA 3  5  len_2_UTR5_junc_2 0   +
VEGFA 4  6  len_2_UTR5_junc_3 0   +
VEGFA 10  12  len_2_CDS_1 0   +
VEGFA 16  18  len_2_UTR3_junc_1 0   +
VEGFA 3  6  len_3_UTR5_junc_1 0   +
VEGFA 3  6  len_3_UTR5_junc_2 0   +
VEGFA 3  6  len_3_UTR5_junc_3 0   +
VEGFA 4  7  len_3_UTR5_junc_4 0   +
VEGFA 4  7  len_3_UTR5_junc_5 0   +
VEGFA 5  8  len_3_UTR5_junc_6 0   +
VEGFA 10  13  len_3_CDS_1 0   +
VEGFA 10  13  len_3_CDS_2 0   +
VEGFA 16  19  len_3_UTR3_junc_1 0   +
VEGFA 16  19  len_3_UTR3_junc_2 0   +
VEGFA 3  7  len_4_UTR5_junc_1 0   +
VEGFA 3  7  len_4_UTR5_junc_2 0   +
VEGFA 3  7  len_4_UTR5_junc_3 0   +
VEGFA 3  7  len_4_UTR5_junc_4 0   +
VEGFA 4  8  len_4_UTR5_junc_5 0   +
VEGFA 4  8  len_4_UTR5_junc_6 0   +
VEGFA 4  8  len_4_UTR5_junc_7 0   +
VEGFA 5  9  len_4_UTR5_junc_8 0   +
VEGFA 5  9  len_4_UTR5_junc_9 0   +
VEGFA 10  14  len_4_CDS_1 0   +
VEGFA 10  14  len_4_CDS_2 0   +
VEGFA 10  14  len_4_CDS_3 0   +
VEGFA 3  8  len_5_UTR5_junc_1 0   +
VEGFA 10  15  len_5_CDS_1 0   +
VEGFA 10  15  len_5_CDS_2 0   +
VEGFA 10  15  len_5_CDS_3 0   +
VEGFA 10  15  len_5_CDS_4 0   +
VEGFA 16  21  len_5_UTR3_junc_1 0   +
VEGFA 16  21  len_5_UTR3_junc_2 0   +
MYC 10  12  len_2_UTR3_junc_1 0   +
MYC 10  12  len_2_UTR3_junc_2 0   +
MYC 0  3  len_3_UTR5_junc_1 0   +
MYC 6  9  len_3_CDS_1 0   +
MYC 10  13  len_3_UTR3_junc_1 0   +
MYC 6  10  len_4_CDS_1 0   +
MYC 6  10  len_4_CDS_2 0   +
MYC 10  14  len_4_UTR3_junc_1 0   +
MYC 10  14  len_4_UTR3_junc_2 0   +
MYC 10  14  len_4_UTR3_junc_3 0   +
MYC 0  5  len_5_UTR5_junc_1 0   +
MYC 10  15  len_UTR3_junc_1 0   +
MYC 10  15  len_UTR3_junc_2 0   +
MYC 10  15  len_UTR3_junc_3 0   +
MYC 10  15  len_UTR3_junc_4 0   +
MYC 10  15  len_UTR3_junc_5 0   +
MYC 10  15  len_UTR3_junc_6 0   +
MYC 10  15  len_UTR3_junc_7 0   +"""


#########################################


#########################################

### Expected Region Counts for READ_SET_1
# Order is UTR5, UTR5_junc, CDS, UTR3_junc, UTR3

## LENGTH 2

GAPDH_counts_length_2 = ( 1, 5, 2, 4, 2 )
VEGFA_counts_length_2 = ( 0, 3, 1, 1, 0 )
MYC_counts_length_2   = ( 0, 0, 0, 2, 0 )

expected_counts_length_2 = np.array( (GAPDH_counts_length_2, VEGFA_counts_length_2, 
                             MYC_counts_length_2) )

ACTUAL_START_SITE_COVERAGE_length_2 = np.array( [
                                [1, 2, 0, 1, 1], 
                                [0, 2, 1, 0, 0],
                                [0, 0, 0, 0, 0]] )

ACTUAL_STOP_SITE_COVERAGE_length_2 = np.array( [ 
                                [2, 1, 1, 0, 0],
                                [0, 0, 1, 0, 0],
                                [0, 0, 2, 0, 0]] )

#########################################

## LENGTH 3

GAPDH_counts_length_3 = ( 2, 7, 3, 9 , 0 )
VEGFA_counts_length_3 = ( 0, 6, 2, 2, 0 )
MYC_counts_length_3   = ( 0, 1, 1, 1, 0 )

expected_counts_length_3 = np.array(( GAPDH_counts_length_3, VEGFA_counts_length_3, 
                             MYC_counts_length_3) )

ACTUAL_START_SITE_COVERAGE_length_3 = np.array( [
                                [3, 2, 1, 1, 0], 
                                [0, 3, 2, 1, 0],
                                [0, 0, 0, 0, 0]] )

ACTUAL_STOP_SITE_COVERAGE_length_3 = np.array( [ 
                                [4, 2, 2, 1, 0],
                                [0, 0, 2, 0, 0],
                                [0, 0, 1, 0, 0]] )

#########################################

## LENGTH 4

GAPDH_counts_length_4 = ( 0, 10, 6, 9, 0 )
VEGFA_counts_length_4 = ( 0, 9, 3, 0, 0 )
MYC_counts_length_4   = ( 0, 0, 2, 3, 0 )

expected_counts_length_4 = np.array(( GAPDH_counts_length_4, VEGFA_counts_length_4, 
                             MYC_counts_length_4))

ACTUAL_START_SITE_COVERAGE_length_4 = np.array( [
                                [2, 3, 4, 0, 1], 
                                [0, 4, 3, 2, 0],
                                [0, 0, 0, 0, 0]] )

ACTUAL_STOP_SITE_COVERAGE_length_4 = np.array( [ 
                                        [4, 2, 2, 1, 0],
                                        [0, 0, 0, 0, 0],
                                        [0, 0, 3, 0, 0]] )

#########################################

## LENGTH 5

GAPDH_counts_length_5 = ( 1, 2, 2, 1, 0 )
VEGFA_counts_length_5 = ( 0, 1, 4, 2, 0 )
MYC_counts_length_5   = ( 0, 1, 0, 7, 0 )

expected_counts_length_5 = np.array(( GAPDH_counts_length_5, VEGFA_counts_length_5, 
                             MYC_counts_length_5))

ACTUAL_START_SITE_COVERAGE_length_5 = np.array( [
                                [1, 0, 1, 0, 0], 
                                [0, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0]] )

ACTUAL_STOP_SITE_COVERAGE_length_5 = np.array( [ 
                                        [0, 1, 0, 0, 0],
                                        [0, 0, 2, 0, 0],
                                        [0, 0, 7, 0, 0]] )

#########################################



LEFT_SPAN  = 3
RIGHT_SPAN = 2
METAGENE_RADIUS = 2

ANNOTATION = [\
               ([0, 5], [5, 15], [15, 20]),\
               ([0, 4], [4, 16], [16, 22]),\
               ([0, 3], [3, 10], [10, 17]),\
              ]

#########################################
#####  Expected EXTENDED ANNOTATION 

GAPDH_regions = ( [0, 2], [2, 7], [7, 12], [12, 18], [18, 20] )
VEGFA_regions = ( [0, 1], [1, 7], [7, 13], [13, 19], [19, 22] )
MYC_regions   = ( [0, 0], [0, 6], [6, 7],  [7, 13],  [13,17] )

extended_boundary_regions = ( GAPDH_regions, VEGFA_regions,
                               MYC_regions )
#########################################

METADATA_EXPERIMENT_DICT = {
    "cell_line": "HeLa",
    "digestion_enzyme": "HindIII",
    "digestion_duration": "5 min",
    "external_link" : "https://www.encodeproject.org/",
    "link" : "https://www.encodeproject.org/"
}

METADATA_EXPERIMENT_STR_1 =\
"""cell_line: HeLa
digestion_enzyme: HindIII
digestion_duration: 5 min
link: https://www.encodeproject.org/"""

RIBO_METADATA_STR_1 = \
"""aligner: bowtie2
mapq_threshold : 3
deduplicated: True"""


############################################
############################################

READ_SET_2=\
"""GAPDH 8  11  len_3_CDS_1 0   +
GAPDH 8  11  len_3_CDS_2 0   +
GAPDH 8  11  len_3_CDS_3 0   +
GAPDH 8  11  len_3_CDS_4 0   +
GAPDH 8  11  len_3_CDS_5 0   +
VEGFA 10  13  len_3_CDS_1 0   +
VEGFA 10  13  len_3_CDS_2 0   +
VEGFA 10  13  len_3_CDS_3 0   +
MYC 0  3 len_3_UTR5_junc_1 0   +
MYC 0  3 len_3_UTR5_junc_2 0   +
MYC 0  3 len_3_UTR5_junc_3 0   +
MYC 0  3 len_3_UTR5_junc_4 0   +"""

# Expected Counts for 
## LENGTH 3

READ_SET_2_GAPDH_counts_length_3 = ( 0, 0, 5, 0, 0 )
READ_SET_2_VEGFA_counts_length_3 = ( 0, 0, 3, 0, 0 )
READ_SET_2_MYC_counts_length_3   = ( 4, 0, 0, 0, 0 )

expected_counts_length_3 = np.array(( GAPDH_counts_length_3, VEGFA_counts_length_3, 
                             MYC_counts_length_3) )


def _get_transcripts( file_in_string ):
    rows = file_in_string.split("\n")
    pairs = tuple( map( lambda x: x.split(), rows ) )
    ref_names = tuple( map( lambda x: x[0], pairs ) )
    ref_lengths = tuple( map( lambda x: x[1], pairs ) )
    ref_lengths = tuple (map( int, ref_lengths ) )
    return (ref_names, ref_lengths)

###############################################

"""
RNASEQ_DATA_1 = \
"GAPDH\t2.89\nVEGFA\t15.46\nMYC\t8"

RNASEQ_DATA_2 = \
"GAPDH\t6.42\nVEGFA\t1.37\nMYC\t0.06"
"""
