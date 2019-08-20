# -*- coding: utf-8 -*-



"""
Target RNA-Seq Data
Region  UTR5  UTR5_Junction  CDS  UTR3_junction   UTR3 
GAPDH   1        0            2       1            3
VEGFA   0        1            3       0            1
MYC     0        1            0       0            2
"""

## Note that we assume LEFTSPAN=3
## and RIGHSPAN=2

RNASEQ_READS = \
"""GAPDH\t1\t3\tUTR5_1\t0\t+
GAPDH\t9\t11\tCDS_1\t0\t+
GAPDH\t10\t11\tCDS_2\t0\t+
GAPDH\t13\t16\tUTR3_junction_1\t0\t+
GAPDH\t18\t19\tUTR3_1\t0\t+
GAPDH\t18\t19\tUTR3_2\t0\t+
GAPDH\t18\t19\tUTR3_3\t0\t+
VEGFA\t2\t4\tUTR5_junction_1\t0\t+
VEGFA\t7\t10\tCDS_1\t0\t+
VEGFA\t7\t11\tCDS_2\t0\t+
VEGFA\t8\t12\tCDS_3\t0\t+
VEGFA\t20\t21\tUTR3_1\t0\t+
MYC\t1\t4\tUTR5_junction_1\t0\t+
MYC\t14\t16\tUTR3_1\t0\t+
MYC\t14\t16\tUTR3_2\t0\t+"""


####################################################################

"""
Target RNA-Seq Data (RNASEQ_READS_2)
Region  UTR5  UTR5_Junction  CDS  UTR3_junction   UTR3 
GAPDH   0        0            0       0            0
VEGFA   0        0            0       2            0
MYC     0        0            0       0            0
"""

RNASEQ_READS_2 = \
"""VEGFA\t14\t16\tUTR3_junction_1\t0\t+
VEGFA\t16\t18\tUTR3_junction_2\t0\t+"""

#######################################################################

RNASEQ_tsv_1 = \
"\tUTR5\tUTR5_junction\tCDS\tUTR3_junction\tUTR3\n"+\
"GAPDH\t0\t7.2\t3\t0\t2\n"+\
"VEGFA\t1\t0\t5\t3\t1\n"+\
"MYC\t0\t0\t2\t0\t0"

RNASEQ_tsv_2 = \
"VEGFA\t1.9\t4.8\t5.7\t2.6\t1.5\n"+\
"GAPDH\t10\t17.2\t27\t37\t47"

RNASEQ_tsv_3 = \
"VEGFA\t10\t20\t30\n"+\
"GAPDH\t100\t200\t300"
