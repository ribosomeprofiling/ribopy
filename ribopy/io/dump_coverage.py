# -*- coding: utf-8 -*-

###########################################################

def _print_coverage_bg(cov_dict, file_handle):
    """
    Print cov_dict in bedgraph format.
    """

    for transcript, t_cov in cov_dict.items():
        for i, i_cov in enumerate(t_cov):
            if i_cov == 0:
                continue
            line = "\t".join( 
                         tuple( 
                            map( str, (transcript, 
                                       i, 
                                       i + 1, 
                                       i_cov)  )  ) )
            print(line, file = file_handle)


def _print_coverage_tsv(cov_dict, file_handle, 
                                   sep = ","):
    """
    Print cov_dict in tab separated file.

    The first column is the transcript name and the
    second column is the coverage values separated by commas.
    """

    for transcript, t_cov in cov_dict.items():
        
        line = transcript + "\t" + ",".\
                  join(tuple( map(str, t_cov) ))                              
        print(line, file = file_handle)


def print_coverage(  cov_dict, 
                     file_handle, 
                     file_format = "bg", 
                     sep         = ","  ):
    """
    Prints the given transcript coverages to a given file handle

    The keys of the coverage dictionary (cov_dict)
    are transcript names. 
    The values are numpy arrays giving the coverage of the transcript, 
    """

    if file_format == "bg":
        _print_coverage_bg( cov_dict    = cov_dict, 
                            file_handle = file_handle)
    else:
        _print_coverage_tsv(cov_dict    = cov_dict, 
                            file_handle = file_handle, 
                            sep         = sep)
