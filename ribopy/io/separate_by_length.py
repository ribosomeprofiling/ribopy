# -*- coding: utf-8 -*-
import os
##############################################################################

def separate_by_length(input_stream, 
                       length_max, length_min, file_prefix = ""):
    """
    Given a bed file handle and a length range [min length,  max length] 
    it partitions the file into bed files for each length in the given range
    file_prefix is pre_appended to each file name.
    So, for example, for length_min = 3, length_max=5, 
    3 files for lengths 3,4,5 are going to be created

    Parameters
    ----------
    input_stream : A stream coming from a bed file
    length_max   : maximum length (inclusive)
    length_min   : minimum length (inclusive)
    file_prefix  : Will be preappended to output file names

    Returns
    -------
    A list of separated file paths.

    """

    output_file_base = file_prefix + ".len_"

    assert length_min > 0
    assert length_max >= length_min

    file_handles = list()
    separated_files = list()

    for i in range(length_min, length_max + 1):
        this_file_path = output_file_base + str(i)
        separated_files.append(this_file_path)
        file_handles.append( open(this_file_path, "w") )  

    try:
        for this_line in input_stream:
            this_line = this_line.strip()
            line_contents = this_line.split()
            
            if len(line_contents) < 6:
                continue

            read_length = int( line_contents[2]) - int( line_contents[1] )
            if read_length < length_min or read_length > length_max:
                continue

            output_handle = file_handles[read_length - length_min]
            print(this_line, file=output_handle)

    except Exception as e:
        [ os.remove( output_file_base + str(i) ) \
             for i in range(length_min, length_max + 1)  ] 
        raise IOError("An error happened while separating reads by length.\n" + str(e))

    finally:
        input_stream.close()
        [ f.close() for f in file_handles ]

    return separated_files
