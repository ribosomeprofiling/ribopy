# -*- coding: utf-8 -*-

def read_line_chunks( file_object, max_bytes ):
    """
    Function that lazily reads file in chunks of lines
    it reads upto max_bytes and then till the end of the existing line
    returns the aray of lines read 
    """
    while True:
        data = file_object.readlines( max_bytes)
        if not data:
            break
        yield data