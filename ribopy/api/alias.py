# -*- coding: utf-8 -*-
from collections import defaultdict, OrderedDict

from ..core.exceptions import AliasError


def apris_human_alias(x):
    """
    A predefined function for convenience
    It extracts transcript IDs from the full names of appris genomes. 
    """
    return x.split("|")[4]

class ReferenceAlias:
    def __init__(self, renaming_function, reference_names):
        self.reference_names = reference_names
        
        # Mapping dictionary maps original reference names
        # to aliases
        
        # Reverse Mapping dictionary maps aliases
        # to original names
        
        if issubclass(type(renaming_function), dict):
            self.mapping_dict      = renaming_function
            self.aliases           = [ renaming_function[x] \
                                       for x in reference_names ]
        else:
            self.aliases           = tuple( map( renaming_function, 
                                                 reference_names ) )
            self.mapping_dict      = dict( tuple( zip( reference_names,
                                                    self.aliases ) ) )
                                        
        # Checks the renaming.
        # Raises an exception if it encounters any errors.                                            
        self._check_renaming()
        
        # Reverse mapping dict maps aliases to original names
        self.reverse_mapping_dict = dict( tuple( zip( self.aliases,
                                                      reference_names) ) )
        
    def _check_renaming(self):
        if len(set(self.reference_names)) < len(self.reference_names):
            raise AliasError("Reference names must be distinct.")
            
        if len(set(self.aliases)) < len(self.aliases):
            raise AliasError("Renaming function must be one-to-one!")
            
        for x in self.aliases:
            if not issubclass(type(x), str):
                print(type(x))
                raise AliasError("The aliases must be strings! {}".format(x))
                
            if len(x) == 0:
                raise AliasError("The length of each alias must be positive.")
        
        
        
    def get_original_name(self, alias):
        original_name = self.reverse_mapping_dict.get(alias, None)
        
        if original_name is None:
            raise AliasError("No such alias exists: " + alias)
        else:
            return original_name
        
        
    def get_alias(self, original_name):
        alias = self.mapping_dict.get(original_name, None)
        
        if alias is None:
            raise AliasError("No such reference exists: " + original_name)
        else:
            return alias

############################################################################
            
        
        
        
