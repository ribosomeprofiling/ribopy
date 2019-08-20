# -*- coding: utf-8 -*-
class RiboBaseError(Exception):
    pass

class AnnotationError(RiboBaseError):
	pass

class ExperimentDoesntExist(RiboBaseError):
    pass

class CoverageDoesntExist(RiboBaseError):
    pass

class InvalidName(RiboBaseError):
    pass

class MissingExperiment(RiboBaseError):
    pass

class InvalidLengthRange(RiboBaseError):
    pass
    
class NORNASEQ(RiboBaseError):
    pass

class ReferecenDoesntExist(RiboBaseError):
    pass
    
class TranscriptDoesntExist(RiboBaseError):
    pass

class AliasError(RiboBaseError):
    pass
