# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:54:22 2022

@author: willi
"""

class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class UnknownElementError(Error):
    """Exception raised for errors in converting models to c.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class SnapGeneMissingError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class InvalidCharError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message,index):
        self.message = message

class FileTypeNotRecognizedError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class PathDoesNotExistError(Error):
    """Exception raised for when trying to save a GB file to a directory 
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class AscNumDoesNotExistError(Error):
    """Exception raised for when trying to pull a gb from an ascession number
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class InvalidSequenceLengthError(Error):
    """Exception raised for when a sequence is not a multiple of 3

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class MismatchError(Error):
    """Exception raised for when trying to pull a gb from an ascession number
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class UnrecognizedAAError(Error):
    """Exception raised for when trying to pull a gb from an ascession number
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class UnrecognizedCodonError(Error):
    """Exception raised for when trying to pull a gb from an ascession number
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class UnrecognizedFlagError(Error):
    """Exception raised for when trying to pull a gb from an ascession number
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class UnrecognizedNormalizationError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
                


class NegativeRateError(Error):
    """Exception raised for errors in the simulation where negative rates appear

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class NegativeResourcesError(Error):
    """Exception raised for errors in the simulation where negative resources appear

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class StatesError(Error):
    """Exception raised for errors in the simulation where states are not equal 
    to 0 or 1

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class SmallTimestepError(Error):
    """Exception raised for errors in the simulation where time steps are less than 1e-8

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message



class EigenMissingError(Error):
    """Exception raised for when an eigen instillation cannot be found

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class UnknownElementError(Error):
    """Exception raised for when an unknown element to convert to C++ was used

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class ExistenceError(Error):
    """Exception raised for when a requesting making a model that already
    exists without overwite == True

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class ModelNameError(Error):
    """Exception raised for when a requesting making a model that already
    exists without overwite == True

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message



class MisMatchedBrackets(Error):
    """Exception raised for when a requesting making a model that already
    exists without overwite == True

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class MisMatchedParenthesis(Error):
    """Exception raised for when a requesting making a model that already
    exists without overwite == True

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
