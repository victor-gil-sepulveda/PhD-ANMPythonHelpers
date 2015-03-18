'''
Created on 18/03/2015

@author: user
'''
import numpy
from anmichelpers.tools.tools import is_int

def single_string_parser(parts):
    return parts[0]

def string_array_parser(parts):
    return parts
    
def float_array_parser(parts):
    return  numpy.array([float(part) for part in parts])

def int_array_parser(parts):
    return  [int(part) for part in parts]

def mode_parser(parts):
    offset = 0
    if is_int(parts[0]):
        offset = 1
    return float(parts[offset]), float_array_parser(parts[offset+1:])