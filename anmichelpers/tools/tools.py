"""
Created on 9/2/2015

@author: victor
"""

import  numpy
import math

def norm(v):
    """
    Numpy compliant norm implementation.
    
    @param v: The vector used to calculate the norm.
    
    @return: A norm or an array of norms.
    """
    return numpy.sqrt(numpy.dot(v,v))

def frec_from_eigenvalue(e_val):
    """
    Calculates the proportional frequency of a given eigenvalue (if it comes from a 
    vibrational study).
    
    @param e_val: The eigenvalue.
    
    @return: The computed frequency (no units).
    """
    return e_val / (2*math.pi)

def ensure_mode_layout(modes):
    """
    If the layout of the modes is flat, it converts it to a (M,N,3) layout.
    
    @param modes: [In/Out] A numpy array containing the modes. 
    
    @return: The same numpy array with a (M,N,3) layout.
    """
    if len(modes.shape) == 3:
        return modes
    else:
        number_of_modes = len(modes)
        number_of_nodes = modes.shape[1] / 3
        modes.resize((number_of_modes, number_of_nodes, 3))
        return modes