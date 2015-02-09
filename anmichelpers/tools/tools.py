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
    if len(v.shape) == 1:
        return numpy.sqrt(numpy.dot(v,v))
    elif len(v.shape) == 2:
        norms = []
        for i in range(len(v)):
            norms.append(norm(v[i]))
        return norms
    else:
        return None

def frec_from_eigenvalue(e_val):
    """
    Calculates the proportional frequency of a given eigenvalue (if it comes from a 
    vibrational study).
    
    @param e_val: The eigenvalue.
    
    @return: The computed frequency (no units).
    """
    return e_val / (2*math.pi)

def ensure_modes_layout(modes):
    """
    If the layout of the modes is flat, it converts it to a (M,N,3) layout.
    
    @param modes: [In/Out] A numpy array containing all the modes. 
    
    @return: The same numpy array with a (M,N,3) layout or (N,3) 
    """
    if len(modes.shape) == 3:
        return modes
    elif len(modes.shape) == 2:
        number_of_modes = len(modes)
        number_of_nodes = modes.shape[1] / 3
        return numpy.reshape(modes, (number_of_modes, number_of_nodes, 3))
    else:
        raise ValueError("The array has an unexpected size")
        