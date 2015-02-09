"""
Created on 9/2/2015

@author: victor
"""
import  anmichelpers.tools.tools as tools
from anmichelpers.tools.tools import ensure_mode_layout
import numpy

def mean_square_fluctuations(eigenvalues, eigenvectors, number_of_nodes, prefactor = 1.8):
    """
    prefactor = 3kBT/ gamma (see Biophys J. 2001 Jan; 80(1): 505-515. and 
    Biophys J. 2005 Feb; 88(2): 1291-1299) (formula from the first, eq 8)
    
    @param eigenvalues: SIGNIFICANT (non zero) eigenvalues (3N-6 in case of cartesian ANM)
    @param eigenvectors: Eigenvectors related to that eigenvalues (indeed eigenvectors of the 
    Kirchoff matrix).
    
    @return: Per node square distance increment (mean square fluctuation) in sq Angstroms.
    """
    ensure_mode_layout(eigenvectors)
    
    nodes_msq = numpy.array([0.] * (number_of_nodes))
    for k in range(number_of_nodes):
        for j in range(len(eigenvalues)): #mode number
            mode = eigenvectors[j]
            nodes_msq[k] += numpy.dot(mode[k], mode[k]) / eigenvalues[j]
    
    return nodes_msq
    
    