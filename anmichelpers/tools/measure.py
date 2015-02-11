"""
Created on 9/2/2015

@author: victor
"""

import numpy
from anmichelpers.tools.tools import ensure_modes_layout
from pyRMSD.RMSDCalculator import RMSDCalculator

def mean_square_fluctuations(eigenvalues, eigenvectors, number_of_nodes, prefactor = 1.8):
    """
    prefactor = 3kBT/ gamma (see Biophys J. 2001 Jan; 80(1): 505-515. and 
    Biophys J. 2005 Feb; 88(2): 1291-1299) (formula from the first, eq 8)
    
    @param eigenvalues: SIGNIFICANT (non zero) eigenvalues (3N-6 in case of cartesian ANM)
    @param eigenvectors: Eigenvectors related to that eigenvalues (indeed eigenvectors of the 
    Kirchoff matrix).
    
    @return: Per node square distance increment (mean square fluctuation) in sq. Angstroms.
    """
    eigenvectors = ensure_modes_layout(eigenvectors)
    
    nodes_msq = numpy.array([0.] * (number_of_nodes))
    for k in range(number_of_nodes):
        for j in range(len(eigenvalues)): #mode number
            mode = eigenvectors[j]
            nodes_msq[k] += numpy.dot(mode[k], mode[k]) / eigenvalues[j]
    
    return nodes_msq

def rmsf(prody_pdb):
    """
    Calculates CA Root Mean Square Fluctuation.
    
    @param prody_pdb: A prody pdb data structure.
    
    @return: An array with the per-residue CA rmsf.
    """
    # Calculate the superposition
    ca_coords = prody_pdb.select("name CA").getCoordsets()

    calculator = RMSDCalculator(calculatorType = "QTRFIT_SERIAL_CALCULATOR",
                                fittingCoordsets = ca_coords)
    
    calculator.iterativeSuperposition()
    
    # Calculate the actual rmsf
    mean_conformation = ca_coords.mean(0)
    
    ssqf = numpy.zeros(mean_conformation.shape)

    for conf in ca_coords:
            ssqf += (conf - mean_conformation) ** 2
            
    return (ssqf.sum(1) / ca_coords.shape[0])**0.5

