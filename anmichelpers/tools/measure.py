"""
Created on 9/2/2015

@author: victor
"""

import numpy
from anmichelpers.tools.tools import ensure_modes_layout
import math
try:
    from pyRMSD.RMSDCalculator import RMSDCalculator
except:
    print "[WARNING] pyRMSD was not found. Some functions cannot be used"
    

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

def calculate_angle(v1, v2):
    norm_v1 = math.sqrt(numpy.dot(v1,v1))
    norm_v2 = math.sqrt(numpy.dot(v2,v2))
    
    _cos =  numpy.dot(v1,v2)/(norm_v1*norm_v2)
    
    if _cos > 1: _cos = 1
    if _cos < -1: _cos = -1
    
    return math.acos(_cos)

def calculate_mode_angles(mode1,mode2):
    mode1_3t = numpy.resize(mode1, (len(mode1)/3, 3))
    mode2_3t = numpy.resize(mode2, (len(mode2)/3, 3))
    assert len(mode1_3t) == len(mode2_3t),\
     "[ERROR] calculate_mode_angles - mode lengths must be equal (%s vs %s)."%(str(mode1_3t.shape), str(mode2_3t.shape))
#    numpy.savetxt("mode1.txt", mode1_3t)
#    numpy.savetxt("mode2.txt", mode1_3t)
    print len(mode1_3t)
    angles_1 = []
    for i in range(len(mode1_3t)):
        angles_1.append(calculate_angle(mode1_3t[i], mode2_3t[i]))
    
    angles_2 = []
    mode2_3t *= -1
    for i in range(len(mode1_3t)):
        angles_2.append(calculate_angle(mode1_3t[i], mode2_3t[i]))
        
    if abs(numpy.sum(angles_1)) > abs(numpy.sum(angles_2)):
        return numpy.array(angles_2)
    else:
        return numpy.array(angles_1)
    
def calculate_mode_magnitudes(mode):
    mode_3t = numpy.resize(mode, (len(mode)/3, 3))
    
    magnitudes = []
    for v in mode_3t:
        magnitudes.append(math.sqrt(numpy.dot(v,v)))
    
    return numpy.array(magnitudes)

