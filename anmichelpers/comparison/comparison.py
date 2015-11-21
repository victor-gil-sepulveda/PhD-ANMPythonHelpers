import numpy
import anmichelpers.tools.tools as tools
import math

# For all measures see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2350220/

def overlap(mode_i, mode_j):
    """
    Measure of the similarity between 2 modes.
    Overlap value is in the range [0,1] (1 means maximum overlap / similarity)
    Tama and Sanejouand 2001.
    """
    return numpy.abs(numpy.dot(mode_i, mode_j)) / (tools.norm(mode_i) * tools.norm(mode_j))

def cumulative_overlap(mode, mode_range):
    """
    Measure of similarity between one mode and a range of modes.
    
    A value of 1 means a perfect match.
    """
    cum_overlap = 0
    for i in range(len(mode_range)):
        o = overlap(mode, mode_range[i])
        cum_overlap += o*o
    return math.sqrt(cum_overlap)

def rmsip(first_mode_range, second_mode_range):
    """
    Root mean square inner product. 
    Indicates how well the motion space spanned by the first range of modes is represented by the
    second range of modes.
    
    "An Analysis of Core Deformations in Protein Superfamilies" Alejandra Leo-Macias,
    Pedro Lopez-Romero, Dmitry Lupyan, Daniel Zerbino and Angel R. Ortiz
    Biophys J. 2005 Feb; 88(2): 1291-1299.
    
    """
    D = len(first_mode_range)
    K = len(second_mode_range)
    
    rmsip_val = 0.
    for i in range(D):
        for j in range(K):
            print " REPASAAAR!!!!!!!!!!"
            o = numpy.dot(first_mode_range[i], second_mode_range[j])
            rmsip_val += o*o
    return math.sqrt(rmsip_val / D)
