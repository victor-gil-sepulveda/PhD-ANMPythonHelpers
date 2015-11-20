import numpy
import anmichelpers.tools.tools as tools
import math

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
        cum_overlap += overlap(mode, mode_range[i])
    return math.sqrt(cum_overlap)

def rmsip(first_mode_range, second_mode_range):
    """
    Root mean square inner product. 
    Indicates how well the motion space spanned by the first range of modes is represented by the
    second range of modes.
    
    "An Analysis of Core Deformations in Protein Superfamilies" Alejandra Leo-Macias,
    Pedro Lopez-Romero, Dmitry Lupyan, Daniel Zerbino and Angel R. Ortiz
    Biophys J. 2005 Feb; 88(2): 1291-1299.
    
    Must vary between 0 and 1
    http://www.biomedcentral.com/1471-2105/15/399
    """
    D = len(first_mode_range)
    K = len(second_mode_range)
    rmsip_val = 0.
    for i in range(D):
        for j in range(K):
#             dot_val = numpy.dot(first_mode_range[i], second_mode_range[j])
#             rmsip_val += dot_val*dot_val
            ovp = overlap(first_mode_range[i], second_mode_range[j])
            rmsip_val += ovp*ovp
    return math.sqrt(rmsip_val / float(D))


def degree_of_collectivity(mode):
    """
    http://peds.oxfordjournals.org/content/14/1/1.long
    """
    alpha = None
    # Calculate "displacements of the mode"
    
    # Calculate alpha
    
    # Calculate the degree of collectivity
