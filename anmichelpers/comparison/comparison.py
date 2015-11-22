import numpy
import anmichelpers.tools.tools as tools
import math
from anmichelpers.tools.tools import norm
from math import exp

# For all measures see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2350220/

def overlap(mode_i, mode_j):
    """
    Measure of the similarity between 2 modes.
    Overlap value is in the range [0,1] (1 means maximum overlap / similarity)
    Tama and Sanejouand 2001. Conformational change of proteins arising from normal mode calculations.
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
    
    Must vary between 0 and 1
    http://www.biomedcentral.com/1471-2105/15/399
    Also:  http://www.biomedcentral.com/1471-2105/14/183
    """
    D = len(first_mode_range)
    K = len(second_mode_range)
    rmsip_val = 0.
    for i in range(D):
        for j in range(K):
            ovp = overlap(first_mode_range[i], second_mode_range[j])
            rmsip_val += ovp*ovp
    return math.sqrt(rmsip_val / float(D))

def degree_of_collectivity(mode, normalize = False):
    """
    http://peds.oxfordjournals.org/content/14/1/1.long
    Tama and Sanejouand 2001. Conformational change of proteins arising from normal mode calculations.

    From th article:
    A measure of how collective a protein motion is was proposed by Bruschweiler (1995). In the present study, it was used in order 
    to estimate the degree of collectivity of each conformational change considered, reflecting the number of atoms which are 
    significantly affected during the conformational change. This degree of collectivity, k, is defined as being proportional to 
    the exponential of the 'information entropy' embedded in vector inc_R.
    It is confined to the interval between N^-1 and 1. If k = 1, the conformational change is maximally collective.
    
    
    Bruschweiler (1995): http://scitation.aip.org/content/aip/journal/jcp/102/8/10.1063/1.469213
    
    
    \kappa = \frac{1}{N} exp \left( -\sum^N_i \alpha \Delta R_i^2 log \left( \alpha \Delta R_i^2 \right) \right)
    """

    # Calculate "displacements of the mode"
    R = numpy.reshape(mode,(len(mode)/3,3))
    N = len(R)
#     print "Degree of collectivity. Max. value: 1 Min. value: ",1./N
    inc_R = norm(R)
    # Calculate alpha
    sum_ri_sq = numpy.dot(inc_R, inc_R)
    alpha = 1./sum_ri_sq
    # Calculate the degree of collectivity
    alpha_ri2 = alpha*(inc_R*inc_R)
    log_alpha_ri2 = numpy.log(alpha_ri2)
    k = (1./N)*exp(-(alpha_ri2*log_alpha_ri2).sum())
    if not normalize:
        return k
    else:
        min_val = 1./N
        return  (k-min_val) / (1-min_val)

    
# a = [1,2,3,4,5,6] -> degree_of_collectivity = 0.76810859305 (range: 1, 0.5)
# normed -> 0.5362171861
