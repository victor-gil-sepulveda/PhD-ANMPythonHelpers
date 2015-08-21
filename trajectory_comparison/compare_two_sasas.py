"""
Created on Aug 21, 2015

@author: victor
"""
import sys
import numpy
import matplotlib.pyplot as plt
from scipy.stats import entropy
from numpy.linalg import norm
import math
import os

def load_values(files):
    """
    Loads the sasa values from the files in the files dictionary. Returns the maximum and minimum values.
    Values are loaded into the "files" structure.
    """
    min_val = float("inf")
    max_val =  0.0
    for filename in files:
        files[filename]["values"] = numpy.loadtxt(filename)
        max_val = max(numpy.max(files[filename]["values"]), max_val)
        min_val = min(numpy.min(files[filename]["values"]), min_val) 
    
    return max_val, min_val

def JSD(P, Q):
    """
    Calculates the Jensen-Shannon divergence as a metric (sq_root)
    See: http://www.researchgate.net/publication/3084774_A_new_metric_for_probability_distributions
    """
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return math.sqrt(0.5 * (entropy(_P, _M) + entropy(_Q, _M)))

def smoothed(distribution,small_value = 1.0e-8):
    """
    Applies a smoothing process to the distribution.
    See http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
    for an explanation about the problem and the solution.
     
    @param distribution: distribution to be smoothed
    @param small_value: value to be set to those bins with 0 probability
     
    @return: The smoothed distribution.
    """
    total_number_of_samples = len(distribution)
    samples_in_distrib = numpy.count_nonzero(distribution)
    pc = small_value * (total_number_of_samples - samples_in_distrib) / samples_in_distrib
    smoothed_distrib = numpy.empty(len(distribution))
    for i in range(len(distribution)):
        if distribution[i] == 0:
            smoothed_distrib[i] = small_value
        else:
            smoothed_distrib[i] = distribution[i] - pc
    return numpy.array(smoothed_distrib)

def calculate_distributions(input_files, max_val,  min_val):
    """
    Calculates the histogram (distribution) for the values inside the "files" structure.
    """
    for filename in input_files:
        input_files[filename]["distrib"] = smoothed(numpy.histogram(input_files[filename]["values"], 
                                                            bins = 100, range = (min_val, max_val), 
                                                            normed=True)[0])

if __name__ == '__main__':
    
    files = {
        sys.argv[1]:{},
        sys.argv[2]:{}         
    }
    
    max_val, min_val = load_values(files)
    
    calculate_distributions(files, max_val,  min_val)
    
    print "-----------------------"
    print " JSD: ", JSD(files[files.keys()[0]]["distrib"], 
                           files[files.keys()[1]]["distrib"])
    print "-----------------------"
    
    for i, filename in enumerate(files.keys()):
        alpha = 1.0
        if i == 1:
            alpha = 0.8
                    
        plt.hist(files[filename]["values"], bins = 100, 
                 range = (min_val, max_val), 
                 label = os.path.basename(filename),
                 alpha = alpha)
    
    plt.legend()
    plt.show()
    
    