"""
Created on 18/05/2015

Calculates the dPCA values from the angles got from "calculate structurally aligned angles"

@author: vgil
"""
import numpy
from optparse import OptionParser
from numpy.core.numeric import inf
import math
from pca.dpca import dPCA
from convert_geo_pca_to_nmd import load_sequences,\
    treat_prolines_in_eigenvectors
from anmichelpers.writers.pronmd import ProdyNMDWriter

def angle_quantization(angles, bin_size = 0.1):
    """
    given an array of angles, does the histogram and returns the value of the most
    dense bin.
    """
    if angles != []:
        hist, boundaries = numpy.histogram(angles, 
                                           bins=int((math.pi*2)/bin_size),
                                           range = (-math.pi, math.pi))
        # Find maximum value
        max_index = hist.tolist().index(numpy.max(hist))
        value = (boundaries[max_index] + boundaries[max_index+1]) / 2.
        return value
    else:
        return inf

def non_inf_angles(angles):
    ret_angles = []
    for angle in angles:
        if angle != inf:
            ret_angles.append(angle)
    return ret_angles

def change_inf_by_val(angles, val):
    for i in range(len(angles)):
        if angles[i] == inf:
            angles[i] = val

def recover_gaps(all_angles):
    for angles in all_angles.T:
        angle_mode = angle_quantization(non_inf_angles(angles))
        change_inf_by_val(angles, angle_mode)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-o", dest="output")
    parser.add_option("--geopca", dest="geopca_format", action="store_true", default=False)
    parser.add_option("--ref", dest="ref", type="int")
    
    (options, args) = parser.parse_args()

    # it reads files with phi, psi pairs, even if the first phi is not defined (we expect it)

    if options.geopca_format:    
        #    angles layout:
        #    [[angle 1 str 1, angle 1 str 2, ..., angle 1 str n],
        #     [angle 2 str 1, angle 2 str 2, ..., angle 2 str n],
        #     ...
        #     [angle m str 1, angle m str 2, ..., angle m str n]]
        all_angles = numpy.loadtxt(options.input+".angles", delimiter = ",").T
    else:
        all_angles = numpy.loadtxt(options.input+".angles")
        
    # Recover gaps by angle mode quantization
    recover_gaps(all_angles)
    
    evals, evecs = dPCA.calc(all_angles[:,1:]) # removing first phi column (must be undefined)
        
    sequences = load_sequences(options.input+".seq")
    
    pro_evecs = treat_prolines_in_eigenvectors(sequences[options.ref], evecs)
    
    header = {"type": "ic:pca", 
              "title": options.output,
              "resnames":sequences[options.ref] }
    ProdyNMDWriter.write(options.output, evals, pro_evecs, header)
