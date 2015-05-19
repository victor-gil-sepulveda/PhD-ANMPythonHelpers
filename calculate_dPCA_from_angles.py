"""
Created on 18/05/2015

Calculates the dPCA values from the angles got from "calculate structurally aligned angles"

@author: vgil
"""
import numpy
from pca.dpca import dPCA
from anmichelpers.writers.pronmd import ProdyNMDWriter
from optparse import OptionParser
from convert_geo_pca_to_nmd import treat_prolines_in_eivgenvectors,\
    load_sequences

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-s", dest="sequence")
    parser.add_option("-o", dest="output")
    (options, args) = parser.parse_args()
    
    angles = numpy.loadtxt(options.input, delimiter = ",")
#    angles layout:
#    [[angle 1 str 1, angle 1 str 2, ..., angle 1 str n],
#     [angle 2 str 1, angle 2 str 2, ..., angle 2 str n],
#     ...
#     [angle m str 1, angle m str 2, ..., angle m str n]]


    evals, evecs = dPCA.calc(angles.T)

    sequences = load_sequences(options.sequence)
    
    pro_evecs = treat_prolines_in_eivgenvectors(sequences[0], evecs)
    
    header = {"type": "ic:pca"}
    ProdyNMDWriter.write(options.output, evals, pro_evecs, header)
