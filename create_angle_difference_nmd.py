'''
Created on 04/05/2015

@author: user
'''
from optparse import OptionParser
import numpy
from anmichelpers.tools.measure import calculate_angle_differences
from anmichelpers.writers.pronmd import ProdyNMDWriter

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f", dest="_from")
    parser.add_option("-t", dest="_to")
    
    parser.add_option("-o", dest="output")
    (options, args) = parser.parse_args()
    
    from_file_v = numpy.loadtxt(options._from)
    to_file_v = numpy.loadtxt(options._to)
    
    difference = calculate_angle_differences(to_file_v, from_file_v)
    
    ##############
    for i in range(len(difference)):
        if not i in range(98*2, 110*2):
            difference[i] = 0.
    ##############
    ProdyNMDWriter.write(options.output, [1.0], numpy.array([difference]),  {"type":"ic:pca"})
    