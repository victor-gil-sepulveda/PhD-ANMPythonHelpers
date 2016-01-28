'''
Created on Jan 28, 2016

@author: victor
'''
from optparse import OptionParser
from anmichelpers.parsers.pronmd import ProdyNMDParser
import numpy
from anmichelpers.comparison.comparison import overlap

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--ic", dest="ic_modes")
    parser.add_option("--cc",  dest="cc_modes")
    
    (options, args) = parser.parse_args()
    
    _, cc_eigenvectors, _ = ProdyNMDParser.read(options.cc_modes)
    _, ic_eigenvectors, _ = ProdyNMDParser.read(options.ic_modes)
    
    num_cc, num_ic = len(cc_eigenvectors), len(ic_eigenvectors)
    
    overlaps = numpy.zeros((num_cc, num_ic))
    for i in range(num_cc):
        print "%d\t"%i,
        for j in range(num_ic):
            overlaps[i,j] = overlap(cc_eigenvectors[i], ic_eigenvectors[j])
            print "%.3f\t"%overlaps[i,j],
        print
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.heatmap(overlaps,annot=True)
    plt.show()