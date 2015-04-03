"""
Created on 18/03/2015

@author: vgil
"""

from optparse import OptionParser
import matplotlib.pyplot as plt
from anmichelpers.parsers.pronmd import ProdyNMDParser
import anmichelpers.tools.measure as measure
from anmichelpers.tools.tools import normalize
from anmichelpers.tools.atoms import get_CA_modes
import numpy

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--i1", dest="input1")
    parser.add_option("--i2", dest="input2")
    parser.add_option("-o", dest="output")
    parser.add_option("-f","--from", type= "int", dest="_from")
    parser.add_option("-t","--to", type= "int", dest="to")
    parser.add_option("-s",action= "store_true", default=False, dest="score")
    (options, args) = parser.parse_args()

    _, eigenvectors, header = ProdyNMDParser.read(options.input1)
    header1, eigenvectors1 = get_CA_modes(header, eigenvectors)
    _, eigenvectors, header = ProdyNMDParser.read(options.input2)
    header2, eigenvectors2 = get_CA_modes(header, eigenvectors)
    
#    numpy.savetxt("arrays.txt", numpy.array([eigenvectors1[0], eigenvectors2[0]]))
    
    assert options.to >= options._from,  "[ERROR] 'from' value is bigger than 'to'. "
    assert eigenvectors1.shape == eigenvectors2.shape,  "[ERROR] eigenvectors must share the same shape(1: %s,2: %s)."%(str(eigenvectors1.shape),str(eigenvectors2.shape))
    
    scores = []
    for i in range(len(eigenvectors1)):
        if i >= options._from and i<= options.to:
            mode1 = eigenvectors1[i]
            mode2 = eigenvectors2[i]
            
            angles = measure.calculate_mode_angles(mode1, mode2)
#            numpy.savetxt("angles.txt", angles)
    
            magnitudes1 = measure.calculate_mode_magnitudes(mode1)
            magnitudes2 = measure.calculate_mode_magnitudes(mode2)
            if not options.score:
                plt.title('----')
                ax = plt.subplot2grid((4, 1), (0,0), rowspan = 2)
                plt.plot(angles ,magnitudes1,  'r.')
                ax.set_xlabel("Angles")
                ax.set_ylabel("Magnitude")
                
                plt.subplot2grid((4, 1), (2,0))
                plt.plot(range(len(magnitudes1)), magnitudes1, 'y.-')
                plt.ylabel('Magnitude')
                plt.xlabel('Residue ID')
                
                plt.subplot2grid((4, 1), (3,0))
                plt.plot(range(len(magnitudes2)), magnitudes2, 'g.-')
                plt.ylabel('Magnitude')
                plt.xlabel('Residue ID')
                
                plt.show()
            else:
                """
                Perform the scoring of this guys (mean error of the normalized magnitudes and
                mean angle error weighted by its correct magnitude value )
                """
                norm_magnitudes1 = normalize(magnitudes1)
                norm_magnitudes2 = normalize(magnitudes2)
                ang_error = (norm_magnitudes1*angles).sum() / len(angles)
                scores.append(ang_error)
                
    if options.score:
        if options.output is not None:
            numpy.savetxt(options.output, scores)
        else:
            print "Scores:"
            for s in scores:
                print s
                