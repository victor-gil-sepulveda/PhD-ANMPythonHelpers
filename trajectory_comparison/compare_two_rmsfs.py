'''
Created on Aug 21, 2015

@author: victor
'''

import sys
import math
import numpy
import matplotlib.pyplot as plt
import os

def rms(one, two):
    return math.sqrt( numpy.dot(one-two, one-two)/len(one))

if __name__ == '__main__':
    
    rmsf1 = numpy.loadtxt(sys.argv[1])
    rmsf2 = numpy.loadtxt(sys.argv[2])

    common_length = min(len(rmsf1),len(rmsf2))
    rmsf1 = rmsf1[0:common_length]
    rmsf2 = rmsf2[0:common_length]

    print "-----------------------"
    print " RMS: ",rms(rmsf1, rmsf2)
    print "-----------------------"

    plt.plot(rmsf1, label = os.path.basename(sys.argv[1]))
    plt.plot(rmsf2, label = os.path.basename(sys.argv[2]))
    plt.legend(loc=9, ncol = 4)
    plt.show()
        