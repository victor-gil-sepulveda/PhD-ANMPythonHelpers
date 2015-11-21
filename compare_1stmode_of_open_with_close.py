"""
Created on Oct 19, 2015

@author: victor
"""

from anmichelpers.parsers.pronmd import ProdyNMDParser
from anmichelpers.comparison.comparison import cumulative_overlap, overlap
from anmichelpers.tools.tools import norm
import numpy

if __name__ == '__main__':
    
    OPEN_NMD = "/home/victor/Desktop/Escritorio/OpenCloseANM/open_anm.nmd"
    CLOSED_NMD = "/home/victor/Desktop/Escritorio/OpenCloseANM/close_anm.nmd"
    
    o_evals, o_evecs, o_header = ProdyNMDParser.read(OPEN_NMD)
    c_evals, c_evecs, c_header = ProdyNMDParser.read(CLOSED_NMD)
    
    
    print "Number of modes Open:", len(o_evecs), " Closed: ",len(c_evecs)
    print "Avg. mode norms. Open: ", numpy.mean(norm(o_evecs)), " Closed: ", numpy.mean(norm(c_evecs))
    
    print "1st open with others", cumulative_overlap(o_evecs[0],o_evecs[1:])
    print "1st open with all closed", cumulative_overlap(o_evecs[0],c_evecs)
    print "1st open with first 3", cumulative_overlap(o_evecs[0],c_evecs[0:3])
    print "1st open with first 6", cumulative_overlap(o_evecs[0],c_evecs[0:6])
    print "Single overlaps of first mode:"
    for i in range(len(c_evecs)):
        print "\t ",i, overlap(o_evecs[0],c_evecs[i])