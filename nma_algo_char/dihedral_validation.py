'''
Created on Nov 26, 2015

@author: victor
'''
from prody.proteins.pdbfile import parsePDB
import sys
from calculate_dihedrals import get_dihedrals_for_conformation
import numpy
import math

if __name__ == '__main__':
    pdb = parsePDB(sys.argv[1])
    
    coordsets = pdb.getCoordsets()
    for i in range(len( coordsets)-1):
        pdb.setCoords(coordsets[i])
        dihedrals_i = get_dihedrals_for_conformation(pdb)
        pdb.setCoords(coordsets[i+1])
        dihedrals_ip1 = get_dihedrals_for_conformation(pdb)
        # sin(\alpha - \beta) = sin \alpha cos \beta - cos \alpha sin \beta.
        sub = numpy.arcsin(numpy.sin(dihedrals_ip1)*numpy.cos(dihedrals_i)-
                           numpy.cos(dihedrals_ip1)*numpy.sin(dihedrals_i))
        print numpy.max(numpy.abs(sub))
        