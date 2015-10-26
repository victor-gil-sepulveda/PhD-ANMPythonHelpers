'''
Created on 25/9/2015

@author: victor
'''
import sys
import numpy
from prody.proteins.pdbfile import parsePDB
from prody.measure.measure import calcPhi, calcPsi

if __name__ == '__main__':
    trajectory_path = sys.argv[1]
    
    structure = parsePDB(trajectory_path)
    all_angles = []
    for i, coords in enumerate(structure.getCoordsets()):
        angles = []
        structure.setCoords(coords)
        hv = structure.getHierView()
        for residue in hv.iterResidues():
            try:
                angles.append(calcPhi(residue, radian=True))
            except:
                angles.append(0.)
            try:
                angles.append(calcPsi(residue, radian=True))
            except:
                angles.append(0.)
        all_angles.append(angles)
    all_angles = numpy.array(all_angles)
    numpy.savetxt(sys.argv[1]+".ang", all_angles, fmt = "%.4f")
