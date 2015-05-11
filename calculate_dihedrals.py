"""
Created on 07/05/2015

@author: user
"""

import numpy
from optparse import OptionParser
from prody.measure.measure import calcPsi, calcPhi, calcOmega
from prody.proteins.pdbfile import parsePDB

def get_dihedrals_for_conformation(conformation):
    dihedral_angles = []
    for residue in conformation.iterResidues():
        if residue.getResname() != "PRO":
            try:
                dihedral_angles.append(calcPhi(residue, radian=True))
            except ValueError:
                dihedral_angles.append(0)
        try:
            dihedral_angles.append(calcPsi(residue, radian=True))
        except ValueError:
            dihedral_angles.append(0)
    return numpy.array(dihedral_angles[1:-1])

def get_omegas_for_conformation(conformation):
    omega_angles = []
    for residue in conformation.iterResidues():
        try:
            omega_angles.append((residue.getResname(),calcOmega(residue, radian=True)))
        except ValueError:
            omega_angles.append(("--",0))
    return numpy.array(omega_angles)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-o", dest="output")
    parser.add_option("-w", action= "store_true", dest="omega")
    (options, args) = parser.parse_args()
    
    pdb = parsePDB(options.input)
    
    if not options.omega:
        angles = get_dihedrals_for_conformation(pdb)
        numpy.savetxt(options.output, angles, delimiter = "\n")
    else:
        angles = get_omegas_for_conformation(pdb)
        open(options.output,"w").write("\n".join([str(omega) for omega in angles]))
    
    
#    
#open = numpy.loadtxt("5XHK_helix.dih")
#closed = numpy.loadtxt("9WVG_helix.dih")
#numpy.savetxt("diff.dih", open-closed)

       
    