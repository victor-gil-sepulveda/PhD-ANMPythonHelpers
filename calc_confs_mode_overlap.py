'''
Created on 26/10/2015

@author: victor
'''

# one time script!

import os
from prody.proteins.pdbfile import parsePDB
from prody.dynamics.anm import ANM
import numpy

if __name__ == '__main__':
    distances = [ 31, 29, 24, 17, 13, 07, 04 ]
    base_folder = "/home/victor/Escritorio/1AKE_dyn"
    
    # Calculate eigenvectors
    all_eigenvectors = {}
    NUM_MODES = 8
    for distance in distances:
        pdb_file = "%02d.pdb"%distance
        pdb_path = os.path.join(base_folder, pdb_file)
        print pdb_path
        pdb_struct = parsePDB(pdb_path)
        pdb_struct_ca = pdb_struct.ca
        pdb_struct_ca_anm = ANM(pdb_file) 
        pdb_struct_ca_anm.buildHessian(pdb_struct_ca) #cutoff 15
        pdb_struct_ca_anm.calcModes(n_modes=NUM_MODES) #cutoff 15
        
        eigenvectors = []
        for i in range(NUM_MODES):
            mode = pdb_struct_ca_anm[i]
            eigenvectors.append(mode.getEigvec().round(3))

        all_eigenvectors[pdb_file] = numpy.array(eigenvectors)
    
    # Calculate overlaps of first modes vs the others
    overlaps = []
    reference = "31.pdb"
    for mode_i in range(NUM_MODES):
        for distance in distances:
            pdb_file = "%02d.pdb"%distance
            cumulative_overlap(all_eigenvectors[reference][mode_i], all_eigenvectors[pdb_file])
        
    
            