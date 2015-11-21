'''
Created on 26/10/2015

@author: victor
'''

# one time script!

import os
from prody.proteins.pdbfile import parsePDB
from prody.dynamics.anm import ANM
import numpy
from anmichelpers.comparison.comparison import cumulative_overlap, overlap
from anmichelpers.tools.tools import norm

if __name__ == '__main__':
    distances = [ 33, 29, 24, 17, 13, 07, 04 ]
    base_folder = "/home/victor/Desktop/1AKE_dyn/"
    reference_open = "33.pdb"
    reference_closed = "04.pdb"
    
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

    others  = list(sorted(all_eigenvectors.keys(), reverse=True))
    
    
    # Calculate transition vectors
    open_struct = parsePDB(os.path.join(base_folder,reference_open)).ca
    closed_struct = parsePDB(os.path.join(base_folder,reference_closed)).ca
    trans_vectors = (open_struct.getCoords() - closed_struct.getCoords()).flatten()

    # Calculate maximum single overlap with translation
    most_overlapping_modes_with_translation = []
    for pdb_file in others:
        overlaps = []
        for j in range(NUM_MODES):
            o = overlap(trans_vectors, all_eigenvectors[pdb_file][j])
            overlaps.append((o,j))
        overlaps.sort(reverse=True)
        most_overlapping_modes_with_translation.append(overlaps[0])
    most_overlapping_modes_with_translation = numpy.array(most_overlapping_modes_with_translation)
    
    # Calculate maximum single overlap with first
    most_overlapping_modes_with_reference = []
    for pdb_file in others:
        overlaps = []
        for j in range(NUM_MODES):
            o = overlap(all_eigenvectors[reference_open][0], all_eigenvectors[pdb_file][j])
            overlaps.append((o,j))
        overlaps.sort(reverse=True)
        most_overlapping_modes_with_reference.append(overlaps[0])
    most_overlapping_modes_with_reference = numpy.array(most_overlapping_modes_with_reference)
    
    # Calculate cumulative overlaps of first reference modes vs the others
    results = [[i] for i in range(len(others))]
    for mode_i in range(NUM_MODES):
        for j, other in enumerate(others):
            c_ov = cumulative_overlap(all_eigenvectors[reference_open][mode_i], all_eigenvectors[other])
            results[j].append(c_ov)
    
    # Calculate overlaps of translation vs all modes
    results_cum_trans = [[i] for i in range(len(others))]
    for j, other in enumerate(others):
        c_ov = cumulative_overlap(trans_vectors, all_eigenvectors[other])
        results_cum_trans[j].append(c_ov)
    
    numpy.savetxt("max_sing_ov_with_trans.txt", numpy.array(most_overlapping_modes_with_translation),fmt="%.3f")
    numpy.savetxt("max_sing_ov_with_first_ref.txt", numpy.array(most_overlapping_modes_with_reference),fmt="%.3f")
    numpy.savetxt("cum_overlap_with_trans.txt", numpy.array(results_cum_trans),fmt="%.3f")
    numpy.savetxt("cum_overlap_with_first_ref.txt", numpy.array(results),fmt="%.3f")
    numpy.savetxt("distances.txt", numpy.array(distances),fmt="%.1f")
     
            