'''
Created on 25/9/2015

@author: victor
'''
import sys
from prody.proteins.pdbfile import parsePDB
from prody.measure.measure import calcPhi, calcPsi
import pca.dpca
import numpy
import math



def center_angles_in_most_dense_peak(cls, angles, hist_bin_size = 0.1):
    """
    Improved angle preconditioning
    """
    
    hist, boundaries = numpy.histogram(angles, 
                                       bins=int((math.pi*2)/ hist_bin_size),
                                       range = (-math.pi, math.pi))

    print "***********************" 
    print "***********************"
    print "***********************"       
    print hist
    print "***********************"
    
    # Detect base value
    min_val = numpy.min(hist)
    
    # Shift histogram down
    shift_hist = hist - min_val
    # Detect zeros
    min_indices = numpy.where(shift_hist == 0)[0]    
    
    # Remove 0s that have 0s around
    to_be_removed = []
    for i in min_indices:
        left = (i-1) % len(hist)
        right = (i+1) %len(hist)
        if shift_hist[left] == 0 and shift_hist[right]==0:
            to_be_removed.append(i)
            
    # Suavizar hasta que solo quede un pico
    
    
    min_indices = numpy.array(sorted(list(set(min_indices)-set(to_be_removed))))
    
    # We circularly calculate the densities between zeros
    densities = []
    for  i in range(len(min_indices)-1):
        d = numpy.sum(hist[min_indices[i]:min_indices[i+1]+1])
        densities.append((d, min_indices[i], min_indices[i+1]))
    
    # Last can cross boundaries
    d = numpy.sum(hist[min_indices[-1]:len(min_indices)])+ numpy.sum(hist[0:min_indices[0]])
    densities.append((d, min_indices[-1], min_indices[0]))
    
    # Pick the biggest zone and center the angles there
    most_dense_zone =  max(densities)
    
    # Then pick the peak :D
    left, right = most_dense_zone[1], most_dense_zone[2]+1
    if left > right:
        left_peak_value = numpy.max(hist[0:right])
        right_peak_value = numpy.max(hist[left:len(hist)])
        if left_peak_value > right_peak_value:
            peak_index = hist[0:right].tolist().index(left_peak_value)
        else:
            peak_index = hist[left:len(hist)].tolist().index(right_peak_value)
    else:
        peak_max = numpy.max(hist[left:right])
        peak_index = left + hist[left:right].tolist().index(peak_max)
        
    peak_value = (boundaries[peak_index] + boundaries[peak_index+1]) / 2.
    # Now center everything in that value
    shifted_angles = angles - peak_value
    # And put the angles that are out of boundaries in the same -pi, pi range
    corrected_angles = []
    for angle in shifted_angles:
        corrected_angles.append(to_pi_mpi_range(angle))
    
    hist, _ = numpy.histogram(corrected_angles, 
                              bins=int((math.pi*2)/ hist_bin_size),
                              range = (-math.pi, math.pi))
    print hist
    print "**********************"
    return numpy.array(corrected_angles)


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
                angles.append(0)
            try:
                angles.append(calcPsi(residue, radian=True))
            except:
                angles.append(0)
        all_angles.append(angles)
    all_angles = numpy.array(all_angles)
    numpy.savetxt(sys.argv[2]+".ang", all_angles)
    
    # Shift all angles in order to avoid having them in boundaries
    var_angles = []
    all_shifted_angles = []
    for angle_observations in all_angles.T:
        shifted_angles = pca.dpca.dPCA.center_angles_in_most_dense_peak(angle_observations)
        all_shifted_angles.append(shifted_angles)
        _, var_a = pca.dpca.circular_mean_and_variance(shifted_angles)
        var_angles.append(var_a)
    
    numpy.savetxt(sys.argv[2]+".ang.sh", all_angles)
    numpy.savetxt(sys.argv[2], var_angles)