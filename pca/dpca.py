"""
Created on 18/05/2015

@author: user
"""
import numpy
import math

def min_density_angle(angles, numbins = 20):
    """
    Returns the angle of less density as the less populated range in the  
    angle histogram (population != 0) 
    """
    his, edges  = numpy.histogram(angles, bins = numbins)
    
    # Get the minimum value != 0
    min_val = (len(angles)+1,-1)
    for i in range(len(his)):
        if his[i] != 0:
            min_val = min(min_val, (his[i],i))
    return (edges[min_val[1]]+edges[min_val[1]])/2.

def shift_all_above(angles, min_dens_val):
    """
    Shifts -360 degrees all the angles above the value 'min_dens_val' where above means 
    'in the range  [min_dens_val, pi] given that the values in the angles array are in the
    range [-pi, pi].
    """
    new_angles = []
    for angle in angles:
        if angle > min_dens_val:
            new_angles.append(angle - (2*math.pi))
        else:
            new_angles.append(angle)
    return new_angles

def circular_mean_and_variance(angles):
    """
    Fisher 1993, Mardia & Jupp 2000
    """
    n = len(angles)
    Cm = numpy.mean(numpy.cos(angles))
    Sm = numpy.mean(numpy.sin(angles))
    Rp = math.sqrt(Cm*Cm + Sm*Sm)
    Rpm= Rp /n;
    
    m = 0
    if Cm > 0.0 and Sm > 0.0 :
        m = math.atan2(Sm, Cm)
    elif Cm < 0.0:
        m = math.atan2(Sm, Cm) + math.pi
    elif Cm < 0.0 and Sm < 0.0:
        m = math.atan2(Sm, Cm) + 2*math.pi
    else: # c == 0
        m = math.copysign((math.pi / 2.0), Sm)
        
    return m, 1-Rpm

def two_angle_circular_correlation_coef(angles1, angles2, mean1, mean2):
    """
    Circular correlation measure. SenGupta 2001
    """
    centered_a = angles1-mean1
    centered_b = angles2-mean2
    sin_centered_a = numpy.sin(centered_a)
    sin_centered_b = numpy.sin(centered_b)
    sin2_a = sin_centered_a*sin_centered_a
    sin2_b = sin_centered_b*sin_centered_b
    
    return numpy.dot(sin_centered_a, sin_centered_b) / math.sqrt(numpy.dot(sin2_a, sin2_b))
   
def calculate_var_covar_matrix(angles, mean_a, var_a):
    """
    Layout:
    [[angles for struct 1],
        [angles for struct 2],
        ...
        [angles for struct n]]
    """
    num_angles = len(angles[0])
    cov_matrix = numpy.zeros((num_angles,)*2)
    angles_per_residue = angles.T
    for i in range(num_angles-1):
        cov_matrix[i][i] = var_a[i]
        for j in range(i+1, num_angles):
            angles_i = angles_per_residue[i]
            angles_j = angles_per_residue[j]
#            print angles_i
#            print angles_j
            cov_matrix[i][j] = two_angle_circular_correlation_coef(angles_i,
                                                                   angles_j,
                                                                   mean_a[i],
                                                                   mean_a[j])
            cov_matrix[j][i] = cov_matrix[i][j]
    return cov_matrix

def to_pi_mpi_range(angle):
    """
    Puts an angle in the -pi, pi range
    """
    if angle > math.pi:
        return angle - math.pi
    elif angle < - math.pi:
        return math.pi + angle
    else:
        return angle

class dPCA:
    def __init__(self):
        pass
    
    @classmethod
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
    
    @classmethod
    def calc(cls, angles):
        """
        Given a matrix of angles for different structures with this layout:
        [[angles for struct 1],
        [angles for struct 2],
        ...
        [angles for struct n]]
        
        Calculates the direct dPCA of the dihedrals.
        """
        mean_angles = []
        var_angles = []
        # Shift all angles in order to avoid having them in boundaries
        all_shifted_angles = []
        for angle_observations in angles.T:
            
            shifted_angles = cls.center_angles_in_most_dense_peak(angle_observations)
            
            all_shifted_angles.append(shifted_angles)
            mean_a, var_a = circular_mean_and_variance(shifted_angles)
            mean_angles.append(mean_a)
            var_angles.append(var_a)
        all_shifted_angles = numpy.array(all_shifted_angles).T
        var_angles = numpy.array(var_angles)
        mean_angles = numpy.array(mean_angles)
        
        # Calculate the covariance matrix
        cov = calculate_var_covar_matrix(all_shifted_angles, 
                                         mean_angles, 
                                         var_angles)
        
        # Calculate the components
        # TODO! Order is not guaranteed.  Order the eigenvectors by eigenvalue!
        w, v = numpy.linalg.eig(cov)
        
        return w[0:10].astype(float), v[:, 0:10].T.astype(float)
        
        
            