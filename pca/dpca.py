"""
Created on 18/05/2015

@author: user
"""
import numpy
import math

def circular_mean(angles):
    """
    Calculates the circular mean as described by Altis 2008 (thesis).
    """
    c = numpy.cos(angles).sum()
    s = numpy.sin(angles).sum()
    m = 0
    if c > 0.0:
        m = math.atan2(s,c)
    elif c < 0.0:
        m = math.atan2(s,c) + math.pi
    else: # c == 0
        m = math.copysign((math.pi / 2.0), s)
    return m

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
    
def covariance_matrix(data, d_mean):
    """
    Calculates the 'variance-covariance' matrix for a data set of the type:
    [[observation1],
    [observation 2],
    ...
    [observation n]]
    """
    num_variables = len(data)
    X = data - d_mean
    return numpy.dot(X.T, X.conj()) / (num_variables-1)
            
class dPCA:
    def __init__(self):
        pass
    
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
        new_angles = []
        
        # Shift all angles in order to avoid having them in boundaries
        for angle_observations in angles.T:
            phi_0 =  min_density_angle(angle_observations)
            new_angles.append(shift_all_above(angle_observations, phi_0))
        
        # Calculate the circular mean for each of the dihedrals
        mean_vector = [circular_mean(a) for a in angle_observations.T]
        
        # Calculate the covariance matrix
        cov = covariance_matrix(angles, mean_vector)
        
        # Calculate the components
        w, v = numpy.linalg.eig(cov)
        
        return w, v.T
        
        
            