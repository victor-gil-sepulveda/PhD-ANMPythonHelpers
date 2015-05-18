'''
Created on 18/05/2015

@author: user
'''
import unittest
import numpy
from pca.dpca import covariance_matrix, min_density_angle, shift_all_above
#import pylab

class TestdPCA(unittest.TestCase):


    def test_cov_matrix(self):
        x_mean = numpy.array([4.1, 2.08, 0.604])
        # 5 observations of 3 vars
        X = numpy.array([   [ 4. , 2. , 0.6 ],
                            [ 4.2, 2.1, 0.59],
                            [ 3.9, 2.,  0.58],
                            [ 4.3, 2.1, 0.62],
                            [ 4.1, 2.2, 0.63]])
        expected_result = [[ 0.025,    0.0075,   0.00175],
                           [ 0.0075,   0.007,    0.00135],
                           [ 0.00175,  0.00135,  0.00043]]
                                    
        cov = covariance_matrix(X, x_mean)
        
        numpy.testing.assert_array_almost_equal(expected_result, cov, 8)
    
    def test_min_density(self):
        angles = [-0.8, -0.75, -0.7, -0.65, -0.6, -0.5, -0.5, -0.4, -0.4, -0.4, -0.3, 0, 0.1, 0.2, 0.2, 0.2, 0.3, 0.4]  
        print 
#        pylab.hist(angles)
#        pylab.show()
        self.assertAlmostEqual(-0.32, min_density_angle(angles, numbins = 10), 10)
    
    def test_shift_all_above(self):
        angles = [-0.8, -0.75, -0.7, -0.65, -0.6, -0.5, -0.5, -0.4, -0.4, -0.4, -0.3, 0, 0.1, 0.2, 0.2, 0.2, 0.3, 0.4]  
        expected = [-0.8, -0.75, -0.7, -0.65, -0.6, -0.5, -0.5, -0.4, -0.4, -0.4, -6.583, -6.283, -6.183, -6.083, -6.083, -6.083, -5.983, -5.883]

        min_dens = -0.32
        numpy.testing.assert_array_almost_equal(expected, shift_all_above(angles, min_dens), 3)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_cov_matrix']
    unittest.main()