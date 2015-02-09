'''
Created on 9/2/2015

@author: victor
'''
import unittest
from anmichelpers.tools.measure import mean_square_fluctuations
import numpy


class TestMeasure(unittest.TestCase):


    def test_mean_square_fluctuations(self):
        number_of_nodes = 3
        number_of_modes = 9
        eigenvalues = numpy.array([0., 0., 0., 0., 0., 0., 1., 2., 3.])
        significative = 6
        eigenvectors = numpy.array([
                        [1.]*number_of_modes,
                        [2.]*number_of_modes,
                        [3.]*number_of_modes,
                        [4.]*number_of_modes,
                        [5.]*number_of_modes,
                        [6.]*number_of_modes,
                        range(1,10),
                        range(2,11),
                        range(3,12)
                        ])
        numpy.testing.assert_array_almost_equal([45.16666667, 181.66666667, 417.16666667],
                              mean_square_fluctuations(eigenvalues[significative:], 
                                       eigenvectors[significative:], 
                                       number_of_nodes),
                                8)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()