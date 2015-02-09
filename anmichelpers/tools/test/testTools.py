'''
Created on 9/2/2015

@author: victor
'''
import unittest
import numpy
from anmichelpers.tools.tools import norm, ensure_modes_layout

class Test(unittest.TestCase):


    def test_norm(self):
        v1 = numpy.array([1.,2.,3.])
        v2 = numpy.array([[1.,2.,3.],
                          [4.,5.,6.],
                          [7.,8.,9.]])
        
        self.assertAlmostEqual(3.74165738677,
                               norm(v1))
        
        numpy.testing.assert_array_almost_equal(norm(v2),
                                                [3.7416573867739413, 8.7749643873921226, 13.928388277184119])
        
    def test_ensure_mode_layout(self):
        twod_eigen = numpy.array([[1., 2., 3., 4., 5., 6., 7., 8., 9.],
                                  [1., 2., 3., 4., 5., 6., 7., 8., 9.]])
        threed_eigen = numpy.array([
                                    [[1., 2., 3.],[ 4., 5., 6.], [7., 8., 9.]],
                                    ])
        twod_eigen = ensure_modes_layout(twod_eigen)
        threed_eigen = ensure_modes_layout(threed_eigen)
        self.assertTupleEqual((2, 3, 3), twod_eigen.shape)
        self.assertTupleEqual((1, 3, 3), threed_eigen.shape)
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_norm']
    unittest.main()