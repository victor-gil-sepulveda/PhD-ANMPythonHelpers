"""
Created on Feb 5, 2015

@author: victor
"""
import unittest
import os
import numpy
import anmichelpers.parsers.test
from anmichelpers.parsers.imods import ImodServerFilesParser
from anmichelpers.parsers.bahar import BaharServerFilesParser

class TestParsers(unittest.TestCase):

    def test_imods(self):
        expected_eigenvalues = [5.9864473576e-06, 9.6676034548e-06]
        expected_eigenvectors = [
                                 [-0.0015498569747, -0.00019775087154000001, -0.002623811054, -0.0015369823802999999, -0.00023667061667, -0.0026169987178000001, -0.0015035583555, -0.0025828151242000001, 0.00017850569692000001, -0.0060579396906999998, -0.0026160250078999998, -0.0012606014619000001, -0.0059246699867999998, -0.0027030145931000002, 0.0028110840724000002, -0.0068931796958000002], 
                                 [-4.7378807072000002e-08, 1.8344767783999999e-07, 1.7676846446000001e-07, -1.0385215654e-07, 1.3044119216e-07, -6.0814995075000004e-08, -1.6662504585999999e-07, 3.0547344434999999e-07, -5.3739931160999996e-07, 2.2262028449e-07, -3.8639963417e-07, 2.1906448182e-07, 2.0558392479999999e-08, 1.3043272768000001e-07, -2.0228736609e-07, -2.5653969897999998e-07]
                                 ]

        file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","imods1.evec")
        eigenvalues, eigenvectors = ImodServerFilesParser.read(file_path)
        numpy.testing.assert_array_almost_equal(expected_eigenvalues, eigenvalues, 12)
        numpy.testing.assert_array_almost_equal(expected_eigenvectors, eigenvectors, 12)
    
    def test_bahars(self):
        expected_eigenvalues = [-2.33561e-07, -1.19815e-07, -7.76845e-08, 8.21996e-09, 4.55134e-08, 2.5498e-07, 0.00356931, 0.00543829, 0.0269361, 0.0501762, 0.0844318, 0.133118, 0.158139, 0.218802, 0.260598, 0.310018, 0.346588, 0.474664, 0.633895, 0.711878, 0.948204, 1.00375, 1.07904, 1.14168, 1.3709, 1.39127, 1.56822, 1.62059, 1.64899, 1.76275, 1.84827, 1.86695, 2.03717, 2.14711, 2.26446, 2.38432]
        expected_eigenvector_0 = [-0.038114, -0.021852, -0.04286, -0.03734, -0.019061, -0.031715, -0.047092, -0.020868, -0.032067, -0.042022, -0.022238, -0.041522, -0.033184, -0.019075, -0.034744, -0.037456, -0.017095, -0.023331, -0.047211, -0.017938, -0.019741, -0.05187, -0.016575, -0.010641, -0.06111, -0.017474, -0.007758, -0.052899, -0.015021, -0.004071]

        evec_file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","bahars.slwevs")
        eval_file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","bahars.eigvals")
        eigenvalues, eigenvectors = BaharServerFilesParser.read(eval_file_path, evec_file_path)
        
        numpy.testing.assert_array_almost_equal(expected_eigenvalues, eigenvalues, 12)
        numpy.testing.assert_array_almost_equal(expected_eigenvector_0, eigenvectors[0], 12)
        
    def test_prodynmd(self):
        file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","prodynmd.nmd")
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()