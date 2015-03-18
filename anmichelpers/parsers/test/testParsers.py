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
from anmichelpers.parsers.pronmd import ProdyNMDParser
from anmichelpers.writers.pronmd import ProdyNMDWriter
import prody 

class TestParsers(unittest.TestCase):

#     def test_imods(self):
#         expected_eigenvalues = [5.9864473576e-06, 9.6676034548e-06]
#         expected_eigenvectors = [
#                                  [-0.0015498569747, -0.00019775087154000001, -0.002623811054, -0.0015369823802999999, -0.00023667061667, -0.0026169987178000001, -0.0015035583555, -0.0025828151242000001, 0.00017850569692000001, -0.0060579396906999998, -0.0026160250078999998, -0.0012606014619000001, -0.0059246699867999998, -0.0027030145931000002, 0.0028110840724000002, -0.0068931796958000002], 
#                                  [-4.7378807072000002e-08, 1.8344767783999999e-07, 1.7676846446000001e-07, -1.0385215654e-07, 1.3044119216e-07, -6.0814995075000004e-08, -1.6662504585999999e-07, 3.0547344434999999e-07, -5.3739931160999996e-07, 2.2262028449e-07, -3.8639963417e-07, 2.1906448182e-07, 2.0558392479999999e-08, 1.3043272768000001e-07, -2.0228736609e-07, -2.5653969897999998e-07]
#                                  ]
# 
#         file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","imods1.evec")
#         eigenvalues, eigenvectors = ImodServerFilesParser.read(file_path)
#         numpy.testing.assert_array_almost_equal(expected_eigenvalues, eigenvalues, 12)
#         numpy.testing.assert_array_almost_equal(expected_eigenvectors, eigenvectors, 12)
#     
#     def test_bahars(self):
#         expected_eigenvalues = [-2.33561e-07, -1.19815e-07, -7.76845e-08, 8.21996e-09, 4.55134e-08, 2.5498e-07, 0.00356931, 0.00543829, 0.0269361, 0.0501762, 0.0844318, 0.133118, 0.158139, 0.218802, 0.260598, 0.310018, 0.346588, 0.474664, 0.633895, 0.711878, 0.948204, 1.00375, 1.07904, 1.14168, 1.3709, 1.39127, 1.56822, 1.62059, 1.64899, 1.76275, 1.84827, 1.86695, 2.03717, 2.14711, 2.26446, 2.38432]
#         expected_eigenvector_0 = [-0.038114, -0.021852, -0.04286, -0.03734, -0.019061, -0.031715, -0.047092, -0.020868, -0.032067, -0.042022, -0.022238, -0.041522, -0.033184, -0.019075, -0.034744, -0.037456, -0.017095, -0.023331, -0.047211, -0.017938, -0.019741, -0.05187, -0.016575, -0.010641, -0.06111, -0.017474, -0.007758, -0.052899, -0.015021, -0.004071]
# 
#         evec_file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","bahars.slwevs")
#         eval_file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","bahars.eigvals")
#         eigenvalues, eigenvectors = BaharServerFilesParser.read(eval_file_path, evec_file_path)
#         
#         numpy.testing.assert_array_almost_equal(expected_eigenvalues, eigenvalues, 12)
#         numpy.testing.assert_array_almost_equal(expected_eigenvector_0, eigenvectors[0], 12)
#         
    def test_prodynmd(self):
        file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","prodynmd.nmd")
        expected_eigenvalues = [ 2.36446078,  1.72973341,  1.70067417,  1.12434682,  1.03043392, 0.9941346] 
        expected_eigenvector_0 = [0.039, 0.009, 0.058, 0.038, -0.011, 0.052, 0.043, -0.027, 0.058, 0.052, -0.048, 0.062, 0.049, -0.056, 0.042, 0.057, -0.077, 0.036, 0.055, -0.088, 0.016, 0.059, -0.1, 0.004, 0.057, -0.106, -0.017, 0.055, -0.105, -0.024, 0.051, -0.11, -0.046, 0.052, -0.11, -0.048, 0.039, -0.09, -0.037, 0.042, -0.086, -0.021, 0.04, -0.074, -0.007, 0.041, -0.065, 0.009, 0.039, -0.052, 0.021, 0.042, -0.045, 0.04, 0.053, -0.056, 0.056, 0.045, -0.042, 0.057, 0.042, -0.049, 0.038, 0.052, -0.066, 0.04, 0.056, -0.08, 0.027, 0.047, -0.076, 0.007, 0.046, -0.085, -0.011, 0.039, -0.083, -0.032, 0.027, -0.072, -0.045, 0.019, -0.07, -0.065, 0.016, -0.072, -0.067, 0.004, -0.067, -0.081, -0.011, -0.062, -0.063, -0.026, -0.043, -0.032, -0.002, -0.047, -0.046, 0.015, -0.059, -0.049, 0.018, -0.053, -0.035, 0.03, -0.061, -0.02, 0.03, -0.055, -0.006, 0.036, -0.054, 0.016, 0.038, -0.047, 0.038, 0.033, -0.033, 0.045, 0.04, -0.036, 0.063, 0.034, -0.019, 0.066, 0.025, -0.009, 0.061, 0.029, -0.018, 0.066, 0.023, -0.015, 0.047, 0.024, -0.03, 0.029, 0.02, -0.027, 0.014, 0.017, -0.031, -0.005, 0.02, -0.036, -0.007, 0.011, -0.036, -0.024, 0.013, -0.042, -0.028, -0.002, -0.035, -0.035, -0.001, -0.049, -0.049, -0.01, -0.041, -0.049, -0.018, -0.023, -0.044, -0.025, -0.018, -0.051, -0.032, -0.029, -0.063, -0.044, -0.017, -0.071, -0.046, 0.001, -0.069, -0.044, -0.003, -0.07, -0.038, -0.014, -0.056, -0.033, -0.002, -0.05, -0.035, 0.011, -0.051, -0.032, 0.003, -0.05, -0.024, -0.001, -0.038, -0.024, 0.015, -0.034, -0.024, 0.019, -0.033, -0.02, 0.012, -0.028, -0.013, 0.014, -0.02, -0.016, 0.028, -0.019, -0.014, 0.025, -0.017, -0.01, 0.02, -0.014, -0.007, 0.026, -0.007, -0.009, 0.032, -0.01, -0.009, 0.027, -0.011, -0.006, 0.025, -0.007, -0.01, 0.023, -0.015, -0.011, 0.019, -0.02, -0.011, 0.019, -0.021, -0.01, 0.021, -0.017, -0.008, 0.022, -0.015, -0.004, 0.018, -0.006, 0.003, 0.012, 0.002, 0.013, 0.005, 0.015, 0.013, -0.0, 0.017, 0.01, -0.003, 0.006, 0.016, -0.013, 0.009, 0.017, -0.022, 0.003, 0.028, -0.036, 0.012, 0.019, -0.027, 0.004, 0.023, -0.018, 0.017, 0.01, -0.005, 0.005, -0.002, 0.003, -0.01, -0.003, -0.008, -0.017, 0.005, -0.024, -0.015, 0.002, -0.017, -0.016, -0.0, -0.026, -0.026, 0.006, -0.033, -0.021, 0.004, -0.021, -0.017, 0.012, -0.026, -0.011, 0.006, -0.014, -0.008, 0.012, -0.015, 0.002, 0.007, -0.003, -0.0, 0.007, -0.003, 0.004, 0.01, -0.01, -0.004, 0.002, 0.015, -0.018, 0.0, 0.024, -0.022, -0.002, 0.039, -0.031, -0.005, 0.036, -0.024, -0.008, 0.035, -0.022, -0.006, 0.042, -0.017, -0.005, 0.046, -0.024, -0.009, 0.04, -0.029, -0.007, 0.037, -0.023, -0.005, 0.047, -0.019, -0.007, 0.05, -0.028, -0.008, 0.041, -0.026, -0.009, 0.03, -0.026, -0.012, 0.027, -0.033, -0.013, 0.019, -0.039, -0.012, 0.009, -0.036, -0.014, 0.009, -0.042, -0.014, 0.018, -0.037, -0.011, 0.014, -0.029, -0.012, 0.005, -0.033, -0.012, 0.011, -0.035, -0.011, 0.017, -0.029, -0.01, 0.011, -0.023, -0.01, 0.005, -0.028, -0.01, 0.012, -0.028, -0.009, 0.014, -0.02, -0.008, 0.006, -0.019, -0.008, 0.007, -0.023, -0.009, 0.014, -0.02, -0.008, 0.011, -0.014, -0.007, 0.005, -0.015, -0.007, 0.012, -0.016, -0.008, 0.016, -0.014, -0.008, 0.01, -0.01, -0.006, 0.01, -0.012, -0.007, 0.021, -0.016, -0.008, 0.016, -0.014, -0.012, 0.019, -0.014, -0.01, 0.016, -0.009, -0.009, 0.02, -0.007, -0.006, 0.022, 0.002, -0.009, 0.027, 0.002, -0.008, 0.024, -0.004, -0.008, 0.029, -0.005, -0.007, 0.033, -0.012, -0.006, 0.037, -0.012, -0.007, 0.028, -0.013, -0.007, 0.028, -0.019, -0.007, 0.03, -0.024, -0.011, 0.031, -0.028, -0.017, 0.032, -0.031, -0.021, 0.035, -0.038, -0.016, 0.027, -0.043, -0.014, 0.027, -0.037, -0.014, 0.023, -0.032, -0.011, 0.024, -0.026, -0.01, 0.024, -0.021, -0.009, 0.023, -0.017, -0.008, 0.021, -0.018, -0.009, 0.018, -0.021, -0.011, 0.022, -0.011, -0.013, 0.025, -0.005, -0.02, 0.021, -0.018, -0.028, 0.005, -0.046, -0.026, 0.008, -0.035, -0.012, 0.022, 0.001, -0.01, 0.036, 0.016, -0.007, 0.057, 0.041, -0.008, 0.065, 0.045, -0.012, 0.054, 0.026, -0.006, 0.05, 0.036, -0.005, 0.057, 0.046, -0.005, 0.059, 0.047, -0.007, 0.059, 0.037, -0.008, 0.054, 0.028, -0.008, 0.047, 0.016, -0.009, 0.041, 0.015, -0.005, 0.037, 0.024, -0.004, 0.033, 0.019, -0.005, 0.025, 0.01, -0.002, 0.022, 0.017, -0.001, 0.012, 0.023, 0.002, 0.011, 0.036, 0.002, 0.007, 0.039, 0.002, 0.015, 0.029, 0.003, 0.023, 0.036, 0.004, 0.021, 0.051, 0.003, 0.017, 0.048, 0.004, 0.005, 0.056, 0.004, -0.001, 0.049, 0.002, 0.002, 0.039, -0.001, 0.001, 0.024, -0.002, -0.002, 0.013, -0.005, 0.002, 0.002, -0.005, -0.002, 0.005, -0.004, 0.007, 0.01, -0.006, 0.01, -0.0, -0.006, 0.003, -0.001, -0.004, 0.005, 0.007, -0.005, 0.015, 0.002, -0.007, 0.011, -0.007, -0.006, 0.007, -0.003, -0.005, 0.016, 0.001, -0.006, 0.02, -0.008, -0.007, 0.013, -0.013, -0.007, 0.014, -0.007, -0.006, 0.024, -0.009, -0.008, 0.023, -0.019, -0.009, 0.017, -0.019, -0.007, 0.022, -0.014, -0.007, 0.031, -0.015, -0.006, 0.033, -0.004, -0.005, 0.033, 0.005, -0.004, 0.023, 0.01, -0.002, 0.026, 0.022, -0.004, 0.038, 0.025, -0.002, 0.039, 0.038, -0.002, 0.05, 0.046, -0.001, 0.045, 0.051, 0.002, 0.034, 0.051, 0.005, 0.03, 0.061, 0.003, 0.033, 0.054, 0.001, 0.028, 0.043, 0.003, 0.018, 0.048, 0.005, 0.016, 0.052, 0.002, 0.018, 0.039, 0.001, 0.01, 0.034, 0.004, 0.001, 0.043, 0.003, 0.002, 0.041, 0.0, -0.001, 0.028, 0.002, -0.009, 0.029, 0.005, -0.015, 0.042, 0.008, -0.015, 0.055, 0.008, -0.022, 0.056, 0.011, -0.036, 0.06, 0.014, -0.042, 0.069, 0.012, -0.044, 0.057, 0.01, -0.03, 0.052, 0.012, -0.027, 0.066, 0.013, -0.034, 0.068, 0.01, -0.026, 0.056, 0.009, -0.014, 0.063, 0.01, -0.011, 0.071, 0.01, -0.001, 0.082, 0.014, -0.006, 0.096, 0.012, 0.007, 0.097, 0.01, 0.004, 0.08, 0.012, -0.01, 0.082, 0.014, -0.006, 0.096, 0.011, 0.004, 0.088, 0.01, -0.004, 0.076, 0.014, -0.013, 0.088, 0.013, -0.002, 0.095, 0.011, -0.002, 0.081, 0.012, -0.009, 0.08, 0.011, -0.019, 0.071, 0.008, -0.017, 0.057, 0.007, -0.025, 0.046, 0.003, -0.021, 0.033, -0.0, -0.015, 0.024, -0.003, -0.018, 0.013, -0.004, -0.011, 0.004, -0.006, -0.013, -0.005, -0.005, -0.007, 0.002, -0.004, 0.002, 0.004, -0.006, 0.004, -0.007, -0.007, 0.003, -0.012, -0.009, 0.0, -0.022, -0.009, 0.0, -0.023, -0.01, -0.008, -0.028, -0.009, -0.017, -0.022, -0.009, -0.017, -0.023, -0.009, -0.007, -0.02, -0.007, -0.011, -0.012, -0.006, -0.02, -0.01, -0.007, -0.013, -0.012, -0.006, -0.007, -0.005, -0.003, -0.016, 0.002, -0.004, -0.018, 0.0, -0.004, -0.008, 0.001, -0.002, -0.006, 0.01, 0.0, -0.015, 0.017, 0.002, -0.012, 0.029, 0.004, -0.02, 0.03, 0.002, -0.014, 0.021, 0.005, -0.025, 0.02, 0.003, -0.03, 0.016, -0.002, -0.017, 0.008, -0.005, -0.017, -0.001, -0.005, -0.012, -0.006, -0.006, -0.005, -0.013, -0.006, -0.011, -0.018, -0.007, -0.019, -0.014, -0.007, -0.013, -0.016, -0.008, -0.01, -0.025, -0.009, -0.02, -0.026, -0.009, -0.02, -0.024, -0.01, -0.023, -0.031, -0.01, -0.013, -0.031, -0.01, -0.01, -0.036, -0.01, -0.017, -0.043, -0.011, -0.01, -0.051, -0.01, -0.005, -0.044, -0.009, -0.011, -0.041, -0.009, -0.011, -0.044, -0.009, -0.013, -0.037, -0.008, -0.008, -0.039, -0.006, 0.001, -0.037, -0.006, 0.003, -0.029, -0.003, 0.01, -0.025, -0.005, 0.012, -0.019, -0.005, 0.021, -0.015, -0.005, 0.03, -0.016, -0.007, 0.037, -0.02, -0.013, 0.038, -0.019, -0.022, 0.048, -0.024, -0.025, 0.041, -0.031, -0.035, 0.052, -0.04, -0.039, 0.045, -0.042, -0.036, 0.034, -0.048, -0.048, 0.045, -0.06, -0.052, 0.049, -0.058, -0.063, 0.044, -0.067, -0.053, 0.036, -0.061, -0.05, 0.023, -0.064, -0.041, 0.007, -0.058, -0.033, 0.012, -0.048, -0.036, 0.028, -0.045, -0.032, 0.021, -0.043, -0.023, 0.009, -0.036, -0.02, 0.023, -0.025, -0.021, 0.031, -0.025, -0.015, 0.017, -0.022, -0.007, 0.014, -0.01, -0.007, 0.031, -0.001, -0.006, 0.028, -0.004, 0.002, 0.013, 0.003, 0.007, 0.02, 0.016, 0.004, 0.035, 0.021, 0.007, 0.026, 0.019, 0.012, 0.031, 0.031, 0.017, 0.021, 0.034, 0.018, 0.025, 0.032, 0.029, 0.027, 0.047, 0.038, 0.022, 0.06, 0.046, 0.042, 0.08]
        parser = ProdyNMDParser()
        
        eigenvalues, eigenvectors, header = parser.read(file_path)
        numpy.testing.assert_array_almost_equal(expected_eigenvalues, eigenvalues, 7)
        numpy.testing.assert_array_almost_equal(expected_eigenvector_0, eigenvectors[0], 7)
         
         
#    def test_convert_with_only_cas(self):
#        evec_file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","cr","imode_cart.evec")
#        prot_file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","cr","imode_model.pdb")
#        eigenvalues, eigenvectors = ImodServerFilesParser.read(evec_file_path)
#        atoms = prody.parsePDB(prot_file_path)
#        ProdyNMDWriter.write_CA("test", "test", atoms, eigenvalues, eigenvectors, center = False)
#        ProdyNMDWriter.write("test2", "test2", 
#                             prody.parsePDB(prot_file_path), 
#                             eigenvalues, eigenvectors, center = False)
#        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()