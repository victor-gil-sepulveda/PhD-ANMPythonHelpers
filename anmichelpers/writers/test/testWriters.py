"""
Created on Feb 6, 2015

@author: victor
"""
import unittest
import anmichelpers.writers.test
import os
from anmichelpers.parsers.pronmd import ProdyNMDParser
from anmichelpers.writers.pronmd import ProdyNMDWriter
import prody
import numpy

class Test(unittest.TestCase):


    def testProdyNMDWriter(self):
        evec_path = os.path.join(anmichelpers.writers.test.__path__[0], "data","eigenvectors.npy")
        eval_path = os.path.join(anmichelpers.writers.test.__path__[0], "data","eigenvalues.npy")
        
        output_path = os.path.join(anmichelpers.writers.test.__path__[0], "data","calculated.nmd")
        expected_path = os.path.join(anmichelpers.writers.test.__path__[0], "data","expected.nmd")
        
        pdb_file_path = os.path.join(anmichelpers.writers.test.__path__[0], "data","protein.pdb")
        
        eigenvalues = numpy.load(eval_path)
        eigenvectors = numpy.load(evec_path)
        
        ProdyNMDWriter.write(output_path, 
                             "TEST",
                             prody.parsePDB(pdb_file_path), 
                             eigenvalues, 
                             eigenvectors)
        
        self.assertEqual("".join(open(expected_path).readlines()[1:]) , "".join(open(output_path).readlines())[1:])
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testProdyNMDWriter']
    unittest.main()