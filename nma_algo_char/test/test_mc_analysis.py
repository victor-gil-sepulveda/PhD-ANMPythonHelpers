'''
Created on 16/12/2015

@author: victor
'''
import unittest
from nma_algo_char.mode_application_analysis import process_energy_differences
import numpy


class Test(unittest.TestCase):


    def test_MC_analysis(self):
        
        initial_u = [-14080.7, -14290.6, -14295.4, -14300.8, -14312.3, -14317.9, -14319.1, -14329.3, -14332.2, -14339.8, -14342.8, -14341.1]
        final_u   = [-14290.6, -14295.4, -14300.8, -14312.3, -14317.9, -14319.1, -14329.3, -14332.2, -14339.8, -14342.8, -14341.1, -14340.1]
        
        data = {"e_after":numpy.array(final_u), "e_before":numpy.array(initial_u)}
        
        print process_energy_differences(data)
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_MC_analysis']
    unittest.main()