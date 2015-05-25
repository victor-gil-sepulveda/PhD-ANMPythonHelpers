'''
Created on 22/05/2015

@author: user
'''
import unittest
from calculate_PCA_from_coordinates import get_seq_positions_with_known_residues,\
    extract_coordinates_of_known_residues, subtract_mean, calculate_weights,\
    calc_mean_of_known_atoms, change_unknown_by_mean
import numpy


class TestPCAScriptStuff(unittest.TestCase):


    def test_get_seq_positions_with_known_residues(self):
        
        sequences = [
                     ["AAA","AAA","AAA","DDD","AAA","GAP","AAA"],
                     ["BBB","AAA","GAP","AAA","EEE","AAA","AAA"],
                     ["AAA","CCC","AAA","AAA","GAP","AAA","FFF"]
                     ]
        expected_known_pos = [0, 1, 3, 6]
        
        known_pos = get_seq_positions_with_known_residues(sequences)
        numpy.testing.assert_equal(expected_known_pos, known_pos)
    
    def test_extract_coordinates_of_known_residues(self):
        known_pos = [0,  2]
        coords = [ 
                  [1,2,3,     4,5,6,    7,8,9,   10,11,12],
                  [13,14,15, 16,17,18, 19,20,21, 22,23,24],
                  [25,26,27, 28,29,30, 31,32,33, 34,35,36]
                 ]
        
        expected_coords = [  [ 1,  2,  3,  7,  8,  9],
                             [13, 14, 15, 19, 20, 21],
                             [25, 26, 27, 31, 32, 33]]
                                    
        extracted_coords = extract_coordinates_of_known_residues(known_pos,coords)
        numpy.testing.assert_equal(expected_coords, extracted_coords)
        
    def test_subtract_mean(self):
        coords = [  [ 1,  2,  3,  7,  8,  9],
                    [13, 14, 15, 19, 20, 21],
                    [25, 26, 27, 31, 32, 33]]
        coords = numpy.array(coords)
        exp_mean = [ 13., 14., 15., 19., 20., 21.]
        exp_coords = [[-12]*6,[0]*6,[12]*6]
        subtract_mean(coords, exp_mean)
        numpy.testing.assert_equal( exp_coords, coords)
    
    def test_calculate_weights(self):
        sequences = [
                     ["AAA","AAA","AAA","DDD","AAA","GAP","AAA"],
                     ["BBB","AAA","GAP","AAA","EEE","AAA","AAA"],
                     ["AAA","CCC","AAA","AAA","GAP","AAA","FFF"]
                     ]
        expected_weights = [
                            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 
                            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 
                            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                            ]

        numpy.testing.assert_equal( expected_weights, calculate_weights(sequences).tolist())
    
    def test_calc_mean_of_known_atoms(self):
        sequences = [
                     ["AAA","AAA"],
                     ["BBB","GAP"],
                     ["GAP","CCC"]
                     ]
        coords = [  [ 1,  2,  3,  7,  8,  9],
                    [13, 14, 15, 19, 20, 21],
                    [25, 26, 27, 31, 32, 33]]
        expected_mean = [7.0, 8.0, 9.0, 19.0, 20.0, 21.0]
        calculated_mean =  calc_mean_of_known_atoms(sequences, coords)
        numpy.testing.assert_equal( expected_mean, calculated_mean)
        
    def test_change_unknown_by_mean(self):
        sequences = [
                     ["AAA","AAA"],
                     ["BBB","GAP"],
                     ["GAP","CCC"]
                     ]
        coords = [  [ 1,  2,  3,  7,  8,  9],
                    [13, 14, 15, 0.,0.,0.],
                    [0.,0.,0., 31, 32, 33]]
        
        mean_coords = [7.0, 8.0, 9.0, 19.0, 20.0, 21.0]
        
        expected_coords = [[1, 2, 3, 7, 8, 9], 
                           [13, 14, 15, 19.0, 20.0, 21.0], 
                           [7.0, 8.0, 9.0, 31, 32, 33]]
        
        change_unknown_by_mean(sequences, coords, mean_coords)
        
        numpy.testing.assert_equal( coords, expected_coords)
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()