'''
Created on 22/05/2015

@author: user
'''
import unittest
from calculate_PCA_from_coordinates import get_seq_positions_with_known_residues,\
    extract_coordinates_from_known_residues, subtract_mean, calculate_weights,\
    calc_mean_of_known_atoms, change_unknown_by_mean
import numpy
from local_to_global_cc_pca import fill_with_zeros, extract_evecs
from obtain_structurally_aligned_angles import get_angles_and_sequence_for_pdb_with_mapping
from numpy.core.numeric import inf
from local_to_global_ic_pca import fill_with_zeros as fill_with_zeros_ic
from local_to_global_ic_pca import extract_evecs as extract_evecs_ic
import os
import test.data
import prody

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
                                    
        extracted_coords = extract_coordinates_from_known_residues(known_pos,coords)
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
        
    def test_eigvec_padding(self):
        evecs = [[1,2,3,4,5,6],
                 [7,8,9,10,11,12]]
        
        expected_evecs =  [
         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
         ]
        
        #[X,X,X, A, B, X, X]
        padded_evecs =  fill_with_zeros(evecs, 4, 5, 7)
        numpy.testing.assert_equal(expected_evecs, padded_evecs)
        
        expected_evecs = [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0], 
                          [7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 0.0, 0.0, 0.0]]
        numpy.testing.assert_equal(expected_evecs,fill_with_zeros(evecs, 1, 2, 3))
        
        expected_evecs = [[0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 
                          [0.0, 0.0, 0.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]]
        numpy.testing.assert_equal(expected_evecs,fill_with_zeros(evecs, 2, 3, 3))
        
    def test_eigvec_ic_padding(self):
        evecs = [[1,2,3,4,],
                 [7,8,9,10]]
        
        expected_evecs =  [
         [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0], 
         [0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 8.0, 9.0, 10.0, 0.0, 0.0, 0.0, 0.0]
         ]
        
        #[X,X,X, A, B, X, X]
        padded_evecs =  fill_with_zeros_ic(evecs, 4, 5, 7).tolist()
        numpy.testing.assert_equal(expected_evecs, padded_evecs)
        
        initial_evecs = [[1,2,3],
                         [7,8,9]]
        expected_evecs = [[1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                          [7.0, 8.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        numpy.testing.assert_equal(expected_evecs, fill_with_zeros_ic(initial_evecs, 1, 2, 5))
        
        expected_evecs = [[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0], 
                          [0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 8.0, 9.0, 10.0]]
        numpy.testing.assert_equal(expected_evecs,  fill_with_zeros_ic(evecs, 4, 5, 5))
        
    
    def test_extract_evecs(self):
        evecs =  [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ]
        expected_evecs = [[1,2,3,4,5,6],
                          [7,8,9,10,11,12]]
        
        extracted_evecs =  extract_evecs(evecs, 4, 5)
        numpy.testing.assert_equal(expected_evecs, extracted_evecs)
    
    def test_extract_evecs_ic(self):
        evecs =  [
        # psi  phi  psi  phi  psi  phi  psi  phi  psi   phi  psi  phi  psi
         [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0,  0.0, 0.0, 0.0, 0.0], 
         [0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 8.0, 9.0, 10.0, 0.0, 0.0, 0.0, 0.0]
         ]
        sequence = ["XXX","XXX"]
        expected_evecs = [[2, 3, 4],
                          [8, 9, 10]]
        
        extracted_evecs =  extract_evecs_ic(evecs, 4, 5, sequence)
        numpy.testing.assert_equal(expected_evecs, extracted_evecs)

        evecs =  [
        # psi  phi  psi   phi   psi   phi   psi
         [1.0, 2.0,  3.0, 4.0,  5.0,  6.0,  7.0], 
         [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
         ]
        
        expected_evecs = [[5.0,   6.0,  7.0], 
                         [ 12.0, 13.0, 14.0]]
        numpy.testing.assert_equal(expected_evecs,  extract_evecs_ic(evecs, 3, 4, sequence))
        
        expected_evecs =  [[1.0, 2.0,  3.0], 
                           [8.0, 9.0, 10.0]]
        numpy.testing.assert_equal(expected_evecs, extract_evecs_ic(evecs, 1, 2, sequence))
        
        
        sequence = ["XXX","PRO"]
        expected_evecs = [[2,  4],
                          [8, 10]]
        evecs =  [
        # psi  phi  psi  phi  psi  phi  psi  phi  psi   phi  psi  phi  psi
         [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0,  0.0, 0.0, 0.0, 0.0], 
         [0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 8.0, 9.0, 10.0, 0.0, 0.0, 0.0, 0.0]
         ]
        extracted_evecs =  extract_evecs_ic(evecs, 4, 5, sequence)
        numpy.testing.assert_equal(expected_evecs, extracted_evecs)

    
    def test_angle_extraction(self):

        expected_angles = numpy.loadtxt(os.path.join(test.data.__path__[0], "helix.angles")) 
        expected_sequence = ['PHE', 'ASN', 'GLU', 'GLU', 'LYS', 'MET', 'HID', 'GLN', 'SER', 'ALA', 'MET', 'TYR', 'GLU']
        expected_nums = range(1,14)
        
        class options:
            use_degrees = False
            initial_res = 1
            final_res = 13
        
        ### One to one mapping
        pdb_path = os.path.join(test.data.__path__[0], "helix.pdb")
        pdb = prody.parsePDB(pdb_path)
        res_mapping = dict(zip(range(1,14),range(1,14)))
        angles, sequence, res_selection = get_angles_and_sequence_for_pdb_with_mapping(pdb, res_mapping, "**", options)
        
        numpy.testing.assert_almost_equal(angles, expected_angles, 10)
        self.assertSequenceEqual(expected_sequence, sequence)
        numpy.testing.assert_equal(res_selection, expected_nums)
        
        ### One residue less in the target sequence (due to a gap in the 
        ### structure)
        
        ['PHE', 'ASN', 'GLU', 'GLU', 'LYS', 'MET', 'HID', 'GLN', 'SER', 'ALA', 'MET', 'TYR', 'GLU']
        ['PHE', 'GLU', 'GLU', 'LYS', 'MET', 'HID', 'GLN', 'SER', 'ALA', 'MET', 'TYR', 'GLU', 'GAP']
        
        expected_angles = numpy.loadtxt(os.path.join(test.data.__path__[0], "helix_gap_1.angles")) 
        expected_sequence = ['PHE', 'GLU', 'GLU', 'LYS', 'MET', 'HID', 'GLN', 'SER', 'ALA', 'MET', 'TYR', 'GLU', 'GAP']
        expected_nums = range(1,14)
        res_mapping = dict(zip(range(1,14),[1,3,4,5,6,7,8,9,10,11,12,13]))
        
        pdb_path = os.path.join(test.data.__path__[0], "helix_gap_1.pdb")
        pdb = prody.parsePDB(pdb_path)
        angles, sequence, res_selection = get_angles_and_sequence_for_pdb_with_mapping(pdb, res_mapping, "**", options)
        
        numpy.testing.assert_almost_equal(angles, expected_angles, 10)
        self.assertSequenceEqual(expected_sequence, sequence)
        
        ### One residue less due to an alignment gap (due to a gap 
        ### in the alignment)
        expected_angles = numpy.loadtxt(os.path.join(test.data.__path__[0], "helix_gap_2.angles")) 
        expected_sequence = ['PHE', 'ASN', 'GLU', 'LYS', 'MET', 'HID', 'GLN', 'SER', 'ALA', 'MET', 'TYR', 'GLU', 'GAP']
        expected_nums = range(1,14)
        res_mapping = dict(zip(range(1,14),[1,2,4,5,6,7,8,9,10,11,12,13]))
        
        pdb_path = os.path.join(test.data.__path__[0], "helix.pdb")
        pdb = prody.parsePDB(pdb_path)
        angles, sequence, res_selection = get_angles_and_sequence_for_pdb_with_mapping(pdb, res_mapping, "**", options)
        numpy.testing.assert_almost_equal(angles, expected_angles, 10)
        self.assertSequenceEqual(expected_sequence, sequence)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()