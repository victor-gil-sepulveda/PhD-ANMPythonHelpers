'''
Created on 12/1/2016

@author: victor
'''
from anmichelpers.parsers.pronmd import ProdyNMDParser
from prody.proteins.pdbfile import parsePDB
from anmichelpers.comparison.comparison import cumulative_overlap
import sys

def calculate_transition_vecs(first_conf, second_conf):
    # From first conf to second conf
    # Calculate transition vectors
    print first_conf, second_conf
    first_struct = parsePDB(first_conf).ca
    second_struct = parsePDB(second_conf).ca
    return (second_struct.getCoords() - first_struct.getCoords()).flatten()

if __name__ == '__main__':
    num_modes = int(sys.argv[1])
    types  = ["cc_closed", "ic_closed", "cc_open", "ic_open"]
    files = [(t, "%s_norm.nmd"%t) for t in types]
    all_eigenvectors = {}
    open_to_vopen = calculate_transition_vecs("open.pdb", "very_open.pdb")
    open_to_vclosed = calculate_transition_vecs("open.pdb", "very_closed.pdb")
    closed_to_vopen = calculate_transition_vecs("closed.pdb", "very_open.pdb")
    closed_to_vclosed = calculate_transition_vecs("closed.pdb", "very_closed.pdb")
    print files
    for t,f in files:
        _, _eigenvectors, _ = ProdyNMDParser.read(f)
        
        c_ov = cumulative_overlap(open_to_vopen, _eigenvectors[0:num_modes])
        print t, "w open to v. open %.3f"%c_ov
        c_ov = cumulative_overlap(open_to_vclosed, _eigenvectors[0:num_modes])
        print t, "w open to v. closed %.3f"%c_ov
        
        c_ov = cumulative_overlap(closed_to_vopen, _eigenvectors[0:num_modes])
        print t, "w closed to v. open %.3f"%c_ov
        c_ov = cumulative_overlap(closed_to_vclosed, _eigenvectors[0:num_modes])
        print t, "w closed to v. closed %.3f"%c_ov
