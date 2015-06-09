"""
Created on 15/05/2015

@author: vgil
"""

from optparse import OptionParser
import numpy
from pca.geopca import GeoPCA
from anmichelpers.writers.pronmd import ProdyNMDWriter

def load_sequences(file_path):
    """
    Get sequences from a file:
    [ struct 1 seq (master),
    struct 2 seq,
    struct 3 seq,
    ...
    struct n seq] 
    """
    sequences = []
    for line in open(file_path):
        sequences.append(line.split())
    return sequences

def treat_prolines_in_eigenvectors(master_seq, old_eigenvectors):
    """
    Whenever a proline is encountered in the sequence, only the psi value of the 
    eigenvector is used, as in our model prolines do not have phi movement.
    
    We have to take into account that angles is (r*2)-1 (first residue only 
    has psi angle).
    
    """
    old_eigenvectors = old_eigenvectors.T
    new_eigenvectors = []
    for i, res in enumerate(master_seq):
        if i == 0:
            #only psi if in first position
            new_eigenvectors.append(old_eigenvectors[0])
        else:
            offset = i*2
            if res == "PRO":
                #only psi
                new_eigenvectors.append(old_eigenvectors[offset])
            else:
                # phi
                new_eigenvectors.append(old_eigenvectors[offset - 1])
                # psi
                new_eigenvectors.append(old_eigenvectors[offset])
    return numpy.array(new_eigenvectors).T

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-s", dest="sequence")
    parser.add_option("-o", dest="output")
    (options, args) = parser.parse_args()
    
    eigenvalues, eigenvectors = GeoPCA.load(options.input)
    
    sequences = load_sequences(options.sequence)
    
    pro_eigenvectors = treat_prolines_in_eigenvectors(sequences[0], eigenvectors)
    
    # And save it as nmd
    header = {"type": "ic:pca"}
    ProdyNMDWriter.write(options.output, eigenvalues, pro_eigenvectors, header)
    