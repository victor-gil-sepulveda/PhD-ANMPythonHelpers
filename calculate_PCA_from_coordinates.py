'''
Created on 22/05/2015

@author: user
'''
from optparse import OptionParser
import numpy
from convert_geo_pca_to_nmd import load_sequences
from pyRMSD.RMSDCalculator import RMSDCalculator
from pca.emPCA import empca
from anmichelpers.writers.pronmd import ProdyNMDWriter

def get_seq_positions_with_known_residues(sequences):
    known_residues = []
    for i,residues in enumerate(numpy.array(sequences).T):
        if not "GAP" in residues:
            known_residues.append(i)
    return known_residues

def extract_coordinates_of_known_residues(known, coords):
    all_new_coords = []
    for coord_set in coords:
        new_coords = []
        for res_id in known:
            offset = res_id*3
            new_coords.append(coord_set[offset])
            new_coords.append(coord_set[offset+1])
            new_coords.append(coord_set[offset+2])
        all_new_coords.append(new_coords)
    return numpy.array(all_new_coords)

def subtract_mean(coords, coords_mean):
    for i in range(len(coords)):
        coords[i] = coords[i] - coords_mean

def calculate_weights(sequences):
    all_weights = []
    for seq in sequences:
        weights = []
        for res in seq:
            if res == "GAP":
                weights.extend([0.0, 0.0, 0.0])
            else:
                weights.extend([1.0, 1.0, 1.0])
        all_weights.append(weights)
    all_weights = numpy.array(all_weights)
    return all_weights

def calc_mean_of_known_atoms(sequences, coords):
    sequences = numpy.array(sequences)
    pos_mean = []
    for i, res_seqs in enumerate(sequences.T):
        atom_values = []
        for j, res_type in enumerate(res_seqs):
            if res_type != "GAP":
                atom_values.append(coords[j][i*3:(i+1)*3])
        pos_mean.extend(numpy.mean(atom_values,0))
    return pos_mean

def change_unknown_by_mean(sequences, coords, mean_coords):
    for i, res_seqs in enumerate(sequences):
        for j, res_type in enumerate(res_seqs):
            if res_type == "GAP":
                offset = j*3
                coords[i][offset] =  mean_coords[offset]
                coords[i][offset+1] = mean_coords[offset+1]
                coords[i][offset+2] = mean_coords[offset+2]

def calc_regular_pca(mean_shifted_coords):
    cov = numpy.dot(mean_shifted_coords.T, mean_shifted_coords.conj()) / (len(mean_shifted_coords)-1)
    # TODO! Order is not guaranteed.  Order the eigenvectors by eigenvalue!
    w, v = numpy.linalg.eig(cov)
    return w[0:10].astype(float), v[:, 0:10].T.astype(float)

NUMBER_OF_COMPONENTS = 10

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-s", dest="sequence")
    parser.add_option("-o", dest="output")
    parser.add_option("--empca", action = "store_true", default = False, dest="empca")
    (options, args) = parser.parse_args()
    
    coords = numpy.loadtxt(options.input)
    sequences = load_sequences(options.sequence)
    original_coords = numpy.copy(coords[0])
    
    # Extract the coordinates we know
    known_residues = get_seq_positions_with_known_residues(sequences)
    
    # Extraer bien las coordenadas
    known_coords = extract_coordinates_of_known_residues(known_residues, coords)
    
    # Do an iterative superposition of that coordinates, but move all coordinates
    known_coords = numpy.reshape(known_coords, (known_coords.shape[0], 
                                                known_coords.shape[1]/3,
                                                3))
    coords = numpy.reshape(coords, (coords.shape[0], 
                                    coords.shape[1]/3,
                                    3))
    calculator = RMSDCalculator("QTRFIT_SERIAL_CALCULATOR", known_coords, coords)
    calculator.iterativeSuperposition() 

    # Reshape iterposed coordinates    
    coords = numpy.reshape(coords, (coords.shape[0], coords.shape[1]*3))
    
    # Calculate known coordinates mean
    known_mean_coords = calc_mean_of_known_atoms(sequences, coords)
    
    # Change unknown coordinates by mean
    change_unknown_by_mean(sequences, coords, known_mean_coords)
    
    # Subtract mean
    subtract_mean(coords, known_mean_coords)
    
    if options.empca:
        # Generate weights for emPCA
        all_weights =  calculate_weights(sequences)
        
        # Calculate emPCA
        pca_model = empca(coords,
                          weights = numpy.array(all_weights), 
                          nvec = 10)
        
        eigvals = numpy.array([float(NUMBER_OF_COMPONENTS)]*NUMBER_OF_COMPONENTS)- numpy.array(range(NUMBER_OF_COMPONENTS))
        eigvecs = pca_model.eigvec
    else:
        eigvals, eigvecs = calc_regular_pca(coords)
    
    # Write the nmd file
    header = {
              "type":"cc:pca",
              "coordinates":original_coords,
              "resnames":sequences[0]
    }
    
    ProdyNMDWriter.write(options.output, eigvals, eigvecs, header)
    