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

def extract_coordinates_from_known_residues(known, coords):
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
        
        if atom_values == []:
            pos_mean.extend([0.,0.,0.])
        else:
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

def calc_sklearn_pca(mean_shifted_coords):
    """
    Uses scikit to do the calculations. Returns the percentage of explained variance 
    and the components.
    """
    from sklearn.decomposition import PCA
    pca = PCA(n_components=10, copy=True, whiten=False)
    pca.fit(mean_shifted_coords)
    return pca.explained_variance_ratio_,pca.components_ 

def calc_regular_pca(mean_shifted_coords, order = False, cast = False):
    """
    PCA calculated with regular numpy operations.
    """
    cov = numpy.dot(mean_shifted_coords.T, mean_shifted_coords.conj()) / (len(mean_shifted_coords)-1)
    #cov = numpy.cov(mean_shifted_coords.T)
    
    w, v = numpy.linalg.eig(cov)
    
    if order:
        eig_pairs = [(w[i], v[:,i]) for i in range(len(w))]
        eig_pairs.sort()
        eig_pairs.reverse()
        w, v =  numpy.array([w for w,_ in eig_pairs[0:10]]), numpy.array([v for _,v in eig_pairs[0:10]])
        v = v.T
        
    if cast:
        return w[0:10].astype(float), v[:, 0:10].T.astype(float)
    else:
        return w[0:10], v[:, 0:10].T


NUMBER_OF_COMPONENTS = 10

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-o", dest="output")
    parser.add_option("--ref", dest="ref", type = "int")
    parser.add_option("--type", dest="type")
    (options, args) = parser.parse_args()
    
    assert options.type in ["REGULAR", "REGULAR_ORDERED", "SKLEARN", "REGULAR_COV", "EMPCA"]
    
    coords = numpy.loadtxt(options.input+".coords")
    sequences = load_sequences(options.input+".seq")
    original_coords = numpy.copy(coords[options.ref])
    original_sequence = sequences[options.ref]
    
    # Extract the coordinates we know
    known_residues = get_seq_positions_with_known_residues(sequences)
    
    # Extraer bien las coordenadas
    known_coords = extract_coordinates_from_known_residues(known_residues, coords)
    
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
    numpy.savetxt("coords_mean", coords , fmt = "%.4f")
    
    # Recalc mean for all values
    recalcd_mean = numpy.mean(coords, axis = 0)
    numpy.savetxt("mean", recalcd_mean , fmt = "%.4f")
    
    # Subtract mean
    coords = coords - recalcd_mean
    
    if options.type == "EMPCA":
        # Generate weights for emPCA
        all_weights =  calculate_weights(sequences)
        
        # Calculate emPCA
        pca_model = empca(coords,
                          weights = numpy.array(all_weights), 
                          nvec = 10)
        
        eigvals = numpy.array([float(NUMBER_OF_COMPONENTS)]*NUMBER_OF_COMPONENTS)- numpy.array(range(NUMBER_OF_COMPONENTS))
        eigvecs = pca_model.eigvec
    elif options.type == "REGULAR":
        eigvals, eigvecs = calc_regular_pca(coords, order = False, cast = True)
    elif options.type == "REGULAR_ORDERED":
        eigvals, eigvecs = calc_regular_pca(coords, order = True, cast = True)
    elif options.type == "SKLEARN":
        eigvals, eigvecs = calc_sklearn_pca(coords)
        
    # Write the nmd file
    header = {
              "type":"cc:pca",
              "title": options.output,
              "coordinates":original_coords,
              "resnames":original_sequence
    }
    
    ProdyNMDWriter.write(options.output, eigvals, eigvecs, header)
    