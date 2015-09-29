import sys
import os
from optparse import OptionParser
from anmichelpers.parsers.imods import ImodServerFilesParser
from anmichelpers.parsers.bahar import BaharServerFilesParser
from anmichelpers.writers.pronmd import ProdyNMDWriter
import prody

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("--evals", dest="input_evals")
    parser.add_option("--beta", dest="input_beta")
    parser.add_option("--prot", dest="protein_path")
    
    
    parser.add_option("-o", dest="output")
    (options, args) = parser.parse_args()
    
    name, extension = os.path.splitext(options.input)

    betas = None
    
    if extension == ".nmd":
        print "There's no need to convert this file."
        sys.exit()
    
    elif extension == ".evec":
        # Imods file
        eigenvalues, eigenvectors = ImodServerFilesParser.read(options.input)
        
    elif extension == ".slwevs":
        # Bahar server file
        eigenvalues, eigenvectors =  BaharServerFilesParser.read(options.input_evals, 
                                                                options.input)
        if options.input_beta is not None:
            betas  = BaharServerFilesParser.read_beta(options.input_beta)
    else:
        print "I don't know how to convert this file."
        sys.exit()
    
    structure = prody.proteins.pdbfile.parsePDB(options.protein_path)
    ca_indices = ProdyNMDWriter.get_alpha_indices(structure)
    ca_evecs = ProdyNMDWriter.filter_eigvecs(ca_indices, eigenvectors)
    ProdyNMDWriter.write(options.output, "conversion from %s"%options.input, 
                         structure.select("name CA").copy(), 
                         ca_evecs, eigenvectors, betas)