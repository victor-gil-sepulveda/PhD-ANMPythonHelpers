import os
import sys
import numpy
import prody
from optparse import OptionParser
from anmichelpers.parsers.imods import ImodServerFilesParser
from anmichelpers.parsers.bahar import BaharServerFilesParser
from anmichelpers.writers.pronmd import ProdyNMDWriter
from anmichelpers.parsers.pronmd import ProdyNMDParser
from anmichelpers.tools.atoms import atoms_from_header, get_CA_atoms

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("--evals", dest="input_evals")
    parser.add_option("--beta", dest="input_beta")
    parser.add_option("--prot", dest="protein_path")
    parser.add_option("--ca_only", action="store_true", default=False, dest="ca_only")
    
    parser.add_option("-o", dest="output")
    (options, args) = parser.parse_args()
    
    name, extension = os.path.splitext(options.input)

    betas = None
    
    header = {}
    if extension == ".nmd":
        eigenvalues, eigenvectors, header = ProdyNMDParser.read(options.input)
    
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
    
    new_header = {}
    if extension != ".nmd":
        # Get a new header from the structural info
        structure = prody.proteins.pdbfile.parsePDB(options.protein_path)
        header["name"] = options.input
        header["type"] = "cc:pca"
        header["atomnames"] = structure.getNames()
        header["coordinates"] = structure.getCoordsets()[0].flatten()
        header["chainids"] = structure.getChids()
        if betas is not None:
            header["bfactors"] = betas
        else:
            header["bfactors"] = structure.getBetas()
        header["resnames"] = structure.getResnames()
        header["resids"] = structure.getResindices()
        
    final_evecs = eigenvectors
    if options.ca_only:
        # Get only the CAs
        atoms = atoms_from_header(header, eigenvectors)
        ca_atoms = get_CA_atoms(atoms)
        # remount eigenvalues and header
        new_header = {
                      "atomnames":[],
                      "resids":[],
                      "coordinates":[]
                      }
        new_modes = []
        for ca in ca_atoms:
            new_header["atomnames"].append("CA")
            new_header["resids"].append(ca.resid)
            new_header["coordinates"].extend(ca.coords)
            new_modes.append(ca.mode_v)
        
        final_evecs = numpy.array(new_modes).T
        
    ProdyNMDWriter.write(options.output, eigenvalues, final_evecs, header)