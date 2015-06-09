"""
Created on Feb 6, 2015

@author: victor
"""
import numpy

class ProdyNMDWriter(object):
    """
    (from http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/ and 
    http://prody.csb.pitt.edu/manual/reference/dynamics/nmdfile.html#nmd-format)
    
    ----- NMD format--------
    NMD files (extension .nmd) are plain text files that contain at least normal mode and coordinate data. 
    PCA, EDA, NMA, ANM, or GNM data can be stored in NMD files. Recognized data fields are listed below.
    
    Mandatory:
    *coordinates: Coordinates must be provided in one line as a list of decimal numbers. Number of atoms in 
    the system is deduced from size of this data line. 
    *mode: Normal mode array. Each normal mode array must be provided in one line as a list of decimal numbers. 
    Mode array may be preceded by mode index and mode length (square root of variance or inverse frequency). 
    
    Extra:
    *title: A title for the dataset. 
    *name: Atom names. Default is "CA" for all atoms. 
    *resnames: Residue names. Default value is "GLY". 
    *chainids: Chain identifiers. Default value is "A". 
    *resids: Residue numbers. If this data line if not found, residue numbers are started from 1 and incremented 
    by one for each atom. 
    *betas: Beta factors. Default value is 0 (zero). B-factors are used to color the protein representation.
    
    Note that all data types must be listed in a single line. The size of data lines must match the number of atoms 
    in the system (the size of coordinates line).
    """
    
    
    def __init__(self,  params):
        pass

    @classmethod
    def write(cls, name, eigenvalues, eigenvectors, header = None):
        handler = open(name+".nmd","w")
        for tag in header:
            if tag in ["name", "type", "title"]:
                handler.write("%s %s\n"%(tag, header[tag]))
            else:
                handler.write("%s %s\n"%(tag, " ".join([str(h) for h in header[tag]])))
        
        for i in range(len(eigenvalues)):
            handler.write("mode %s %s\n"%(eigenvalues[i], " ".join([str(e) for e in eigenvectors[i]])))

    @classmethod
    def write_CA(cls, file_path, name, atoms, eigenvalues, eigenvectors, betas = None, center = False):
        indices = cls.get_alpha_indices(atoms)
        filtered_eigenvectors = cls.filter_eigvecs(indices, eigenvectors)
        cls.write(file_path, name, atoms.select("name CA").copy(), eigenvalues, filtered_eigenvectors, betas, center)
    
    @classmethod
    def get_alpha_indices(cls, structure):
        return numpy.where(structure.getNames() == "CA")[0]
    
    @classmethod
    def filter_eigvecs(cls, indices, evecs):
        print evecs
        new_evecs = []
        for evec in evecs:
            new_evec = []
            for index in indices:
                new_evec.append(evec[index*3])
                new_evec.append(evec[index*3+1])
                new_evec.append(evec[index*3+2])
            new_evecs.append(new_evec)
        return numpy.array(new_evecs)