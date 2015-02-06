"""
Created on Feb 6, 2015

@author: victor
"""
import prody

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
    def write(cls, file_path, name, atoms, eigenvalues, eigenvectors, betas = None):
        nma = prody.dynamics.anm.ANM()
        nma.setEigens(eigenvectors.T, eigenvalues)
        
        atoms.setTitle(name)
        
        if betas is not None:
            atoms.setBetas(betas)
        else:
            # new beta calculation is weighted by current (experimental)\
            # beta, so we initialize it to 1.
            atoms.setBetas([1.]*atoms.numAtoms()) 
            atoms.setBetas(prody.dynamics.analysis.calcTempFactors(nma, atoms))
        
        prody.dynamics.nmdfile.writeNMD(file_path, nma, atoms)
        
        