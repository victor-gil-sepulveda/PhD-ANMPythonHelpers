"""
Created on Feb 6, 2015

@author: victor
"""

import prody

class ProdyNMDParser(object):
    """
    Uses prody to parse eigenvectors and eigenvalues from a NMD file.
    """
    def __init__(self):
        self.anm = None
        self.atoms = None
    
    def read(self, nmd_file):
        self.anm, self.atoms = prody.parseNMD(nmd_file)
    
        return self.anm.getEigvals(), self.anm.getEigvecs().T