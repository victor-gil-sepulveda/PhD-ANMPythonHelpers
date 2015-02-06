"""
Created on 5/2/2015

@author: victor
"""
import numpy

class BaharServerFilesParser(object):
    """
    Reads eigenvalues/vectors generated from http://anm.csb.pitt.edu/cgi-bin/anm2/anm2.cgi
    """

    def __init__(self):
        pass
    
    @classmethod
    def read(cls, eigenvalues_file_path, eigenvectors_file_path):
        eigenvalues = []
        if eigenvalues_file_path is not None and eigenvalues_file_path != "":
            eigenvalues = cls.__read_eigenvalues(eigenvalues_file_path)
        
        return (eigenvalues, cls.__read_eigenvectors(eigenvectors_file_path))
    
    @classmethod
    def read_beta(cls, beta_factors_file):
        """
        Beta factor files ("*.bfactors") have the same format than eigenvalues file.
        """
        return cls.__read_eigenvalues(beta_factors_file)
    
    @classmethod
    def __read_eigenvalues(cls,  eigenvalues_file_path):
        """
        Reads a file with an enumeration of eigenvalues (6+30). Ex.
            1 Eigval_1
            2 Eigval_2
            ...
            N Eigval_N
        
        First 6 eigenvalues must be 0.
        """
        return numpy.loadtxt(eigenvalues_file_path).T[1]

    @classmethod  
    def __read_eigenvectors(cls, eigenvectors_file_path):
        """
        Reads a file containing eigenvectors (20) written in columns.
        """ 
        return numpy.loadtxt(eigenvectors_file_path).T[1:]
    
    

    