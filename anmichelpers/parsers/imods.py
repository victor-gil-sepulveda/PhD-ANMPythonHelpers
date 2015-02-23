"""
Created on 5/2/2015

@author: victor
"""
import numpy

class ImodServerFilesParser(object):
    """
    Reads eigenvalues/vectors generated from http://imods.chaconlab.org/
    """

    def __init__(self, params):
        pass
    
    @classmethod
    def read(cls, file_path):
        """
        
        TODO: Using streams and extracting functions would improve testability.
        """
        eigenvalues = []
        eigenvectors = []
        
        handler = open(file_path, "r")
        file_lines = handler.readlines()
        
        # We skip first line
        (eigenlength, num_eigenvectors) = cls.__get_number_and_length_of_eigenvectors(file_lines)
        current_line = 3
        while len(eigenvectors) != num_eigenvectors:
            # process eigenvalue
            eigenvalues.append(float(file_lines[current_line].split()[1]))
            current_line += 1
            
            # process eigenvector
            eigenvector = []
            while current_line < len(file_lines) and file_lines[current_line] != "****\n":
                eigenvector.extend(numpy.array(file_lines[current_line].split()).astype(float))
                current_line += 1

            if len(eigenvector) != eigenlength:
                print "[Warning ImodServerFilesParser::read] Read eigenvector does not have the expected length (read: %d, expected: %d)."%(len(eigenvector),eigenlength)
            
            eigenvectors.append(eigenvector)
            
            # skip ****
            current_line += 1
        return numpy.array(eigenvalues), numpy.array(eigenvectors)
    
    @classmethod        
    def __get_number_and_length_of_eigenvectors(self, file_lines):
        parts = file_lines[1].split()
        return (int(parts[0]), int(parts[3]))
    