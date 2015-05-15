'''
Created on 15/05/2015

@author: vgil
'''
import numpy
import math
class GeoPCA(object):

    def __init__(self):
        pass
    
    @classmethod
    def load(cls, path):
        file_lines = open(path).readlines()

        eigenvalues = []
        eval_parts = file_lines[1].split(",")
        for part in eval_parts:
            eigenvalues.append(float(part.split(":")[1]))

        eigenvectors = []
        for line in file_lines[2:]:
            parts = line.split(",")
            for part in parts:
                eigenvectors.append((float(part) * 180.)/math.pi )
        return numpy.array(eigenvalues), numpy.array(eigenvectors).T
            
print GeoPCA.load("geoPCA.out")