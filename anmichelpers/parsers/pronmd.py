"""
Created on Feb 6, 2015

@author: victor
"""
import anmichelpers.parsers.lineParsers as lparsers
import numpy

line_parsers = {
                "type": lparsers.single_string_parser, # in house tag
                "mode": lparsers.mode_parser,
                "coordinates": lparsers.float_array_parser,
                "atomnames": lparsers.string_array_parser,
                "name": lparsers.single_string_parser,
                "chainids": lparsers.string_array_parser,
                "bfactors": lparsers.float_array_parser,
                "resnames": lparsers.string_array_parser,
                "resids": lparsers.int_array_parser,
}

class ProdyNMDParser(object):
    def __init__(self):
        pass
    
    @classmethod
    def read(cls, nmd_file):
        header = {}
        eigenvalues = []
        eigenvectors = []
        for line in open(nmd_file):
            parts = line.split()
            tag = parts[0]
            if tag in line_parsers:
                line_parser = line_parsers[tag]
                if tag == "mode":
                    eigenvalue, eigenvector = line_parser(parts[1:])
                    eigenvalues.append(eigenvalue)
                    eigenvectors.append(eigenvector)
                else:
                    header[tag] = line_parser(parts[1:])
        return numpy.array(eigenvalues), numpy.array(eigenvectors), header
        