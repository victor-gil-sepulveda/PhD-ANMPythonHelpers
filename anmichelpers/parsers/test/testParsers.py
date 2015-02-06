"""
Created on Feb 5, 2015

@author: victor
"""
import unittest
from anmichelpers.parsers.imods import ImodServerFilesParser
import os
import anmichelpers.parsers.test

class TestParsers(unittest.TestCase):

    def test_imods(self):
        file_path = os.path.join(anmichelpers.parsers.test.__path__[0], "data","imods1.evec")
        eigenvalues, eigenvectors = ImodServerFilesParser.read(file_path)
        print eigenvalues, eigenvectors


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()