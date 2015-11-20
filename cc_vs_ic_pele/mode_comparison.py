"""
Created on 20/11/2015

@author: victor
"""
import os
from anmichelpers.parsers.pronmd import ProdyNMDParser
from anmichelpers.comparison.comparison import rmsip

if __name__ == '__main__':
    data_folder = "/home/victor/VBShared/ICvsCC_movement"
    proteins =  ["1ddt",  "1ex6",  "1ggg",  "1su4", "4ake", "src_kin"]
    file_name = {"cc":"original_modes.1.nmd",
                 "ic":"original_modes_cc.1.nmd"}
    
    for protein in proteins:
        _, cc_eigenvectors, _ = ProdyNMDParser.read(os.path.join(data_folder,protein, "cc",file_name["cc"])) 
        _, ic_eigenvectors, _ = ProdyNMDParser.read(os.path.join(data_folder,protein, "ic",file_name["ic"])) 
        print cc_eigenvectors.shape,ic_eigenvectors.shape
        
        print rmsip(cc_eigenvectors, ic_eigenvectors)
    