"""
Created on 20/11/2015

@author: victor
"""
import os
from anmichelpers.parsers.pronmd import ProdyNMDParser
from anmichelpers.comparison.comparison import rmsip, cumulative_overlap,\
    degree_of_collectivity
import pandas as pd

if __name__ == '__main__':
    data_folder = "/home/victor/Desktop/ICvsCC_movement"
    proteins =  ["1ddt",  "1ex6",  "1ggg",  "1su4", "4ake", "src_kin", "src_kin_2",]
    file_name = {"cc":"original_modes.1.nmd",
                 "ic":"original_modes_cc.1.nmd"}
    
    all_overlaps = {}
    all_rmsip = {}
    all_doc_cc = {}
    all_doc_ic = {}
    for protein in proteins:
        _, cc_eigenvectors, _ = ProdyNMDParser.read(os.path.join(data_folder,protein, "cc",file_name["cc"])) 
        _, ic_eigenvectors, _ = ProdyNMDParser.read(os.path.join(data_folder,protein, "ic",file_name["ic"])) 
        print cc_eigenvectors.shape,ic_eigenvectors.shape
        
        cum_overlaps = {}
        deg_of_collectivity_cc = {}
        deg_of_collectivity_ic = {}
        for i,mode in enumerate(ic_eigenvectors):
            cum_overlaps[i]  = cumulative_overlap(ic_eigenvectors[i],cc_eigenvectors)
            deg_of_collectivity_ic[i] = degree_of_collectivity(ic_eigenvectors[i], normalize=True)
            deg_of_collectivity_cc[i] = degree_of_collectivity(cc_eigenvectors[i], normalize=True)
            
            
        all_overlaps[protein] = cum_overlaps
        all_rmsip[protein] = [rmsip(cc_eigenvectors, ic_eigenvectors)]
        all_doc_cc[protein] = deg_of_collectivity_cc
        all_doc_ic[protein] = deg_of_collectivity_ic
    
    pd.DataFrame.from_dict(all_overlaps).to_csv("overlaps")
    pd.DataFrame.from_dict(all_rmsip).to_csv("rmsip")
    pd.DataFrame.from_dict(all_doc_cc).to_csv("doc_cc")
    pd.DataFrame.from_dict(all_doc_ic).to_csv("doc_ic")

    #handler = open("rmsip.","w")