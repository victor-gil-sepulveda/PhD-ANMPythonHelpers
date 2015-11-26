"""
Created on 20/11/2015

@author: victor
"""
import os
from anmichelpers.parsers.pronmd import ProdyNMDParser
from anmichelpers.comparison.comparison import rmsip, cumulative_overlap,\
    degree_of_collectivity
import pandas as pd
from prody.proteins.pdbfile import parsePDB
from docutils.frontend import OptionParser
from nma_algo_char.common import load_control_json, pair_parameter_values
from _collections import defaultdict

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("--workspace", dest="workspace")
    parser.add_option("--plot", action = "store_true", default= False, dest="do_plots")
    
    (options, args) = parser.parse_args()
    
    data_folder = "/home/victor/Desktop/ICvsCC_movement"
    proteins =  ["1ddt.fixed.pdb", 
                 "1ex6.fixed.pdb", 
                 "1ggg.fixed.pdb", 
                 "1su4.fixed.pdb", 
                 "4ake.fixed.pdb", 
                 "src_kin.fixed.pdb", 
                 "src_kin_2.fixed.pdb", 
                 "2lzm.fixed.pdb"]
    sizes = [(parsePDB(os.path.join("structs", pdb_file)).numResidues(), pdb_file) for pdb_file in proteins]
    nmd_file_name = {"CC":, "IC":}
        
    all_overlaps = {}
    all_rmsip = {}
    all_doc_cc = {}
    all_doc_ic = {}
    
    # Load data
    experiment_details = load_control_json(options.experiment)   
    eigenvectors = defaultdict(lambda: defaultdict())
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        for prefix in ["CC","IC"]
            folder_name = "%s_%s_%s_%s_%s"%(prefix,
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
            _, _eigenvectors, _ = ProdyNMDParser.read(os.path.join(data_folder, protein, "ic", nmd_file_name[experiment_details["prefix"]])) 
            eigenvectors[experiment_details["prefix"]][(v1,v2)] = _eigenvectors

    # Do stuff with data
    for eig_key in eigenvectors["CC"]:
        cc_eigenvectors = eigenvectors["CC"][eig_key]
        ic_eigenvectors = eigenvectors["IC"][eig_key]

        cum_overlaps = {}
        deg_of_collectivity_cc = {}
        deg_of_collectivity_ic = {}
        for i,mode in enumerate(ic_eigenvectors):
            cum_overlaps[i]  = cumulative_overlap(ic_eigenvectors[i],cc_eigenvectors)
            deg_of_collectivity_ic[i] = degree_of_collectivity(ic_eigenvectors[i], normalize=True)
            deg_of_collectivity_cc[i] = degree_of_collectivity(cc_eigenvectors[i], normalize=True)
            
            
        all_overlaps[eig_key] = cum_overlaps
        all_rmsip[eig_key] = [rmsip(cc_eigenvectors[0:last_freq], ic_eigenvectors[0:last_freq] for last_freq in range(2, 10))
        all_doc_cc[eig_key] = deg_of_collectivity_cc
        all_doc_ic[eig_key] = deg_of_collectivity_ic
    
#     pd.DataFrame.from_dict(all_overlaps).to_csv("overlaps")
#     pd.DataFrame.from_dict(all_rmsip).to_csv("rmsip")
#     pd.DataFrame.from_dict(all_doc_cc).to_csv("doc_cc")
#     pd.DataFrame.from_dict(all_doc_ic).to_csv("doc_ic")

    #handler = open("rmsip.","w")