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
from nma_algo_char.common import load_control_json, pair_parameter_values,\
    parameter_value_to_string, create_directory
from _collections import defaultdict
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("--plot", action = "store_true", default= False, dest="do_plots")
    
    (options, args) = parser.parse_args()
    
    if not options.experiment:  
        parser.error('Experiment file not given')
    
    if not options.results_folder:  
        parser.error('Results folder not given')
    
    create_directory(options.results_folder)
    
    proteins =  [
#                  "1ddt.fixed.pdb", 
#                  "1ex6.fixed.pdb", 
#                  "1ggg.fixed.pdb", 
#                  #"1su4.fixed.pdb", 
#                  "4ake.fixed.pdb", 
#                  "src_kin.fixed.pdb", 
#                  "src_kin2.fixed.pdb", 
#                  "2lzm.fixed.pdb"
                    "2F4J.pdb"
                 ]
    
    sizes = [parsePDB(os.path.join("structs", pdb_file)).numResidues() for pdb_file in proteins]
    size_per_protein = dict(zip(proteins, sizes)) 
    # order proteins per size
    size_ordered_proteins = [s[1] for s in sorted(zip(sizes,proteins))]
    
    nmd_file_name = {"CC":"normalized_modes.1.nmd", "IC":"normalized_modes_cc.1.nmd"}
    prefix_to_folder = {"CC":"cc", "IC":"ic"}
        
    # Load data
    experiment_details = load_control_json(options.experiment)   
    eigenvectors = defaultdict(lambda: defaultdict())
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        for prefix in ["CC","IC"]:
            
            folder_name = "%s_%s_%s_%s_%s"%(prefix,
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
            
            rev_folder_name = "%s_%s_%s_%s_%s"%(prefix,
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2),
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1))
            
            nmd_file_path = os.path.join(prefix_to_folder[prefix], folder_name, "info",  nmd_file_name[prefix])
            rev_nmd_file_path = os.path.join(prefix_to_folder[prefix], rev_folder_name, "info",  nmd_file_name[prefix])
            
            try:
                _, _eigenvectors, _ = ProdyNMDParser.read(nmd_file_path)
                eigenvectors[prefix][(v1,v2)] = _eigenvectors
            except:
                try:
                    _, _eigenvectors, _ = ProdyNMDParser.read(rev_nmd_file_path)
                    eigenvectors[prefix][(v1,v2)] = _eigenvectors
                except IOError, e:
                    print os.path.join(rev_nmd_file_path), "NOT FOUND"
                    eigenvectors[prefix][(v1,v2)] = None
    
    # Do stuff with data
    # Max_eigen = 30
    MAX_EIGEN = 30
    for eig_key in eigenvectors["CC"]:
        all_ic_cc_overlaps = {}
        all_cc_ic_overlaps = {}
        all_rmsip = {}
        all_doc_cc = {}
        all_doc_ic = {}
        
        cc_eigenvectors = eigenvectors["CC"][eig_key]
        ic_eigenvectors = eigenvectors["IC"][eig_key]
        
        if cc_eigenvectors is not None and ic_eigenvectors is not None:
            cum_overlaps_ic_explained_by_cc = {}
            cum_overlaps_cc_explained_by_ic = {}
            deg_of_collectivity_cc = []
            deg_of_collectivity_ic = []
            for i,mode in enumerate(ic_eigenvectors):
                cum_overlaps_ic_explained_by_cc[i]  = cumulative_overlap(ic_eigenvectors[i],cc_eigenvectors[0:10])
                cum_overlaps_cc_explained_by_ic[i]  = cumulative_overlap(cc_eigenvectors[i],ic_eigenvectors[0:10])
                deg_of_collectivity_ic.append(degree_of_collectivity(ic_eigenvectors[i], normalize=True))
                deg_of_collectivity_cc.append(degree_of_collectivity(cc_eigenvectors[i], normalize=True))
                 
            all_ic_cc_overlaps[eig_key[0]] = cum_overlaps_ic_explained_by_cc
            all_cc_ic_overlaps[eig_key[0]] = cum_overlaps_cc_explained_by_ic
            all_rmsip[eig_key[0]] = [(last_freq, rmsip(cc_eigenvectors[0:last_freq], ic_eigenvectors[0:last_freq])) for last_freq in range(5, MAX_EIGEN+1,5)]
            all_doc_cc[eig_key[0]] = deg_of_collectivity_cc
            all_doc_ic[eig_key[0]] = deg_of_collectivity_ic
     
            if options.do_plots:
                import seaborn as sns
    
                cutt_res_folder = os.path.join(options.results_folder,"cutt_%f"%( eig_key[1]))
                create_directory(cutt_res_folder)
                
                # Plot rmsip for some mode ranges limit
                for protein in size_ordered_proteins:
                    if protein in all_rmsip:
                        x = [] # <- rmsip per protein
                        y = [] # <-frequencies used in rmsip (must be the same for all of them)
                        for last_freq,rmsip_val in all_rmsip[protein]:
                            x.append(last_freq)
                            y.append(rmsip_val)
                        plt.plot(x, y, label = "%s (%d)"%(protein,size_per_protein[protein]))
                lgd = plt.legend()
                plt.tight_layout()
                plt.savefig(os.path.join(cutt_res_folder,"rmsip_per_cutoff.svg"), bbox_extra_artists=(lgd,), bbox_inches='tight')
                
                # plot collectivities
        
                def reform_data(data, cutoff):
                    reformed_data = {"doc":[],"protein":[],"mode":[]}
                    for protein in size_ordered_proteins:
                        if protein in data:
                            for mode in range(cutoff):
                                reformed_data["protein"].append(protein)
                                reformed_data["mode"].append(mode)
                                reformed_data["doc"].append(data[protein][mode])
                    return reformed_data
                 
                for prefix, data, pos in [("cc", all_doc_cc,(0,0)), ("ic", all_doc_ic,(1,0))]:
                    ax = plt.subplot2grid((2,1),pos)
                    reformed_data = reform_data(data, cutoff = 10)
                    pd_data = pd.DataFrame.from_dict(reformed_data)
                    ax = sns.barplot(x="protein", y="doc", hue="mode", data=pd_data)
                    for item in ax.get_xticklabels():
                        item.set_rotation(30)
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
                plt.tight_layout()
                plt.savefig(os.path.join(cutt_res_folder,"collectivity_of_modes.svg"),
                            bbox_extra_artists=(lgd,), 
                            bbox_inches='tight')
          
        #         plt.show()
                
                # Average collectivity for 10 modes, over different cutoffs
                for prefix, data, pos in [("cc", all_doc_cc,(0,0)), ("ic", all_doc_ic,(0,1))]:
                    ax = plt.subplot2grid((1,2),pos)
                    for protein in size_ordered_proteins:
                        x = [] # <- frequencies used 
                        y = [] # <- avg. collectivity per protein
                        y_err = []
                        for cutoff in range(5, MAX_EIGEN+1, 5):
                            x.append(cutoff)
                            y.append(numpy.mean(data[protein][0:cutoff]))
                            y_err.append(numpy.std(data[protein][0:cutoff]))
                        ax = plt.plot(x, y, label = "%s (%d)"%(protein,size_per_protein[protein]))
                        #plt.errorbar(x, y, yerr= y_err) #demasiado grandes
                        plt.ylim((0.15,0.7))
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
                plt.tight_layout()
                plt.savefig(os.path.join(cutt_res_folder,"avg_collectivity.svg"), 
                                    bbox_extra_artists=(lgd,), 
                                    bbox_inches='tight')
                    
                # Cumulative overlap plot for each mode and protein and calculation type (2 plots, modes to 10)
                def reform_data_for_overlap(cum_overlap):
                    data = {"protein":[], "mode":[], "cum. overlap":[]}
                    for protein in size_ordered_proteins:
                        for mode in range(20):
                            data["protein"].append(protein)
                            data["mode"].append(mode)
                            data["cum. overlap"].append(cum_overlap[protein][mode])
                    return data        
                    
                for prefix, data, pos in [("cc", all_cc_ic_overlaps,(0,0)), ("ic", all_ic_cc_overlaps,(1,0))]:
                    ax = plt.subplot2grid((2,1),pos)
                    dict_data  = reform_data_for_overlap(data)
                    df = pd.DataFrame.from_dict(dict_data)
                    sns.barplot("protein","cum. overlap", "mode", df)
                    avgs = []
                    for protein in size_ordered_proteins:
                        co_p_mode = []
                        for mode in range(20): 
                            co_p_mode.append(data[protein][mode])
                        avgs.append(numpy.mean(co_p_mode))
                    plt.plot(range(len(set(dict_data["protein"]))), avgs, marker="o", color="black")
                    plt.ylim((0.0,1.0))
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
                plt.tight_layout()
                plt.savefig(os.path.join(cutt_res_folder,"cumulative_overlap.svg"), 
                                    bbox_extra_artists=(lgd,), 
                                    bbox_inches='tight')
                    
