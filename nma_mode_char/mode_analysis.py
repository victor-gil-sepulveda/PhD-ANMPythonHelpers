"""
Created on 20/11/2015

@author: victor
"""
import os
import numpy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
from prody.proteins.pdbfile import parsePDB
from _collections import defaultdict
from anmichelpers.parsers.pronmd import ProdyNMDParser
from anmichelpers.comparison.comparison import rmsip, cumulative_overlap,\
    degree_of_collectivity
from nma_algo_char.common import load_control_json, pair_parameter_values,\
    parameter_value_to_string, create_directory
import anmichelpers.tools.tools as tools

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--results", dest="results_folder")
    
    (options, args) = parser.parse_args()
    
    if not options.experiment:  
        parser.error('Experiment file not given')
    
    if not options.results_folder:  
        parser.error('Results folder not given')
    
    create_directory(options.results_folder)
    
    proteins =  [
                    "1ddt.fixed.pdb", 
                    "1ex6.fixed.pdb", 
                    "1ggg.fixed.pdb", 
                    "4ake.fixed.pdb", 
                    "src_kin.fixed.pdb", 
                    "src_kin2.fixed.pdb", 
                    "2lzm.fixed.pdb",
                    "ubi_cut.fixed.pdb",
                    "ubi_start.pdb"
                 ]
    protein_ids = {
                    "1ddt.fixed.pdb": "1ddt", 
                    "1ex6.fixed.pdb": "1ex6", 
                    "1ggg.fixed.pdb": "1ggg", 
                    "4ake.fixed.pdb": "4ake", 
                    "src_kin.fixed.pdb": "1y57", 
                    "src_kin2.fixed.pdb": "1y57_MD", 
                    "2lzm.fixed.pdb":"2lzm",
                    "ubi_cut.fixed.pdb":"1ubq_cut",
                    "ubi_start.pdb":"1ubq"
                  }
    
    sizes = [parsePDB(os.path.join("structs", pdb_file)).numResidues() for pdb_file in proteins]
    size_per_protein = dict(zip(proteins, sizes)) 
    # order proteins per size
    size_ordered_proteins = [s[1] for s in sorted(zip(sizes,proteins))]
    
    nmd_file_name = {"CC":"normalized_modes.1.nmd", "IC":"normalized_modes_cc.1.nmd"}
    prefix_to_folder = {"CC":"cc", "IC":"ic"}
        
    # Load data
    experiment_details = load_control_json(options.experiment)   
    eigenvectors = defaultdict(lambda: defaultdict())
    
    prot_keys = []
    cutoff_keys= []
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        for prefix in ["CC","IC"]:
            if p1 =="prot": 
                prot_keys.append(v1)
                cutoff_keys.append(v2)
                key = (v1,v2)
            else:
                prot_keys.append(v2)
                cutoff_keys.append(v1)
                key = (v2,v1)
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
                eigenvectors[prefix][key] = _eigenvectors
            except:
                try:
                    _, _eigenvectors, _ = ProdyNMDParser.read(rev_nmd_file_path)
                    eigenvectors[prefix][key] = _eigenvectors
                except IOError, e:
                    print os.path.join(rev_nmd_file_path), "NOT FOUND"
                    eigenvectors[prefix][key] = None
                    
            ##
            ## OJO! ESTA NORMALIZACION NO SE TENDRIA QUE HACER AQUI!!! DEBERIA VENIR DE SERIE!!
            ##
            for mode in range(len(eigenvectors[prefix][key])):
                norm = tools.norm(numpy.array(eigenvectors[prefix][key][mode]))
                eigenvectors[prefix][key][mode] = numpy.array(eigenvectors[prefix][key][mode]) / norm
                print tools.norm(eigenvectors[prefix][key][mode])
    
    MAX_EIGEN = 30
    NUM_MODES_TO_EXPLAIN_ANOTHER_MODE_CO = 10 # for cumulative overlap
    sns.set_style("whitegrid")
    
    avg_cc_collectivity_per_cutoff = defaultdict(dict)
    avg_ic_collectivity_per_cutoff = defaultdict(dict)
    avg_cc_cum_overlap_per_cutoff = defaultdict(dict)
    avg_ic_cum_overlap_per_cutoff = defaultdict(dict)
    for cutoff in set(cutoff_keys):
        cutt_res_folder = os.path.join(options.results_folder,"cutt_%.4f"%(cutoff))
        create_directory(cutt_res_folder)
        
        all_rmsip = {}
        all_doc_cc = {}
        all_doc_ic = {}
        all_cc_ic_overlaps = {}
        all_ic_cc_overlaps = {}
        for protein in size_ordered_proteins:
            key = (protein, cutoff)
            cc_eigenvectors = eigenvectors["CC"][key]
            ic_eigenvectors = eigenvectors["IC"][key]
            if cc_eigenvectors is not None and ic_eigenvectors is not None:
                cum_overlaps_ic_explained_by_cc = {}
                cum_overlaps_cc_explained_by_ic = {}
                deg_of_collectivity_cc = []
                deg_of_collectivity_ic = []
                for i in range(NUM_MODES_TO_EXPLAIN_ANOTHER_MODE_CO):
                    print protein, key, i, tools.norm(ic_eigenvectors[i]), tools.norm(cc_eigenvectors[i])
                    cum_overlaps_ic_explained_by_cc[i]  = cumulative_overlap(ic_eigenvectors[i],cc_eigenvectors[:NUM_MODES_TO_EXPLAIN_ANOTHER_MODE_CO])
                    cum_overlaps_cc_explained_by_ic[i]  = cumulative_overlap(cc_eigenvectors[i],ic_eigenvectors[:NUM_MODES_TO_EXPLAIN_ANOTHER_MODE_CO])
                    if protein == "ubi_start.pdb" and i == 1 and cutoff == 9:
                        print "***"
                        cumulative_overlap(cc_eigenvectors[i],
                                           ic_eigenvectors[:NUM_MODES_TO_EXPLAIN_ANOTHER_MODE_CO],
                                           verbose = True)
                
                for i in range(MAX_EIGEN):
                    deg_of_collectivity_ic.append(degree_of_collectivity(ic_eigenvectors[i], normalize=True))
                    deg_of_collectivity_cc.append(degree_of_collectivity(cc_eigenvectors[i], normalize=True))
                all_ic_cc_overlaps[protein] = cum_overlaps_ic_explained_by_cc
                all_cc_ic_overlaps[protein] = cum_overlaps_cc_explained_by_ic
                all_rmsip[protein] = [(last_freq, rmsip(cc_eigenvectors[0:last_freq], ic_eigenvectors[0:last_freq])) for last_freq in range(5, MAX_EIGEN+1,5)]
                all_doc_cc[protein] = deg_of_collectivity_cc
                all_doc_ic[protein] = deg_of_collectivity_ic
            avg_cc_collectivity_per_cutoff[protein][cutoff] = numpy.mean(deg_of_collectivity_cc)
            avg_ic_collectivity_per_cutoff[protein][cutoff] = numpy.mean(deg_of_collectivity_ic)
            avg_cc_cum_overlap_per_cutoff[protein][cutoff] = numpy.mean([cum_overlaps_cc_explained_by_ic[k] for k in cum_overlaps_cc_explained_by_ic])
            avg_ic_cum_overlap_per_cutoff[protein][cutoff] = numpy.mean([cum_overlaps_ic_explained_by_cc[k] for k in cum_overlaps_ic_explained_by_cc])
        
        #------------------------------
        # RMSIP
        #------------------------------
        # The space overlap between the modes (first 5, first 10 .. etc)
        with sns.color_palette("hls", len(protein)):                
            # Plot rmsip for some mode ranges limit
            for protein in size_ordered_proteins:
                if protein in all_rmsip:
                    x = [] # <- rmsip per protein
                    y = [] # <- frequencies used in rmsip (must be the same for all of them)
                    for last_freq, rmsip_val in all_rmsip[protein]:
                        x.append(last_freq)
                        y.append(rmsip_val)
                    plt.plot(x, y, label = "%s (%d)"%(protein_ids[protein], size_per_protein[protein]))
            lgd = plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(cutt_res_folder,"rmsip_per_cutoff.svg"), bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.title("RMSIP")
            plt.close()            
   
        #------------------------------
        # Collectivity
        #------------------------------
        def reform_data(data, cut):
            reformed_data = {"DOC":[],"Protein":[],"Mode":[]}
            for protein in size_ordered_proteins:
                if protein in data:
                    for mode in range(cut):
                        reformed_data["Protein"].append(protein_ids[protein])
                        reformed_data["Mode"].append(mode)
                        reformed_data["DOC"].append(data[protein][mode])
            return reformed_data
          
        with sns.color_palette("BuGn_r",10):
            for prefix, data, pos in [("CC", all_doc_cc,(0,0)), ("IC", all_doc_ic,(1,0))]:
                ax = plt.subplot2grid((2,1),pos)
                ax.set_title(prefix)
                reformed_data = reform_data(data, cut = 10)
                pd_data = pd.DataFrame.from_dict(reformed_data)
                ax = sns.barplot(x="Protein", y="DOC", hue="Mode", data=pd_data)
                avgs = []
                errors = []
                for protein in size_ordered_proteins:
                    co_p_mode = []
                    modes_for_stats = 10
                    for mode in range(modes_for_stats): 
                        co_p_mode.append(data[protein][mode])
                    avgs.append(numpy.mean(co_p_mode[:modes_for_stats]))
                    errors.append(numpy.std(co_p_mode[:modes_for_stats]))
                ax.errorbar(range(len(avgs)), avgs, yerr = errors, color="black")
                plt.plot(range(len(avgs)), avgs, marker="o", color="black")
                for item in ax.get_xticklabels():
                    item.set_rotation(30)
            plt.tight_layout()
            lgd = plt.legend(loc='upper left')
              
            plt.savefig(os.path.join(cutt_res_folder,"collectivity_of_modes.svg"),
                        bbox_extra_artists=(lgd,), 
                        bbox_inches='tight')
            plt.close()
          
  
        #------------------------------
        # Cumulative overlap
        #------------------------------
        # Cumulative overlap plot for each mode and protein and calculation type (2 plots, modes to 10)
        def reform_data_for_overlap(cum_overlap):
            data = {"Protein":[], "Mode":[], "Cum. Overlap":[]}
            for protein in size_ordered_proteins:
                for mode in range(NUM_MODES_TO_EXPLAIN_ANOTHER_MODE_CO):
                    data["Protein"].append(protein_ids[protein])
                    data["Mode"].append(mode)
                    data["Cum. Overlap"].append(cum_overlap[protein][mode])
            return data        
          
        with sns.color_palette("BuGn_r",10):  
            for prefix, data, pos in [("CC", all_cc_ic_overlaps,(0,0)), ("IC", all_ic_cc_overlaps,(1,0))]:
                ax = plt.subplot2grid((2,1),pos)
                ax.set_title(prefix)
                dict_data  = reform_data_for_overlap(data)
                df = pd.DataFrame.from_dict(dict_data)
                sns.barplot("Protein","Cum. Overlap", "Mode", df)
                avgs = []
                errors = []
                for protein in size_ordered_proteins:
                    co_p_mode = []
                    modes_for_stats = 10
                    for mode in range(modes_for_stats): 
                        co_p_mode.append(data[protein][mode])
                    avgs.append(numpy.mean(co_p_mode[:modes_for_stats]))
                    errors.append(numpy.std(co_p_mode[:modes_for_stats]))
                plt.errorbar(range(len(set(dict_data["Protein"]))), avgs, yerr = errors, color="black")
                plt.plot(range(len(set(dict_data["Protein"]))), avgs, marker="o", color="black")
                plt.ylim((0.0,1.0))
                for item in ax.get_xticklabels():
                    item.set_rotation(30)
            plt.tight_layout()
            lgd = plt.legend(loc='upper left')
            plt.savefig(os.path.join(cutt_res_folder,"cumulative_overlap.svg"), 
                                bbox_extra_artists=(lgd,), 
                                bbox_inches='tight')
            plt.close()
    with sns.color_palette("hls", len(proteins)):
        for  ylab, title, cc_ic_data in [("Degree of collectivity","avg_collectivities_per_cutoff.svg",(avg_cc_collectivity_per_cutoff,avg_ic_collectivity_per_cutoff)),
                                   ("Cumulative overlap","avg_cum_overlap_per_cutoff.svg",(avg_cc_cum_overlap_per_cutoff, avg_ic_cum_overlap_per_cutoff))]:
            for prefix, pos, data in [("CC", (0,0),cc_ic_data[0]),
                                      ("IC", (0,1), cc_ic_data[1])]:
                ax = plt.subplot2grid((1,2), pos)
                ax.set_title(prefix)
                for protein in size_ordered_proteins:
                    x = []
                    y = []
                    for cutoff in sorted(set(cutoff_keys)):         
                        x.append(cutoff)
                        y.append(data[protein][cutoff])  
                    ax.plot(x,y, label = "%s (%d)"%(protein_ids[protein],size_per_protein[protein]))
                    ax.set_ylim((0.0,1.0))
                    ax.set_xlabel("Cutoff ($\AA$)")
                    ax.set_ylabel(ylab)
            plt.tight_layout()
            lgd = plt.legend()
            plt.savefig(os.path.join(options.results_folder,title),
                                bbox_extra_artists=(lgd,), 
                                bbox_inches='tight')
            plt.close()
        