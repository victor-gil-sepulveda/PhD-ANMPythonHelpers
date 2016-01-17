'''
Created on 14/12/2015

@author: victor
'''
from optparse import OptionParser
from nma_algo_char.common import load_control_json, create_directory,\
    pair_parameter_values, parameter_value_to_string, MetropolisMCSimulator,\
    get_values_by_hue
import os
from _collections import defaultdict
from nma_algo_char.mode_application_analysis import load_cc_data, load_ic_data,\
    process_energy_differences, process_after_perturb_rmsd, load_data,\
    load_energy
import numpy
import matplotlib.pyplot as plt
import seaborn as sns
import cPickle as pickle
from anmichelpers.tools.measure import coords_rmsf
from trajectory_comparison.compare_two_rmsfs import rms
import scipy.stats

if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--workspace", dest="workspace")
    
    (options, args) = parser.parse_args()
    
    if not options.experiment:  
        parser.error('Experiment file not given')
    
    
    experiment_details = load_control_json(options.experiment)
    if options.workspace is not None:
        workspace = os.path.normpath(options.workspace)
    else:
        workspace = os.path.normpath(experiment_details["workspace"])
    
    # Things to load
    all_anm_energies = {}
    all_final_energies = {}
    
    ## ONLY FOR CC SIMULATIONS
    temperatures = [300, 583, 866, 1150, 1432, 2000, 2568, 300, 300]
    for T in temperatures:
        all_anm_energies[T] = []
        all_final_energies[T] = []
        main_folder = "%s_%d"%(experiment_details["prefix"],T)
        print "Processing %s"%main_folder
        for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
            
            folder_name = "%s_%s_%s_%s_%s"%("CC",
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
            
            data_folder = os.path.join(main_folder, workspace, folder_name,"info")

            beg_ener = load_energy(data_folder, "perturb_energy_before.log")
            anm_ener = load_energy(data_folder, "perturb_energy_after.log")
            fin_ener = load_energy(data_folder, "final_energy.log")
            
            min_len = min(len(anm_ener), len(fin_ener), len(beg_ener))
            
            all_anm_energies[T].extend(anm_ener[:min_len]-beg_ener[:min_len])
            all_final_energies[T].extend(fin_ener[:min_len]-beg_ener[:min_len])
        numpy.savetxt("%s_final.txt"%main_folder, numpy.array(all_final_energies[T]).T)
        numpy.savetxt("%s_anm.txt"%main_folder, numpy.array(all_anm_energies[T]).T)
        
    
    def prepare_subplots(row_len, col_len):
        if row_len > 1 or col_len > 1:
            f, axes = plt.subplots( col_len, row_len, sharey='row', sharex='col')
            f.subplots_adjust(hspace=0.4, wspace=0.3 )
            f.set_size_inches(12, 12, forward=True)
        else:
            f = plt.gcf()
            axes = {(0,0): plt.gca()}
        return f, axes
    
    def filter_outliers(ener1,ener2):
        filter_1 = []
        filter_2 = []
        for i in range(len(ener1)):
            if ener1[i]< -80 and ener2[i]<150:
                pass
            else:
                filter_1.append(ener1[i])
                filter_2.append(ener2[i])
        return numpy.array(filter_1), numpy.array(filter_2)
      
#     row_len = 3
#     col_len = 3
#     f, axes = prepare_subplots(row_len, col_len)
#     for i,T in enumerate(temperatures):
#         ax = axes[i/row_len, i%row_len]
# #         ax = plt.gca()
#         x , y = filter_outliers(all_anm_energies[T], all_final_energies[T])
#         ax.scatter(x, y, marker=',', c = "g", alpha = 0.3, lw = 0, s=0.6)
#         ax.set_title("T = %d"%T)
#         ax.set_xlim((-20,80))
#         ax.set_ylim((-40,60))
#         if i%row_len == 0:
#             ax.set_ylabel("Step U inc.")
#         if i/row_len == 2:
#             ax.set_xlabel("ANM U inc.")
#     
#     plt.show()       

    for T in sorted(list(set(temperatures))):
        ranked_anm_energies = scipy.stats.rankdata(numpy.array(all_anm_energies[T]))
        ranked_final_energies = scipy.stats.rankdata(numpy.array(all_final_energies[T]))
        scaled_ranked_anm_energies = ranked_anm_energies*numpy.max(ranked_final_energies)/numpy.max(ranked_anm_energies)
        rho, p_val =  scipy.stats.spearmanr(scaled_ranked_anm_energies, ranked_final_energies)
        print "%d\t%.3f\t%.3f"%(T, rho, p_val)
    
    def filter_by_final_increment(anm_ener, final_ener, T):
        filt_anm_ener = []
        filt_final_ener = []
        threshold_ener = MetropolisMCSimulator.energy_for_probability(0.1, T)
        
        for i in range(len(anm_ener)):
            if final_ener[i]<threshold_ener:
                filt_anm_ener.append(anm_ener[i])
                filt_final_ener.append(final_ener[i])
        
        return filt_anm_ener, filt_final_ener, threshold_ener
    
    for T in sorted(list(set(temperatures))):
        filt_anm_ener, filt_final_ener, threshold_ener = filter_by_final_increment(all_anm_energies[T], all_final_energies[T], T)
        ranked_anm_energies = scipy.stats.rankdata(numpy.array(filt_anm_ener))
        ranked_final_energies = scipy.stats.rankdata(numpy.array(filt_final_ener))
        scaled_ranked_anm_energies = ranked_anm_energies*numpy.max(ranked_final_energies)/numpy.max(ranked_anm_energies)
        rho, p_val =  scipy.stats.spearmanr(scaled_ranked_anm_energies, ranked_final_energies)
        print "%d\t%d\t%.3f\t%.3f\t%.3f"%(T, len(filt_anm_ener), threshold_ener, rho, p_val)
