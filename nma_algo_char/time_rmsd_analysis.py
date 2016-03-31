'''
Created on Mar 12, 2016

@author: victor
'''
from nma_algo_char.common import create_directory, MetropolisMCSimulator
from optparse import OptionParser
import seaborn as sns
from nma_algo_char.acceptance_and_rmsf_from_logs import load_data_from_multiple_processors,\
    load_single_proc_data
from nma_algo_char.data_retrieval import process_after_perturb_rmsd,\
    process_energy_differences
import os.path
import matplotlib.pyplot as plt


def load_the_data(num_procs, sim_type, full_pele_step, folder):
    if num_procs > 1:
        raw_data, min_len = load_data_from_multiple_processors(sim_type, 
                                                               options.num_procs, 
                                                               folder,
                                                               skip_first = 1,
                                                               #max_samples=500,
                                                               full_pele_energy=full_pele_step)
    else:
        raw_data, min_len = load_single_proc_data(sim_type, 
                                                  folder,
                                                  skip_first = 1,
                                                  #max_samples=500,
                                                  full_pele_energy=full_pele_step)
    return raw_data, min_len

if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("--results", dest="results_folder")
    parser.add_option("-t", type = "float", dest="temperature")
    parser.add_option("-n", type  = "int", dest="num_procs")
    
    (options, args) = parser.parse_args()
    
    create_directory(options.results_folder)
    
    base_folders = {
#                     "/home/victor/Desktop/Desktop/SRCKIN/CC":"CC",
#                     "/home/victor/Desktop/Desktop/SRCKIN/method1":"IC"
                    "/home/victor/Desktop/Desktop/UBI/CC":"CC",
                    "/home/victor/Desktop/Desktop/UBI/IC_method1":"IC"
    }
    
    for folder in base_folders:
        sim_type = base_folders[folder]
        data_folder = "info"
        
        raw_data, min_len = load_the_data(options.num_procs, 
                                          sim_type, 
                                          False, 
                                          os.path.join(folder, data_folder))
        
        energy_increments = process_energy_differences(raw_data)
        mc = MetropolisMCSimulator(energy_increments)
        who_is_accepted = mc.who_is_accepted(options.temperature)
        
        rmsds = process_after_perturb_rmsd(raw_data)
        time = raw_data["time_per_step"]

        if sim_type == "CC" :
            who_is_accepted = range(len(energy_increments))
            label = "CC ANM step (all)"
        else:
            label = "IC NMA step (accepted)"
            
        
        plt.hist((rmsds/time)[who_is_accepted], bins = 70, 
#                  range = [0., 0.11],
                 range = [0., 0.28], #ubi 
                 label = label, alpha = 0.75, normed=True)
        
        if sim_type == "IC":
            label = "IC + Sd. perturb (accepted)"
            plt.hist((rmsds/(time+2.5))[who_is_accepted], bins = 70, 
#                      range = [0., 0.11], 
                        range = [0., 0.28], #ubi
                 label = label, alpha = 0.75, normed=True)
     
        # Then calculate the "full pele stuff too"
        if sim_type == "CC" :
            raw_data, min_len = load_the_data(options.num_procs, 
                                              sim_type, 
                                              True, 
                                              os.path.join(folder, data_folder))
              
            energy_increments = process_energy_differences(raw_data)
            mc = MetropolisMCSimulator(energy_increments)
            who_is_accepted = mc.who_is_accepted(options.temperature)
              
            rmsds = process_after_perturb_rmsd(raw_data)
            time = raw_data["time_per_step"]
              
            plt.hist((rmsds/time)[who_is_accepted], bins = 70, 
#                      range = [0., 0.11], 
                        range = [0., 0.28], #ubi
                     label = "CC Full step (accepted)", alpha = 0.75 , normed=True)
    plt.xlabel("Computational efficiency ($\AA$ / s)")
    plt.legend()
    plt.show()
    
    