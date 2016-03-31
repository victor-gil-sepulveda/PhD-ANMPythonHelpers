"""
Created on Mar 13, 2016

@author: victor
"""
import seaborn as sns
from nma_algo_char.data_retrieval import load_data_from_multiple_processors,\
    load_single_proc_data, process_energy_differences,\
    process_after_perturb_rmsd
from optparse import OptionParser
import matplotlib.pyplot as plt
import os.path
import numpy
from nma_algo_char.common import MetropolisMCSimulator

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

def calc_ratio(energy_increments):
    positives = numpy.count_nonzero(energy_increments>0.)
    negatives = numpy.count_nonzero(energy_increments<0.)
    return float(negatives) / positives

if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("--results", dest="results_folder")
    parser.add_option("-n", type  = "int", dest="num_procs")
    
    (options, args) = parser.parse_args()
    
#     create_directory(options.results_folder)
    
    base_folders = {
#                     "/home/victor/Desktop/Desktop/SRCKIN/CC":"CC",
#                     "/home/victor/Desktop/Desktop/SRCKIN/method1":"IC"
                    "/home/victor/Desktop/Desktop/UBI/CC":"CC",
                    "/home/victor/Desktop/Desktop/UBI/IC_method1":"IC"
    }
    
    f, axes = plt.subplots(3, 1)
    colors = sns.color_palette("hls", 3)
    for folder in base_folders:
        sim_type = base_folders[folder]
        data_folder = "info"
        
        raw_data, min_len = load_the_data(options.num_procs, 
                                          sim_type, 
                                          False, 
                                          os.path.join(folder, data_folder))
        
        energy_increments = process_energy_differences(raw_data)
        
        if sim_type == "CC" :
            label = "CC ANM step"
            data_range = (-1.5, 39)
            data_range = (-2, 60)#ubi
            axis = axes[1]
            color = colors[1]
        else:
            label = "IC"
            data_range = (-140, 177)
            data_range = (-65, 400)#ubi
            axis = axes[0]
            color = colors[0]

        print label,  numpy.mean(energy_increments), numpy.std(energy_increments)
        print "Ratio", calc_ratio(energy_increments)
        print "Acceptance", MetropolisMCSimulator(energy_increments).perform_simulation(5000,200,300)
        
        axis.hist(energy_increments,
                  #energy_increments[energy_increments < numpy.mean(energy_increments)+(numpy.std(energy_increments)*3)], 
                 bins =100, 
                 range = data_range, 
                 label = label, 
                 alpha = 0.75, 
                 normed=True,
                 color = color)
        axis.legend()
#         plt.xlabel("Energy increment (kcal / mol)")
#         plt.legend()
#         plt.show()
    
        # Then calculate the "full pele stuff too"
        if sim_type == "CC" :
            raw_data, min_len = load_the_data(options.num_procs, 
                                              sim_type, 
                                              True, 
                                              os.path.join(folder, data_folder))
             
            energy_increments = process_energy_differences(raw_data)
            print "CC_FULL",  numpy.mean(energy_increments), numpy.std(energy_increments)
            print label,  numpy.mean(energy_increments), numpy.std(energy_increments)
            print "Ratio", calc_ratio(energy_increments)
            print "Acceptance", MetropolisMCSimulator(energy_increments).perform_simulation(5000,200,300)
        
            data_range = (-10.5, 10.5)
            data_range = (-2.1, 6) #ubi
            axis = axes[2]
            color = colors[2]
            axis.hist(energy_increments, bins = 100, 
                     range = data_range, 
                     label = "CC Full step", 
                     alpha = 0.75 , 
                     normed=True,
                     color = color)
            axis.legend()
#             plt.legend()
#             plt.xlabel("Energy increment (kcal / mol)")
#             plt.show()
    plt.legend()
    plt.xlabel("Energy increment (kcal / mol)")
    plt.show()
    
    