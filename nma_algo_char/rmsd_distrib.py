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


if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("-n", type  = "int", dest="num_procs")
    parser.add_option("-t", type = "float", dest="temperature")
    
    (options, args) = parser.parse_args()
    
    base_folders = {
#                     "/home/victor/Desktop/Desktop/SRCKIN/CC":"CC",
#                     "/home/victor/Desktop/Desktop/SRCKIN/method1":"IC"
                    "/home/victor/Desktop/Desktop/UBI/CC":"CC",
                    "/home/victor/Desktop/Desktop/UBI/IC_method1":"IC"
    }
    
    colors = sns.color_palette("hls", 2)
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
        
        if sim_type == "CC" :
            label = "CC ANM step"
#             data_range = (-1.5, 39)
#             data_range = (-2, 60)#ubi
            color = colors[1]
        else:
            label = "IC"
#             data_range = (-140, 177)
#             data_range = (-65, 400)#ubi
            color = colors[0]
        
        print label, numpy.mean(rmsds), numpy.std(rmsds),numpy.mean(rmsds[who_is_accepted]), numpy.std(rmsds[who_is_accepted])
#         plt.hist(rmsds[who_is_accepted],
#                  bins =100, 
# #                  range = data_range, 
#                  label = label, 
#                  alpha = 0.75, 
#                  normed=True,
#                  color = color)
#     plt.legend()
#     plt.xlabel("Step RMSD ($\AA$) ")
#     plt.legend()
#     plt.show()
    