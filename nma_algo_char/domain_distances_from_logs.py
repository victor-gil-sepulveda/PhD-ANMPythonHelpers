'''
Created on Jan 28, 2016

@author: victor
'''
from nma_algo_char.common import MetropolisMCSimulator, prepare_subplots
from nma_algo_char.acceptance_and_rmsf_from_logs import load_single_proc_data
import os.path
from nma_algo_char.mode_application_analysis import process_energy_differences
import numpy
from anmichelpers.tools.tools import norm
from optparse import OptionParser
import matplotlib.pyplot as plt
import seaborn as sns

def calc_distances(in_coords, who_is_accepted):
    coords = in_coords[who_is_accepted]
    distances = []
    for coordset in coords:
        distances.append(norm(coordset[128]-coordset[18]))
    return numpy.array(distances)

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("--type", dest="sim_type")
    parser.add_option("-t", type="int", dest="temperature")
    (options, args) = parser.parse_args()
    
    folders = [
             "IC_dispFact_0.65_dm_1",
             "IC_dispFact_0.65_dm_2",
             "IC_dispFact_0.65_dm_3",
             "IC_dispFact_0.65_dm_4",
             "IC_dispFact_0.65_dm_5",
             "IC_dispFact_0.65_dm_6",
             "IC_dispFact_0.65_dm_7",
             "IC_dispFact_0.65_dm_8",
             "IC_dispFact_0.65_dm_9",
             "IC_dispFact_0.65_dm_10"
             ]
    
    distances = {}
    for folder in folders:
        raw_data, min_len = load_single_proc_data(options.sim_type, os.path.join(folder,"info"))
        
        energy_increments = process_energy_differences(raw_data)
        mc = MetropolisMCSimulator(energy_increments)
        
        who_is_accepted = mc.who_is_accepted(options.temperature)
        
        coords = numpy.reshape(raw_data["coords_after"], (len(raw_data["coords_after"]), 
                                                                           len(raw_data["coords_after"][0])/3, 
                                                                           3))
        distances[folder] =  calc_distances(coords, range(len(coords)))#who_is_accepted[:150])
    
    
    sns.set_style("whitegrid")
    row_len = 4
    col_len = 3
    folders.extend(["IC_dispFact_0.65_dm_10","IC_dispFact_0.65_dm_10"])
    f, axes = prepare_subplots(row_len, col_len)
    for i,folder in enumerate(folders):
        ax = axes[i/row_len, i%row_len]
        ax.set_title(folder)
        if i%row_len == 0:
            ax.set_ylabel("Distance ($\AA$)")
        if i/row_len == col_len-1:
            ax.set_xlabel("Step")
        ax.plot(distances[folder])
    plt.suptitle("Inter domain distance")                   
    plt.show()
    
    