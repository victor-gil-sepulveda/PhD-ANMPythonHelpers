'''
Created on 20/1/2016

@author: victor
'''
from optparse import OptionParser
from nma_algo_char.common import create_directory, MetropolisMCSimulator
from nma_algo_char.acceptance_and_rmsf_from_logs import load_data_from_multiple_processors,\
    load_single_proc_data
import numpy
import os.path
from nma_algo_char.mode_application_analysis import process_energy_differences
import seaborn as sns
import matplotlib.pyplot as plt
from anmichelpers.tools.tools import norm

def add_coords_data(data, folder, filename, num_procs, min_len):
    # name = "coords_after_anm"
    data[filename] = []
    
    if num_procs > 1:
        for i, proc in enumerate(range(1,num_procs)):
            loaded_data = numpy.delete(numpy.loadtxt(os.path.join(data_folder, "%s.p%d.log"%(filename, proc))),0,1) 
            data[filename].extend(loaded_data[:min_len[i]]) 
    else:
        loaded_data = numpy.delete(numpy.loadtxt(os.path.join(data_folder, "%s.log"%filename )),0,1)
        data[filename].extend(loaded_data[:min_len])
    
    return data

if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("--results", dest="results_folder")
    parser.add_option("-t", type = "int", dest="temperature")
    parser.add_option("-n", type  = "int", dest="num_procs")
    
    (options, args) = parser.parse_args()
    
    create_directory(options.results_folder)
    data_folder = "info"
    if options.num_procs > 1:
        raw_data, min_len = load_data_from_multiple_processors("CC", options.num_procs, data_folder)
    else:
        raw_data, min_len = load_single_proc_data("CC", data_folder)

    raw_data = add_coords_data(raw_data, data_folder, "after_anm_cc", options.num_procs, min_len)
    
    energy_increments = process_energy_differences(raw_data)[1:]
    mc = MetropolisMCSimulator(energy_increments)
    who_is_accepted = mc.who_is_accepted(options.temperature)
    
    initial_coords = numpy.reshape(raw_data["coords_before"], (len(raw_data["coords_before"]), 
                                                                           len(raw_data["coords_before"][0])/3, 
                                                                           3))
    
    anm_coords = numpy.reshape(raw_data["after_anm_cc"], (len(raw_data["after_anm_cc"]), 
                                                                           len(raw_data["after_anm_cc"][0])/3, 
                                                                           3))
    
    minim_coords = numpy.reshape(raw_data["coords_after"], (len(raw_data["coords_after"]), 
                                                                           len(raw_data["coords_after"][0])/3, 
                                                                           3))
    
    def calc_distances(in_coords, who_is_accepted):
        coords = in_coords[who_is_accepted]
        distances = []
        for coordset in coords:
            distances.append(norm(coordset[128]-coordset[18]))
        return numpy.array(distances)
    
    initial_distances = calc_distances(initial_coords, who_is_accepted)
    anm_distances = calc_distances(anm_coords, who_is_accepted) 
    min_distances = calc_distances(minim_coords, who_is_accepted) 
    
    anm_increments = initial_distances - anm_distances
    anm_min_increments = min_distances - anm_distances
    min_increments = initial_distances - min_distances
    plt.plot(anm_increments, label = "after ANM")
    plt.plot(min_increments, label = "after Minim.")
    plt.plot(anm_min_increments, label = "Minim. - ANM")
    plt.legend()
    plt.savefig(os.path.join(options.results_folder,"domain_distances.svg"))
    plt.close()
    
    print "Initial to Anm. negative to positive count ratio:  ",
    positives = numpy.count_nonzero(anm_increments >= 0.)
    negatives = len(anm_increments) - positives
    print float(negatives) / positives
    
    print "Initial to min. negative to positive count ratio:  ",
    positives = numpy.count_nonzero(min_increments >= 0.)
    negatives = len(min_increments) - positives
    print float(negatives) / positives
    
    print "Anm to min. negative to positive count ratio:  ",
    positives = numpy.count_nonzero(anm_min_increments >= 0.)
    negatives = len(anm_min_increments) - positives
    print float(negatives) / positives
    
    print "Initial to Anm. negative to positive sum ratio:  ",
    positives = numpy.sum(anm_increments[anm_min_increments >= 0.])
    negatives = numpy.sum(anm_increments[anm_min_increments < 0.])
    print float(abs(negatives)) / positives
    
    print "Initial to min. negative to positive sum ratio:  ",
    positives = numpy.sum(min_increments[anm_min_increments >= 0.])
    negatives = numpy.sum(min_increments[anm_min_increments < 0.])
    print float(abs(negatives)) / positives
    
    print "Anm to min. negative to positive sum ratio:  ",
    positives = numpy.sum(anm_min_increments[anm_min_increments >= 0.])
    negatives = numpy.sum(anm_min_increments[anm_min_increments < 0.])
    print float(abs(negatives)) / positives
    
    
#     plt.hist(anm_min_increments, 100, label = "ANM to Min.", alpha = 0.7)
#     plt.hist(anm_increments, 100, label = "Initial to ANM", alpha = 0.7)
#     plt.hist(min_increments, 100, label = "Initial to Minim.", alpha = 0.7)
    # filter between -1 and 1 
    anm_min_increments = anm_min_increments[anm_min_increments >-1.]
    anm_min_increments = anm_min_increments[anm_min_increments <1.]
    sns.kdeplot(anm_min_increments, label = "ANM to Min.", shade=True)
    anm_increments = anm_increments[anm_increments >-1.]
    anm_increments = anm_increments[anm_increments <1.]
    sns.kdeplot(anm_increments, label = "Initial to ANM", shade=True)
    min_increments = min_increments[min_increments >-1.]
    min_increments = min_increments[min_increments <1.]
    sns.kdeplot(min_increments, label = "Initial to Minim.", shade=True)
    
    plt.legend()
    plt.show()
    
    
    
    
    