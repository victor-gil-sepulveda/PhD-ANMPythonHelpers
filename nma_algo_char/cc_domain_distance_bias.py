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
from nma_algo_char.data_retrieval import add_coords_data

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
        raw_data, min_len = load_data_from_multiple_processors("CC", 
                                                               options.num_procs, 
                                                               data_folder,
                                                               skip_first = 0,
                                                               #max_samples=500,
                                                               full_pele_energy=True)
    else:
        raw_data, min_len = load_single_proc_data("CC", 
                                                  data_folder,
                                                  skip_first = 1,
                                                  #max_samples=500,
                                                  full_pele_energy=True)
        
    raw_data = add_coords_data(raw_data, data_folder, "after_minimization_cc", options.num_procs, min_len)
    raw_data = add_coords_data(raw_data, data_folder, "proposal_cc", options.num_procs, min_len)
    
    energy_increments = process_energy_differences(raw_data)
    mc = MetropolisMCSimulator(energy_increments)
    who_is_accepted = mc.who_is_accepted(options.temperature)
#     who_is_accepted = range(len(energy_increments))
    
    initial_coords = numpy.reshape(raw_data["coords_before"], (len(raw_data["coords_before"]), 
                                                               len(raw_data["coords_before"][0])/3, 
                                                               3))
    
    proposal_coords = numpy.reshape(raw_data["proposal_cc"], (len(raw_data["proposal_cc"]), 
                                                            len(raw_data["proposal_cc"][0])/3, 
                                                            3))
    
    minim_coords = numpy.reshape(raw_data["after_minimization_cc"], (len(raw_data["after_minimization_cc"]), 
                                                           len(raw_data["after_minimization_cc"][0])/3, 
                                                           3))
    
    anm_coords = numpy.reshape(raw_data["coords_after"], (len(raw_data["coords_after"]), 
                                                           len(raw_data["coords_after"][0])/3, 
                                                           3))
    def calc_distances(in_coords, who_is_accepted):
        coords = in_coords[who_is_accepted]
        distances = []
        for coordset in coords:
            distances.append(norm(coordset[128]-coordset[18]))
        return numpy.array(distances)
    
    initial_distances = calc_distances(initial_coords, who_is_accepted)
    proposal_distances = calc_distances(proposal_coords, who_is_accepted)
    anm_distances = calc_distances(anm_coords, who_is_accepted) 
    min_distances = calc_distances(minim_coords, who_is_accepted) 
    
    proposal_increments = proposal_distances - initial_distances 
    proposal_anm_increments = anm_distances - proposal_distances
    anm_increments = anm_distances - initial_distances
    anm_min_increments = min_distances - anm_distances
    min_increments = min_distances - initial_distances
    
    print "Initial to Proposal negative to positive count ratio:  "
    positives = numpy.count_nonzero(proposal_increments > 0.)
    negatives = numpy.count_nonzero(proposal_increments < 0.)
    print "+", positives, "-", negatives,  "ratio:", float(negatives) / positives
    
    print "Proposal to ANM negative to positive count ratio:  "
    positives = numpy.count_nonzero(proposal_anm_increments > 0.)
    negatives = numpy.count_nonzero(proposal_anm_increments < 0.)
    print "+", positives, "-", negatives, "ratio:", float(negatives) / positives
    
    print "Initial to Anm. negative to positive count ratio:  "
    positives = numpy.count_nonzero(anm_increments > 0.)
    negatives = numpy.count_nonzero(anm_increments < 0.)
    print  "+", positives, "-", negatives, "ratio:", float(negatives) / positives
    
    print "Initial to min. negative to positive count ratio:  "
    positives = numpy.count_nonzero(min_increments > 0.)
    negatives = numpy.count_nonzero(min_increments < 0.)
    print  "+", positives, "-", negatives, "ratio:", float(negatives) / positives
    
    print "Anm to min. negative to positive count ratio:  "
    positives = numpy.count_nonzero(anm_min_increments > 0.)
    negatives = numpy.count_nonzero(anm_min_increments < 0.)
    print  "+", positives, "-", negatives, "ratio:", float(negatives) / positives
    
    print "*****************"
    print "Initial to Proposal negative to positive sum ratio:\t",
    positives = numpy.sum(proposal_increments[proposal_increments > 0.])
    negatives = numpy.sum(proposal_increments[proposal_increments < 0.])
    print float(abs(negatives)) / positives
    
    print "Proposal to Anm negative to positive sum ratio:\t",
    positives = numpy.sum(proposal_anm_increments[proposal_anm_increments > 0.])
    negatives = numpy.sum(proposal_anm_increments[proposal_anm_increments < 0.])
    print float(abs(negatives)) / positives
    
    print "Initial to Anm negative to positive sum ratio:\t",
    positives = numpy.sum(anm_increments[anm_increments > 0.])
    negatives = numpy.sum(anm_increments[anm_increments < 0.])
    print float(abs(negatives)) / positives
    
    print "Initial to min. negative to positive sum ratio:\t",
    positives = numpy.sum(min_increments[min_increments > 0.])
    negatives = numpy.sum(min_increments[min_increments < 0.])
    print float(abs(negatives)) / positives
    
    print "Anm to min. negative to positive sum ratio:\t",
    positives = numpy.sum(anm_min_increments[anm_min_increments > 0.])
    negatives = numpy.sum(anm_min_increments[anm_min_increments < 0.])
    print float(abs(negatives)) / positives
    
    
    # filter between -1 and 1 
    
    proposal_increments = proposal_increments[proposal_increments >-1.]
    proposal_increments = proposal_increments[proposal_increments <1.]
    sns.kdeplot(proposal_increments, label = "Intial to proposal", shade=True)
    
    proposal_anm_increments = proposal_anm_increments[proposal_anm_increments >-1.]
    proposal_anm_increments = proposal_anm_increments[proposal_anm_increments <1.]
    sns.kdeplot(proposal_anm_increments, label = "Proposal to ANM", shade=True)
    
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
    
    
    
    
    