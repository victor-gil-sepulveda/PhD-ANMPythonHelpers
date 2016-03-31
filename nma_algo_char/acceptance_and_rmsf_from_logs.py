'''
Created on Jan 15, 2016

@author: victor
'''
import numpy
from nma_algo_char.common import MetropolisMCSimulator, create_directory
from anmichelpers.tools.measure import coords_rmsf
from trajectory_comparison.compare_two_rmsfs import rms
import matplotlib.pyplot as plt
from optparse import OptionParser
from pyRMSD.utils.proteinReading import Reader
import os.path
from anmichelpers.tools.tools import norm
import seaborn as sns
from nma_algo_char.data_retrieval import load_cc_data, load_cc_data_from_proc,\
    load_ic_data_from_proc, load_ic_data, process_energy_differences,\
    load_data_from_multiple_processors, load_single_proc_data


# Reference rmsf
def ref_rmsf_calculation(reftraj):
    coords = Reader().readThisFile(reftraj).gettingOnlyCAs().read()
    return coords_rmsf(coords)
    
if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    data_folder = "info"
    reftraj = "/home/victor/Desktop/Desktop/SRCKIN/MD/pro_noh_md.pdb"
#     reftraj = "/media/victor/c2fe358b-c6f7-4562-b2b5-c8d825cc0ed7/MD/Shaw/pro_noh_md.pdb"
#     reftraj = "/home/victor/Desktop/Desktop/UBI/MD/opls_skip10.pdb"
    
    parser = OptionParser()
    parser.add_option("--type", dest="sim_type")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("-t", type = "int", dest="temperature")
    parser.add_option("-n", default = 1, type  = "int", dest="num_procs")
    parser.add_option("-d", default = False, action="store_true", dest="calc_distances")
    
    (options, args) = parser.parse_args()
    
    create_directory(options.results_folder)
    
    if options.num_procs > 1:
        raw_data, min_lens = load_data_from_multiple_processors(options.sim_type, 
                                                                options.num_procs, 
                                                                data_folder,
                                                                10)
    else:
        raw_data, min_len = load_single_proc_data(options.sim_type, data_folder, 10)
    
    
    ####################
    # Acceptance
    ####################
    energy_increments = process_energy_differences(raw_data)
    mc = MetropolisMCSimulator(energy_increments)
    acc = mc.perform_simulation(min(200,len(energy_increments)), 40, options.temperature)
    print "Acceptance ", acc
    open(os.path.join(options.results_folder,"acceptance.txt"),"w").write("%.3f (%.3f)"%acc)

    print "Energy avg.", numpy.mean(energy_increments), numpy.std(energy_increments)
    print "Energy median",numpy.median(energy_increments) 

#     ####################
#     # Rmsf
#     ####################    
#     rmsf_coords = numpy.reshape(raw_data["coords_after"], (len(raw_data["coords_after"]), 
#                                                                            len(raw_data["coords_after"][0])/3, 
#                                                                            3))
#     who_is_accepted = mc.who_is_accepted(options.temperature)
#     rmsf = coords_rmsf(rmsf_coords[who_is_accepted])
#     numpy.savetxt("ic_m1.rmsf", rmsf)
#     rmsf_ref = ref_rmsf_calculation(reftraj)
#     rmsf_ref = rmsf_ref[:-1] # skip last capping CA
#     rms_rmsf = rms(rmsf, rmsf_ref)
#     print "RMS(RMSF)", rms_rmsf
#     
#     plt.plot(rmsf_ref, label = "MD")
#     plt.plot(rmsf, label="MC T = %d"%options.temperature)
#     plt.legend()
#     plt.savefig(os.path.join(options.results_folder,"rmsf.svg"))
#     plt.close()
#     open(os.path.join(options.results_folder,"rms_rmsf.txt"),"w").write("%.3f"%rms_rmsf)
# 
#     if options.calc_distances:
#         ####################
#         # Domain distance
#         ####################    
#         # Position of the atoms to measure distance
#         # CYS:277:CA -> res = 18
#         # LEU:387:CA -> res = 259
#         coords = rmsf_coords[who_is_accepted]
#         distances = []
#         for coordset in coords:
#             distances.append(norm(coordset[128]-coordset[18]))
#         plt.plot(distances)
#         plt.savefig(os.path.join(options.results_folder,"domain_distances.svg"))
#         numpy.savetxt("ic_m1.dists", distances)
#     
    