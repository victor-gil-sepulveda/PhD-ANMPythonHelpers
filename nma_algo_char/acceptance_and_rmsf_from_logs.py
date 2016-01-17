'''
Created on Jan 15, 2016

@author: victor
'''
from nma_algo_char.mode_application_analysis import load_data, load_ic_data,\
    process_energy_differences
import numpy
from nma_algo_char.common import MetropolisMCSimulator
from anmichelpers.tools.measure import coords_rmsf
from trajectory_comparison.compare_two_rmsfs import rms
import matplotlib.pyplot as plt
from optparse import OptionParser
from pyRMSD.utils.proteinReading import Reader

def load_data_from_multiple_processors(sim_type, num_procs, data_folder, max_samples = numpy.inf):
    # processors from 1 to num_procs-1 have data
    all_raw_data = {}
    for proc in range(1,num_procs):
        if sim_type == "CC":
            raw_data, _ = load_data(data_folder, 
                                 "perturb_energy_before.p%d.log"%proc, 
                                 "final_energy.p%d.log"%proc,  
                                 "initial_cc.p%d.log"%proc, 
                                 "after_minimization_cc.p%d.log"%proc, 
                                 "step_time.p%d.log"%proc,
                                 max_samples)
        if sim_type == "IC":
            raw_data, _ = load_data(data_folder, 
                             "ener_mc_move_before.p%d.log"%proc, 
                             "ener_mc_move_after.p%d.log"%proc,  
                             "ca_mc_move_before.p%d.log"%proc, 
                             "ca_mc_move_after.p%d.log"%proc, 
                             "step_time.p%d.log"%proc,
                             max_samples)
        for key in raw_data:
            if key in all_raw_data:
                all_raw_data[key].extend(raw_data[key])
            else:
                all_raw_data[key] = list(raw_data[key])
    
    for key in all_raw_data:
        all_raw_data[key] = numpy.array(all_raw_data[key])
        print len(all_raw_data[key])
    
    return all_raw_data


# Reference rmsf
def ref_rmsf_calculation(reftraj):
    coords = Reader().readThisFile(reftraj).gettingOnlyCAs().read()
    return coords_rmsf(coords)
    
if __name__ == '__main__':
    data_folder = "info"
    rmsf_reference = "/home/victor/Desktop/NMA_acceptance_T/pro_noh_md.pdb.rmsf"
    reftraj = "/home/victor/Desktop/NMA_acceptance_T/pro_noh_md.pdb"
    
    parser = OptionParser()
    parser.add_option("--type", dest="sim_type")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("-t", type = "int", dest="temperature")
    parser.add_option("-n", type  = "int", dest="num_procs")
    
    (options, args) = parser.parse_args()
    
    raw_data = load_data_from_multiple_processors(options.sim_type, options.num_procs, data_folder)
        
    energy_increments = process_energy_differences(raw_data)[1:]
    mc = MetropolisMCSimulator(energy_increments)
    acc = mc.perform_simulation(min(200,len(energy_increments)), 40, options.temperature)
    print "Acceptance ", acc
    
    rmsf_coords = numpy.reshape(raw_data["coords_after"], (len(raw_data["coords_after"]), 
                                                                           len(raw_data["coords_after"][0])/3, 
                                                                           3))
    
    rmsf = coords_rmsf(rmsf_coords[mc.who_is_accepted(options.temperature)])
#     rmsf_ref = numpy.loadtxt(rmsf_reference)[:-1]
    rmsf_ref = ref_rmsf_calculation(reftraj)
    rmsf_ref = rmsf_ref[:-1] # skip last capping CA
    print "RMS(RMSF)", rms(rmsf, rmsf_ref)
    
    plt.plot(rmsf_ref)
    plt.plot(rmsf)
    plt.show()
    
    numpy.savetxt("rmsf_ref", rmsf_ref)
    numpy.savetxt("rmsf", rmsf)
    