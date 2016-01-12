'''
Created on 14/12/2015

@author: victor
'''
from optparse import OptionParser
from nma_algo_char.common import load_control_json, create_directory,\
    pair_parameter_values, parameter_value_to_string, MetropolisMCSimulator
import os
from _collections import defaultdict
from nma_algo_char.mode_application_analysis import load_cc_data, load_ic_data,\
    process_energy_differences, process_after_perturb_rmsd, load_data
import numpy
import matplotlib.pyplot as plt
import seaborn as sns
import cPickle as pickle
from anmichelpers.tools.measure import coords_rmsf
from trajectory_comparison.compare_two_rmsfs import rms
from operator import itemgetter

if __name__ == '__main__':
    sns.set_style("whitegrid")

    
    parser = OptionParser()
    parser.add_option("--results", dest="results_folder")
    parser.add_option("-d", dest="data")
    parser.add_option("--rmsf-reference", dest="rmsf_ref")
    
    (options, args) = parser.parse_args()
    
    if not options.results_folder:  
        parser.error('Results folder not given')
    
    create_directory(options.results_folder)
    
    IC_params = {
                    (300,0.5,0.05),
                    (583,0.8,0.1),
                    (866,0.8,0.1),
                    (1150,1.1,0.05),
                    (1432,1.1,0.1),
                    (2000,1.1,0.05),
                    (2568,1.4,0.1)
                 }
    
    CC_params = {
                    (300,0.66,0.1),
                    (583,1.5,0.01),
                    (866,1.92,0.05),
                    (1150,2.34,0.01),
                    (1432,4.02,0.1),
                    (2000,3.17,0.01),
                    (2568,4.02,0.01)
                 }
    
    datasets = [{
                 "prefix":"IC",
                 "parameters": IC_params,
                 "source_folder":"IC_1500",
                 "workspace":"ic_open_rmsg_dispf"
                 },
                {
                 "prefix":"CC",
                 "parameters": CC_params,
                 "source_folder":"CC_1500",
                 "workspace":"cc_open_rmsg_dispf"
                 }
                ]
    
    
    step_granularity = 100
    max_samples = 1000
    
    if options.data is None:
        acceptances = {}
        rmsf = {}
        rmsd = {}
        rmsf_ref = numpy.loadtxt(options.rmsf_ref)[:-1]
        for dataset in datasets:
            acceptances[dataset["prefix"]] = {}
            rmsf[dataset["prefix"]] = {}
            rmsd[dataset["prefix"]] = {}
            for T,v1,v2 in dataset["parameters"]:
                folder_name = "%s_%s_%s_%s_%s"%(dataset["prefix"], "dispFact", parameter_value_to_string(v1),"rmsg", parameter_value_to_string(v2))
                main_folder = "%s_%d"%(dataset["prefix"],T)
                data_folder = os.path.join(dataset["source_folder"],
                                           main_folder, 
                                           dataset["workspace"], 
                                           folder_name,
                                           "info")
                
                print data_folder
                if dataset["prefix"] == "CC":

                    raw_data, min_len = load_data(data_folder, 
                                                  "perturb_energy_before.log",
                                                  "final_energy.log",   # Energy of the whole step
                                                  "initial_cc.log", 
                                                  "after_minimization_cc.log", 
                                                  "step_time.log",
                                                  max_samples = max_samples)
                
                if dataset["prefix"] == "IC":
                    raw_data, min_len = load_ic_data(data_folder, max_samples)
                
                # skip first frame (usually an outlayer)
                energy_increments = process_energy_differences(raw_data)[1:]
                rmsd_increments = process_after_perturb_rmsd(raw_data)[1:]
                # Coordinates reshape for rmsf
                rmsf_coords = numpy.reshape(raw_data["coords_after"], (len(raw_data["coords_after"]), 
                                                                               len(raw_data["coords_after"][0])/3, 
                                                                               3))
                
                acceptances[dataset["prefix"]][T,v1,v2] = []
                rmsd[dataset["prefix"]][T,v1,v2] = []
                rmsf[dataset["prefix"]][T,v1,v2] = []
                for num_steps in range(step_granularity, max_samples+1, step_granularity):
                    # Acceptance calc.
                    mc = MetropolisMCSimulator(energy_increments[0:num_steps])
                    acceptances[dataset["prefix"]][T,v1,v2].append(mc.perform_simulation(min(num_steps/2, len(energy_increments)),40,T))
                    # Rmsf calculations
                    _rmsf = coords_rmsf(rmsf_coords[0:num_steps])
                    rmsf[dataset["prefix"]][T,v1,v2].append(rms(rmsf_ref,_rmsf))
                    # RMSD avg. calc
                    rmsd[dataset["prefix"]][T,v1,v2].append((numpy.mean(rmsd_increments[0:num_steps]), 
                                                             numpy.std(rmsd_increments[0:num_steps])))
                
            pickle.dump((acceptances,rmsf,rmsd),
                       open(os.path.join(options.results_folder,"per_step_data"),"w"))
    else:
        acceptances,rmsf,rmsd = pickle.load(open(options.data))
    
    #-----------
    # PLOTS
    #-----------
    all_temperatures = [300, 583, 866, 1150, 1432, 2000, 2568, 300, 300]
    
    x = range(step_granularity, max_samples+1, step_granularity)
    print x
    for T in all_temperatures:
        for prefix in ["IC","CC"]:
            for key in acceptances[prefix]:
                if key[0] == T: 
                    print T, key
                    y = rmsd[prefix][key]
                    print x, y
                    plt.plot(x,y,label = prefix)
        plt.legend()
        plt.show()
    