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
    process_energy_differences, process_after_perturb_rmsd
import numpy
import matplotlib.pyplot as plt


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("--workspace", dest="workspace")
    
    (options, args) = parser.parse_args()
    
    if not options.experiment:  
        parser.error('Experiment file not given')
    
    if not options.results_folder:  
        parser.error('Results folder not given')
    
    experiment_details = load_control_json(options.experiment)
    if options.workspace is not None:
        workspace = os.path.normpath(options.workspace)
    else:
        workspace = os.path.normpath(experiment_details["workspace"])
    
    # create a place for the results
    create_directory(os.path.join(options.results_folder, os.path.basename(workspace)))

    # Things to load
    acceptances = defaultdict(lambda: defaultdict(list))
    avg_energy = defaultdict(list)
    std_energy = defaultdict(list)
    avg_rmsd = defaultdict(list)
    std_rmsd = defaultdict(list)
    p1_keys = []
    p2_keys = []
    acceptance_temperatures = [300, 500, 700, 900, 1000, 1200, 1500, 1700, 1900, 2100, 2300, 2500, 2700, 2900]
    
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        folder_name = "%s_%s_%s_%s_%s"%(experiment_details["prefix"],
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
        
        
        ## CAUTION: HARDCODED FOLDER ('info'). Must be extracted using the control file
        if experiment_details["prefix"] == "CC":
            raw_data, min_len = load_cc_data(os.path.join(workspace, folder_name,"info"))
        
        if experiment_details["prefix"] == "IC":
            raw_data, min_len = load_ic_data(os.path.join(workspace, folder_name,"info"))
        
        # skip first frame (usually an outlayer)
        #print "DBG", min_len-1, len(modes), len(process_modes(experiment_details["prefix"], modes, 10))
        energy_increments = process_energy_differences(raw_data)[1:]
        rmsd_increments = process_after_perturb_rmsd(raw_data)[1:]
        
        mc = MetropolisMCSimulator(energy_increments)
        for T in acceptance_temperatures :
            acceptances[T][v1,v2] = mc.perform_simulation(min(100,len(energy_increments)), 20, T)
        avg_energy[v1,v2] = numpy.mean(energy_increments)
        std_energy[v1,v2] = numpy.std(energy_increments)
        avg_rmsd[v1,v2] = numpy.mean(rmsd_increments)
        std_rmsd[v1,v2] = numpy.std(rmsd_increments)
        p1_keys.append(v1)
        p2_keys.append(v2)
    
    p1_keys = sorted(list(set(p1_keys)))
    p2_keys = sorted(list(set(p2_keys)))   
    acc_colors = {"high": "yellow","good":"green" ,"low":"red"}
    cat_to_float = {"high": 1.0, "good":0.3 ,"low": 0}
    for T in acceptance_temperatures:
        cat_acceptances = {}
        for key in acceptances[T]:
            acceptance_avg, acceptance_std = acceptances[T][key]
            if acceptance_avg <= 1. and acceptance_avg > 0.35:
                cat_acceptances[key] = "high"
            elif acceptance_avg <= 0.35  and acceptance_avg > 0.25:
                cat_acceptances[key] = "good"
            else:
                cat_acceptances[key] = "low"
        hue = []
        x = []
        y = []
        for v1 in p1_keys:
            for v2 in p2_keys:
                x.append(v1)
                y.append(v2)
                hue.append(cat_acceptances[v1,v2])
        vals_by_hue = get_values_by_hue(x,y,hue)
        for hue_val in vals_by_hue:
            plt.scatter(vals_by_hue[hue_val]["x"], 
                        vals_by_hue[hue_val]["y"], 
                        label = hue_val, 
                        color = acc_colors[hue_val])
        
        plt.show()
        
        matrix = numpy.zeros((len(p1_keys), len(p2_keys)))
        for i,v1 in enumerate(p1_keys):
            for j,v2 in enumerate(p2_keys):
                matrix[i][j] = cat_to_float[cat_acceptances[v1,v2]]
        plt.matshow(matrix)
        plt.show()