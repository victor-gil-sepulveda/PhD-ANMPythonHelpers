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
    process_energy_differences, process_after_perturb_rmsd, load_data
import numpy
import matplotlib.pyplot as plt
import seaborn as sns
import cPickle as pickle
from anmichelpers.tools.measure import coords_rmsf
from trajectory_comparison.compare_two_rmsfs import rms

if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("--workspace", dest="workspace")
    parser.add_option("-d", dest="data")
    parser.add_option("--rmsf-reference", dest="rmsf_ref")
    
    
    
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
    results_folder = os.path.join(options.results_folder, os.path.basename(workspace))
    create_directory(results_folder)

    # Things to load
    acceptances = defaultdict(lambda: defaultdict(list))
    def dd_float():
        return defaultdict(float)
    avg_energy = defaultdict(dd_float)
    std_energy = defaultdict(dd_float)
    avg_rmsd = defaultdict(dd_float)
    std_rmsd = defaultdict(dd_float)
    rmsf = defaultdict(dict)
    p1_keys = []
    p2_keys = []
    acceptance_temperatures = [300, 866, 1432, 2000, 2568, 3000]
#     acceptance_temperatures = [3000]
    
    if options.data is None:
        for T in acceptance_temperatures:
            for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
                
                folder_name = "%s_%s_%s_%s_%s"%(experiment_details["prefix"],
                                                experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                                experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
                
                main_folder = "%s_%d"%(experiment_details["prefix"],T)
                data_folder = os.path.join(main_folder, workspace, folder_name,"info")
                
                print data_folder
                
                ## CAUTION: HARDCODED FOLDER ('info'). Must be extracted using the control file
                if experiment_details["prefix"] == "CC":
                    raw_data, min_len = load_data(data_folder, 
                                                  "perturb_energy_before.log",
        #                                           "perturb_energy_after.log", 
                                                "final_energy.log",   # Energy of the whole step
                                                  "initial_cc.log", 
#                                                   "after_anm_cc.log",
                                                  "after_minimization_cc.log", 
                                                  "step_time.log")
                
                if experiment_details["prefix"] == "IC":
                    raw_data, min_len = load_ic_data(os.path.join(workspace, folder_name,"info"))
                
                # skip first frame (usually an outlayer)
                #print "DBG", min_len-1, len(modes), len(process_modes(experiment_details["prefix"], modes, 10))
                energy_increments = process_energy_differences(raw_data)[1:]
                rmsd_increments = process_after_perturb_rmsd(raw_data)[1:]
                
                mc = MetropolisMCSimulator(energy_increments)
                acceptances[T][v1,v2] = mc.perform_simulation(min(200,len(energy_increments)), 20, T)
                avg_energy[T][v1,v2] = numpy.mean(energy_increments)
                std_energy[T][v1,v2] = numpy.std(energy_increments)
                avg_rmsd[T][v1,v2] = numpy.mean(rmsd_increments)
                std_rmsd[T][v1,v2] = numpy.std(rmsd_increments)
                p1_keys.append(v1)
                p2_keys.append(v2)
                
                # Rmsf calculations
                rmsf_coords = numpy.reshape(raw_data["coords_after"], (len(raw_data["coords_after"]), 
                                                                           len(raw_data["coords_after"][0])/3, 
                                                                           3))
                rmsf[T][v1,v2] = coords_rmsf(rmsf_coords)
                
        def default_to_regular(d):
            if isinstance(d, defaultdict):
                d = {k: default_to_regular(v) for k, v in d.iteritems()}
            return d
        pickle.dump((default_to_regular(acceptances), 
                     default_to_regular(avg_energy), 
                     default_to_regular(std_energy), 
                     default_to_regular(avg_rmsd), 
                     default_to_regular(std_rmsd),
                     rmsf, 
                     p1, p1_keys, 
                     p2, p2_keys),
                    open(os.path.join(results_folder,"data"),"w"))
    else:
        (acceptances, 
         avg_energy, std_energy, 
         avg_rmsd, std_rmsd, 
         rmsf,
         p1, p1_keys, 
         p2, p2_keys) = pickle.load(open(options.data))
    
    p1_keys = sorted(list(set(p1_keys)))
    p2_keys = sorted(list(set(p2_keys)))   
    acc_colors = {"high": "yellow","good":"green" ,"low":"red"}
    cat_to_float = {"high": 1.0, "good":0.5 ,"low": 0.0}
    
    good_parameters = defaultdict(list)
    cat_acceptances = defaultdict(dict)
    for T in acceptance_temperatures:
        for key in acceptances[T]:
            acceptance_avg, acceptance_std = acceptances[T][key]
            if acceptance_avg <= 1. and acceptance_avg > 0.35:
                cat_acceptances[T][key] = "high"
            elif acceptance_avg <= 0.35  and acceptance_avg > 0.25:
                cat_acceptances[T][key] = "good"
            else:
                cat_acceptances[T][key] = "low"
        for v1 in p1_keys:
            for v2 in p2_keys:
                if cat_acceptances[T][v1,v2] == "good":
                    good_parameters[T].append((avg_rmsd[T][v1,v2]-std_rmsd[T][v1,v2], 
                                               avg_rmsd[T][v1,v2]+std_rmsd[T][v1,v2],
                                               avg_rmsd[T][v1,v2], 
                                               p1, v1, 
                                               p2, v2))
    
    row_len = 3
    col_len = 2
    f, axes = plt.subplots( col_len, row_len, sharey='row', sharex='col')
    f.subplots_adjust(hspace=0.4, wspace=0.3 )
    f.set_size_inches(10, 6, forward=True)
    matrix = numpy.zeros((len(p1_keys), len(p2_keys)))
    for i,T in enumerate(acceptance_temperatures):
        for j,v1 in enumerate(p1_keys):
            for k,v2 in enumerate(p2_keys):
                matrix[j][k] = acceptances[T][v1,v2][0]#cat_to_float[cat_acceptances[T][v1,v2]]
        ax = axes[i/row_len, i%row_len] 
        sns.heatmap(matrix, #cmap=ListedColormap(['red', 'green', 'yellow']), 
                    ax = ax, square = False, cbar = True, center = 0.5,
                    annot = True)
        ax.set_xticklabels([str(i) for i in p2_keys])
        ax.set_yticklabels([str(i) for i in sorted(p1_keys, reverse=True)])
        ax.set_title("T = %d"%T)
#     plt.show()
    plt.savefig(os.path.join(results_folder,"acceptance_avg.svg"))
    plt.close()
    
    # Do rmsf comparisons
    if options.rmsf_ref:
        f, axes = plt.subplots( col_len, row_len, sharey='row', sharex='col')
        f.subplots_adjust(hspace=0.4, wspace=0.3 )
        f.set_size_inches(10, 6, forward=True)
        rmsf_matrix = numpy.zeros((len(p1_keys), len(p2_keys)))
        rmsf_ref = numpy.loadtxt(options.rmsf_ref)[:-1]
        for i,T in enumerate(acceptance_temperatures):
            for j,v1 in enumerate(p1_keys):
                for k,v2 in enumerate(p2_keys):
                    rmsf_matrix[j][k] = rms(rmsf_ref,rmsf[T][v1,v2])
            ax = axes[i/row_len, i%row_len] 
            sns.heatmap(rmsf_matrix, #cmap=ListedColormap(['red', 'green', 'yellow']), 
                    ax = ax, square = False, cbar = True, #center = 0.5,
                    annot = True)
        plt.savefig(os.path.join(results_folder,"rmsf_avg.svg"))
        plt.close()
        
        sns.set_style("whitegrid", {"lines.linewidth": ".7"})
        f, axes = plt.subplots( col_len, row_len, sharey='row', sharex='col')
        f.subplots_adjust(hspace=0.2, wspace=0.3 )
        f.set_size_inches(10, 6, forward=True)
        for i,T in enumerate(acceptance_temperatures):
            ax = axes[i/row_len, i%row_len]
            ax.plot(rmsf_ref, lw=0.8)
            ax.set_title("T = %d"%T)
            if i%row_len == 0:
                ax.set_ylabel("RMSD($\AA$)")
            if i/row_len == 1:
                ax.set_xlabel("Num. Residue")
                
            for v1 in p1_keys:
                for v2 in p2_keys:
                    ax.plot(rmsf[T][v1,v2], lw=0.5)
        plt.savefig(os.path.join(results_folder,"rmsf.svg"))
        plt.close()
        
    # for each displacement pick smaller temperature
    good_parameters = {}
    for T in acceptance_temperatures:
        good_parameters[T] = []
        for j,v1 in enumerate(p1_keys):
            for k,v2 in enumerate(p2_keys):
                    # range 0.25 0.35 
                    if avg_rmsd[T][v1,v2]+std_rmsd[T][v1,v2] >= 0.25 and avg_rmsd[T][v1,v2]-std_rmsd[T][v1,v2] <= 0.35:
                        good_parameters[T].append((avg_rmsd[T][v1,v2], v1, v2, rms(rmsf_ref,rmsf[T][v1,v2])))
    
    print "T\tAcceptance\tdisp\trmsg\trms(rmsf)"
    for T in good_parameters:
        best_for_T = sorted(good_parameters[T], reverse=True)[0]
        print T,"\t",
        values = ["%.3f"%f for f in  best_for_T]
        for val in values:print val,"\t",
        print