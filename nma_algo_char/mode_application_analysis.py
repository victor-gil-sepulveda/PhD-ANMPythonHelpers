"""
Created on Nov 19, 2015

@author: victor
"""
import numpy
import os
from optparse import OptionParser
from nma_algo_char.common import load_control_json, pair_parameter_values,\
    parameter_value_to_string, create_directory, MetropolisMCSimulator,\
    scatter_plot_by_hue
from collections import defaultdict
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import cPickle as pickle
from nma_algo_char.data_retrieval import load_ic_data, load_cc_data,\
    process_after_perturb_rmsd, process_energy_differences, get_mode_frequencies
from imghdr import what

def remove_energy_outlayers(all_data, ENERGY_LABEL):
    outlayer_margin  = int(len(all_data[ENERGY_LABEL])*0.98)
    ener_low = min(sorted(all_data[ENERGY_LABEL], reverse = True)[:outlayer_margin])
    ener_high = max(sorted(all_data[ENERGY_LABEL])[:outlayer_margin])
    

    indices = numpy.logical_and(all_data[ENERGY_LABEL] < ener_high,
                                all_data[ENERGY_LABEL] > ener_low)
                                 
    for label in all_data:
        print label , "initial len ", len(all_data[label]),
        all_data[label] = numpy.array(all_data[label])
        all_data[label] = all_data[label][indices]
        print "final len ", len(all_data[label])
    
    return ener_low, ener_high

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("--workspace", dest="workspace")
    parser.add_option("--folder", default= "info" ,dest="folder")
    parser.add_option("--plot", action = "store_true", default= False, dest="do_plots")
    
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
    
    # Create a place for the results
    create_directory(os.path.join(options.results_folder, os.path.basename(workspace)))
    
    # Initialize stuff
    all_data = defaultdict(list)
    relax_iterations = []
    acceptances = defaultdict(lambda: defaultdict(list))
    avg_energy = defaultdict(list)
    std_energy = defaultdict(list)
    avg_rmsd = defaultdict(list)
    avg_time = defaultdict(list)
    norm_rmsd = defaultdict(list)
    norm_energy = defaultdict(list)
    mode_frequencies = defaultdict(list)
    modes_p_v = defaultdict(list)
    
    ENERGY_LABEL = "$\Delta$ U"
    RMSD_LABEL = "RMSD"
    nmd_file_name = {"CC":"normalized_modes.1.nmd", 
                     "IC":"normalized_modes_cc.1.nmd"}
    
    v1s = []
    v2s = []
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        v1s.append(v1)
        v2s.append(v2)
        
        folder_name = "%s_%s_%s_%s_%s"%(experiment_details["prefix"],
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
        
        if experiment_details["prefix"] == "CC":
            raw_data, data_len = load_cc_data(os.path.join(workspace, 
                                                          folder_name,
                                                          options.folder),
                                             full_pele_energy = False,
                                             skip_first = 15)
        
        if experiment_details["prefix"] == "IC":
            raw_data, data_len = load_ic_data(os.path.join(workspace, 
                                                          folder_name,
                                                          options.folder),
                                             skip_first = 15)

        # Start processing        
        modes = raw_data["modes"]
        _mode_frequencies = get_mode_frequencies(modes, 
                                                os.path.join(workspace, folder_name, "info",
                                                    nmd_file_name[experiment_details["prefix"]])
                                                )
        
        energy_increments = process_energy_differences(raw_data)
        rmsd_increments = process_after_perturb_rmsd(raw_data)
        acc_mean_and_avg = MetropolisMCSimulator(energy_increments).perform_simulation(
                                                    min(100,len(energy_increments)), 20, 300)
        
        # Fill all data structure
        all_data[ENERGY_LABEL].extend(energy_increments)
        all_data[RMSD_LABEL].extend(rmsd_increments)
        all_data["Mode"].extend(modes)
        all_data["time_per_step"].extend(raw_data["time_per_step"])
        all_data[p1].extend([v1]*data_len)
        all_data[p2].extend([v2]*data_len)
        
        # Fill the other structures
        acceptances[v1,v2] = acc_mean_and_avg
        avg_rmsd[v1,v2] = (numpy.mean(rmsd_increments),
                           numpy.std(rmsd_increments))
        avg_energy[v1,v2] = numpy.mean(energy_increments)
        std_energy[v1,v2] = numpy.std(energy_increments)
        print std_energy[v1,v2]
        avg_time[v1,v2] = (numpy.mean(raw_data["time_per_step"]),
                           numpy.std(raw_data["time_per_step"]))
        norm_rmsd[v1,v2] = numpy.array(rmsd_increments) / numpy.max(rmsd_increments)
        norm_energy[v1,v2] = numpy.array(energy_increments) / numpy.max(numpy.abs(energy_increments))
        modes_p_v[v1,v2] = numpy.array(modes)+1 # Modes will start from index 1 in the plots
        mode_frequencies[v1,v2] = _mode_frequencies
    
    v1s = sorted(list(set(v1s)))
    v2s = sorted(list(set(v2s)))
    
    # Remove outliers
    energy_with_outlayers = numpy.array(all_data[ENERGY_LABEL])
    ener_low, ener_high = remove_energy_outlayers(all_data, ENERGY_LABEL)
    
    # Save all_data in pickled format
    pickle.dump(all_data, open(os.path.join(options.results_folder,os.path.basename(workspace),"all_data.pickle"),"w"))
    
    db = DataFrame.from_dict(all_data, orient="index")
    db.transpose().to_csv(os.path.join(options.results_folder,os.path.basename(workspace),"data.csv"))
    db.columns = db.columns.get_level_values(0)
    
    if options.do_plots:
        import seaborn as sns
        sns.set_style("whitegrid")
        
        # Facet grid to see rmsd vs energy vs dispfactor vs relaxation whatever
        g = sns.FacetGrid(db.transpose(), col=p1, row=p2, hue="Mode", 
                          sharex = True, sharey = True, margin_titles=True,
                          legend_out = True, ylim= (ener_low, ener_high))
        g.map(plt.scatter, RMSD_LABEL, ENERGY_LABEL)
        g.add_legend(label_order = sorted(g._legend_data.keys()))
        g.fig.suptitle("RMSD - $\Delta U$ - displacement - relaxation intensity (colored by mode)")
        g.savefig(os.path.join(options.results_folder,os.path.basename(workspace),os.path.basename(workspace)+"_u_rmsd.svg"))
        plt.close()
 
        # Global rmsd vs time, color by dispfact and relax
        colors = sns.color_palette("hls", len( set(all_data[p1])))
        ax = plt.subplot2grid((2,2), (0,0))
        scatter_plot_by_hue(all_data[RMSD_LABEL], all_data["time_per_step"], all_data[p1], colors)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
        ax.set_ylabel("Tps (s) (hue=%s)"%p1)
         
        ax = plt.subplot2grid((2,2), (1,0))
        scatter_plot_by_hue(all_data[RMSD_LABEL], all_data["time_per_step"], all_data[p2], colors)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
        ax.set_ylabel("Tps (s) (hue=%s)"%p2)
        ax.set_xlabel("RMSD")
 
        # Global energy vs time, color by dispfact and relax
        ax = plt.subplot2grid((2,2), (0,1))
        scatter_plot_by_hue(all_data[ENERGY_LABEL], all_data["time_per_step"], all_data[p1], colors)
         
        ax = plt.subplot2grid((2,2), (1,1))
        scatter_plot_by_hue(all_data[ENERGY_LABEL], all_data["time_per_step"], all_data[p2], colors)
        ax.set_xlabel("Energy" )
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace),os.path.basename(workspace)+"_vs_timespep.svg"))
        plt.close()
         
        # Define categories for acceptance
        # from 100% to 35%, yellow, high
        # from 35% to 25%, green, good acceptance
        # from 0 to 25%, red, bad acceptance
        cat_acceptances = {}
        for key in acceptances:
            acceptance = acceptances[key][0]
            if acceptance <= 1. and acceptance > 0.35:
                cat_acceptances[key] = "high"
            elif acceptance <= 0.35  and acceptance > 0.25:
                cat_acceptances[key] = "good"
            else:
                cat_acceptances[key] = "low"
         
        # rmsd vs time and color by acceptance
        plt.figure()
        keys = sorted(avg_rmsd.keys())
        labels = []
        colors = {"high": "yellow","good":"green" ,"low":"red"}
         
        for  acc_cat in ["high","good","low"]:
            x = []
            y = []
            for x_k in avg_rmsd:
                if cat_acceptances[x_k] == acc_cat:
                    x.append(avg_rmsd[x_k][0])
                    y.append(avg_time[x_k][0])
            plt.scatter(x, y, label = acc_cat, color = colors[acc_cat])
            #plt.errorbar(x, y, xerr=avg_rmsd[key][1], yerr=acceptances[key][1],  fmt='o')
        plt.title("Avg. RMSD vs Avg. Time per Step")
        plt.xlabel("Avg. RMSD (${\AA}$)")
        plt.ylabel("Avg. Time per Step (s)")    
        plt.legend()
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace),os.path.basename(workspace)+"_rmsd_vs_tps_hue_acc.svg"))
        plt.close()
         
        # inc_U vs time and color by acceptance
        plt.figure()
        keys = sorted(avg_rmsd.keys())
        labels = []
        colors = {"high": "yellow","good":"green" ,"low":"red"}
        for  acc_cat in ["high","good","low"]:
            x = []
            y = []
            for x_k in avg_rmsd:
                if cat_acceptances[x_k] == acc_cat:
                    x.append(avg_energy[x_k])
                    y.append(avg_time[x_k][0])
            plt.scatter(x, y, label = acc_cat, color = colors[acc_cat])
            #plt.errorbar(x, y, xerr=avg_rmsd[key][1], yerr=acceptances[key][1],  fmt='o')
        plt.title("Avg. $\Delta U$ vs Time per Step (hue=acceptance)")
        plt.xlabel("Avg. $\Delta U$ (kcal/mol) ")
        plt.ylabel("Avg. Time per Step (s)")    
        plt.legend()
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace),os.path.basename(workspace)+"_energy_vs_tps_hue_acc.svg"))
        plt.close()
         
        # Do rmsd/energy vs acceptance
        plt.figure()
        keys = sorted(avg_rmsd.keys())
        labels = []
        colors = sns.color_palette("hls", len(keys))
        for i, key in enumerate(keys):
            x = avg_rmsd[key][0] / avg_energy[key]
            y = acceptances[key][0]
            label = "%.2f %.2f"%key
            plt.scatter(x, y, label = label, color = colors[i])
            #plt.errorbar(x, y, xerr=avg_rmsd[key][1], yerr=acceptances[key][1],  fmt='o')
            plt.annotate(
                label, 
                xy = (x, y), xytext = (5, 5),
                textcoords = 'offset points', ha = 'right', va = 'bottom', size=6)
        plt.xlabel("RMSD / energy")
        plt.ylabel("Acceptance")    
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace),os.path.basename(workspace)+"_rmsd_energy_vs_acc.svg"))
        plt.close()
         
        # Do rmsd vs acceptance
        plt.figure()
        keys = sorted(avg_rmsd.keys())
        labels = []
        colors = sns.color_palette("hls", len(keys))
        for i, key in enumerate(keys):
            x = avg_rmsd[key][0]
            y = acceptances[key][0]
            label = "%.2f %.2f"%key
            plt.scatter(x, y, label = label, color = colors[i])
            plt.errorbar(x, y, xerr = avg_rmsd[key][1], yerr=acceptances[key][1],  fmt='o')
            plt.annotate(
                label, 
                xy = (x, y), xytext = (5, 5),
                textcoords = 'offset points', ha = 'right', va = 'bottom', size=6)
        plt.xlabel("RMSD")
        plt.ylabel("Acceptance")   
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace),os.path.basename(workspace)+"_rmsd_vs_acc.svg"))
        plt.close()
         
        # Do energy vs acceptance
        plt.figure()
        keys = sorted(avg_rmsd.keys())
        labels = []
        colors = sns.color_palette("hls", len(keys))
        for i, key in enumerate(keys):
            x = avg_energy[key]
            y = acceptances[key][0]
            label = "%.2f %.2f"%key
            plt.errorbar(x, y, xerr = std_energy[key], yerr=acceptances[key][1],  color = colors[i])
            plt.scatter(x, y, label = label, color = colors[i])
            plt.annotate(
                label, 
                xy = (x, y), xytext = (5, 5),
                textcoords = 'offset points', ha = 'right', va = 'bottom', size=6)
        plt.xlabel("Energy")
        plt.ylabel("Acceptance")   
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace),os.path.basename(workspace)+"_energy_vs_acc.svg"))
        plt.close()
        
        # Do energy vs rmsd
        plt.figure()
        keys = sorted(avg_rmsd.keys())
        labels = []
        colors = sns.color_palette("hls", len(keys))
        markers = ["o", "s", "h", "*", "+", "D"]
        colors = sns.color_palette("hls", len(v1s))
        for i, v1 in enumerate(v1s):
            for j,v2 in enumerate(v2s):
                key = (v1,v2)
                x = avg_rmsd[key][0]
                x_err =  avg_rmsd[key][1]
                y = avg_energy[key]
                y_err = std_energy[key]
                label = "%.2f %.2f"%key
                plt.errorbar(x, y, 
                             xerr=x_err, yerr=y_err, 
                             marker=markers[j], 
                             color=colors[i],
                             label = label, mec='black',mew=1)
                plt.scatter(x, y, marker=markers[j], linewidths=1,  edgecolors='black', color=colors[i])
        plt.xlabel("RMSD")
        plt.ylabel("$\Delta U$")    
        plt.legend()
        plt.savefig(os.path.join(options.results_folder,
                                 os.path.basename(workspace),
                                 os.path.basename(workspace)+"_rmsd_vs_energy_avgs.svg"))
        plt.close()
         

        other_results = open(os.path.join(options.results_folder,
                                 os.path.basename(workspace),"other_results.txt"),"w")
        
        def save_result(what):
            print what
            other_results.write(what+"\n")

        #---------------------------------------        
        # Energy vs RMSD correlation
        #---------------------------------------
        ranked_rmsd = scipy.stats.rankdata(numpy.array(all_data[RMSD_LABEL]))
        ranked_energy = scipy.stats.rankdata(numpy.array(all_data[ENERGY_LABEL]))
        min_rank = min(numpy.max(ranked_rmsd),numpy.max(ranked_energy))
        if min_rank == numpy.min(ranked_rmsd):
            ranked_energy *= min_rank /numpy.max(ranked_energy)
        else:
            ranked_rmsd *= min_rank /numpy.max(ranked_rmsd)
        rho, p_val =  scipy.stats.spearmanr(ranked_rmsd, ranked_energy)
        save_result( "Does the RMSD and energy inc. correlate?\t%.2f\t%.2f"%(rho,p_val))
        
        #---------------------------------------        
        # Displacement (normalized) vs Frequency
        #---------------------------------------
        all_freqs = []
        for key in mode_frequencies: all_freqs.extend(mode_frequencies[key])
        all_norm_rmsd = []
        for key in norm_rmsd: all_norm_rmsd.extend(norm_rmsd[key])
        all_modes = []
        for key in modes_p_v: all_modes.extend(modes_p_v[key])
        
        # X is the independent var. In this case a discrete (or even categorical) one.
        # Y is the dependent var, is continuous. 
#         sns.regplot(numpy.array(all_freqs), numpy.array(all_norm_rmsd))
#         plt.show()
        
        # We can use spearman rank (rho, symmetric) to check the "association force" or correlation
        # For this it is needed to convert both to categorical ranked variables.
        # frequencies is already a ranked categorical variable (is it discretized? It may be not
        # because we do not have, and it is not possible to have, values in the full range)
        ranked_norm_rmsd =  scipy.stats.rankdata(all_norm_rmsd)
        # rescale in 1-10 like modes
        scaled_ranked_norm_rmsd = ranked_norm_rmsd*10./numpy.max(ranked_norm_rmsd)
        rho, p_val =  scipy.stats.spearmanr(scaled_ranked_norm_rmsd, all_modes)
        # The smaller the p-value is, the better (evidence agains the hypothesis that variables are uncorrelated)
        #save_result(  "Does the (ANM) norm. RMSD depend on the frequency of the mode?\t %.3f\t (%.3f)"%(rho,p_val))
        
        #--------------------------------------------        
        # Energy increment (normalized) vs Frequency
        #--------------------------------------------
        all_norm_energies = []
        for key in norm_energy: all_norm_energies.extend(norm_energy[key])
        sns.regplot(numpy.array(all_freqs), numpy.array(all_norm_energies))
#         plt.show()
        ranked_norm_energies =  scipy.stats.rankdata(all_norm_energies)
        scaled_ranked_norm_energies = ranked_norm_energies*10./numpy.max(ranked_norm_energies)
        rho, p_val =  scipy.stats.spearmanr(ranked_norm_energies, all_modes)
        #save_result( "Does the (ANM) norm. energy increment depend on the frequency of the mode?\t %.3f\t (%.3f)"%(rho,p_val))
        
        #-----------------------------------------------------------
        # Time per step (cont. dep.) and p1, p2 (cat. ranked ind.)
        #-----------------------------------------------------------
        # rank them
        for p in [p1,p2]:
            ranked_time = scipy.stats.rankdata(numpy.array(all_data["time_per_step"]).astype(float))
            ranked_p = scipy.stats.rankdata(numpy.array(all_data[p]))
            # norm range of time to range of p
            ranked_time *= numpy.max(ranked_p)/numpy.max(ranked_time)
            rho, p_val =  scipy.stats.spearmanr(ranked_time, ranked_p)
            save_result( "Does the time per (ANM) step depend on %s?\t%.3f\t%.3f"%(p,rho,p_val))
        
        
        #-----------------------------------------------------------
        # RMSD (cont. dep.) and p1, p2 (cat. ranked ind.)
        #-----------------------------------------------------------
        #in this case rmsd is not normed
        for p in [p1,p2]:
            ranked_p = scipy.stats.rankdata(numpy.array(all_data[p]))
            ranked_rmsd = scipy.stats.rankdata(numpy.array(all_data[RMSD_LABEL]))
            scaled_ranked_rmsd = ranked_rmsd*numpy.max(ranked_p)/numpy.max(ranked_rmsd)
            rho, p_val =  scipy.stats.spearmanr(scaled_ranked_rmsd, ranked_p)
            save_result( "Does the (ANM) RMSD depend on %s?\t%.3f\t%.3f"%(p,rho,p_val))
        
        #-----------------------------------------------------------
        # Energy increment (cont. dep.) and p1, p2 (cat. ranked ind.)
        #-----------------------------------------------------------
        for p in [p1,p2]:
            ranked_p = scipy.stats.rankdata(numpy.array(all_data[p]))
            ranked_energies = scipy.stats.rankdata(numpy.array(all_data[ENERGY_LABEL]))
            scaled_ranked_energies = ranked_energies*numpy.max(ranked_p)/numpy.max(ranked_energies)
            rho, p_val =  scipy.stats.spearmanr(scaled_ranked_energies, ranked_p)
            save_result( "Does the (ANM) Energy increment depend on %s?\t%.3f\t%.3f"%(p,rho,p_val))
        
        # Energy absolute values
        smaller_than_0 = []
        bigger_than_0 = []
        for energy in energy_with_outlayers:
            if energy <= 0:
                smaller_than_0.append(energy)
            else:
                bigger_than_0.append(energy)
                
        save_result( "[Energy] Number of <0: %.2f Mean of > 0: %.2f (%.2f)"%(100*float(len(smaller_than_0)) /len(all_data[ENERGY_LABEL]), 
                                                                    numpy.mean(bigger_than_0), numpy.std(bigger_than_0)))
        
        save_result( "Avg. time per step %.2f (%.2f)"%(numpy.mean(all_data["time_per_step"]), numpy.std(all_data["time_per_step"])))
        
        other_results.close()
