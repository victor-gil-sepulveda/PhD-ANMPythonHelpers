"""
Created on Nov 19, 2015

@author: victor
"""
import numpy
import os
from min_bias.min_bias_calc import calculate_rmsds
from optparse import OptionParser
from nma_algo_char.common import load_control_json, pair_parameter_values,\
    parameter_value_to_string, create_directory, LineCounter, LineParser,\
    MetropolisMCSimulator, get_values_by_hue, scatter_plot_by_hue
from collections import defaultdict
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
import pandas as pd
from anmichelpers.parsers.pronmd import ProdyNMDParser
import math

def load_data(data_folder, e_before, e_after,  coords_before, coords_after, step_time):
    data = {}
    # energies
    print "e_after folder" , os.path.join(data_folder,e_after)
    data["e_after"] = numpy.loadtxt(os.path.join(data_folder,e_after)).T[1]
    data["e_before"] = numpy.loadtxt(os.path.join(data_folder,e_before)).T[1]
    
    # coordinates
    data["coords_after"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder, coords_after)),0,1) # -> delete first column (index=
    data["coords_before"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder, coords_before)),0,1)
    
    # time per step
    data["time_per_step"] = numpy.loadtxt(os.path.join(data_folder, step_time)).T[1]
    
    # trimm data as there can be an excess of 1
    min_len = min(len(data["e_after"]), len(data["e_before"]), len(data["coords_after"]), len(data["coords_before"]))
    for key in data:
        data[key] = data[key][:min_len]
    
    return data, min_len

def load_ic_data(data_folder):
    return load_data(data_folder, "ener_mc_move_before.log", "ener_mc_move_after.log",  
                     "ca_mc_move_before.log", "ca_mc_move_after.log", "step_time.log")
    
def load_cc_data(data_folder):
    return load_data(data_folder, "perturb_energy_before.log", "perturb_energy_after.log",  
                     "initial_cc.log", "after_anm_cc.log", "step_time.log")

def process_energy_differences(data):
    return data["e_after"] - data["e_before"]

def process_after_perturb_rmsd(data):
    return calculate_rmsds(data["coords_before"], data["coords_after"])

def def_get_data_from_output_file(out_path):
    out_file = open(out_path)
    timesConvergedLineCounter = LineCounter("DBG: RELAX CLASH GRAD has converged")
    relaxIterationsParser = LineParser("DBG: RELAX CLASH GRAD - Iterations performed:", 7, int)
    for line in out_file:
        timesConvergedLineCounter.parse_line(line)
        relaxIterationsParser.parse_line(line)
    out_file.close()
    return relaxIterationsParser.data, timesConvergedLineCounter.counter

def get_data_from_log(log_path):
    PELEStepTimeParser = LineParser("Total step time:", 3, float)
    PickedModeParser = LineParser("Picked mode:", 2, int)
    for line in open(log_path):
        PELEStepTimeParser.parse_line(line)
        PickedModeParser.parse_line(line)
    return PELEStepTimeParser.data , PickedModeParser.data

def process_modes(NMA_type, modes, num_iterations = 10):
    if NMA_type == "CC":
        return modes
    if NMA_type == "IC":
        #print "DBG", "IC"
        # in CC the main picker will output the message once, then the NMA IC pickers. We need to eliminate the first ones.
        return numpy.array(modes)[(numpy.arange(len(modes))%(num_iterations+1) != 0)]

def save_relax_iterations(path, relax_iterations):
    open(path,"w").writelines( [" ".join([parameter_value_to_string(val) for val in item]) +"\n" for item in relax_iterations])

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("--workspace", dest="workspace")
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
    
    # create a place for the results
    create_directory(os.path.join(options.results_folder, os.path.basename(workspace)))
    all_data = defaultdict(list)
    relax_iterations = []
    acceptances = defaultdict(lambda: defaultdict(list))
    avg_energy = defaultdict(list)
    std_energy = defaultdict(list)
    avg_rmsd = defaultdict(list)
    avg_time = defaultdict(list)
    
    ENERGY_LABEL = "$\Delta$ U"
    RMSD_LABEL = "RMSD"
    nmd_file_name = {"CC":"normalized_modes.1.nmd", "IC":"normalized_modes_cc.1.nmd"}
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        folder_name = "%s_%s_%s_%s_%s"%(experiment_details["prefix"],
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
        
        pele_step_time, modes = get_data_from_log(os.path.join(workspace, folder_name, "log.txt"))
        
        ## CAUTION: HARDCODED FOLDER ('info'). Must be extracted using the control file
        if experiment_details["prefix"] == "CC":
            raw_data, min_len = load_cc_data(os.path.join(workspace, folder_name,"info"))
        
        if experiment_details["prefix"] == "IC":
            raw_data, min_len = load_ic_data(os.path.join(workspace, folder_name,"info"))
        
        _, eval, header =  ProdyNMDParser.read(os.path.join(workspace, 
                                                    folder_name,
                                                    "info",
                                                    nmd_file_name[experiment_details["prefix"]]))
        frequencies = numpy.sqrt(eval)/ (2*math.pi)
        mode_to_freq = dict(zip(range(len(frequencies)), frequencies))
        
        # skip first frame (usually an outlayer)
        #print "DBG", min_len-1, len(modes), len(process_modes(experiment_details["prefix"], modes, 10))
        energy_increments = process_energy_differences(raw_data)[1:]
        all_data[ENERGY_LABEL].extend(energy_increments)
        rmsd_increments = process_after_perturb_rmsd(raw_data)[1:]
        all_data[RMSD_LABEL].extend(rmsd_increments)
        modes = process_modes(experiment_details["prefix"], modes, 10)[1:min_len]
        all_data["Mode"].extend(modes)
        all_data["Mode_Freq"].extend([mode_to_freq[m] for m in modes])
        all_data["time_per_step"].extend(raw_data["time_per_step"][1:min_len])
        all_data[p1].extend([v1]*(min_len-1))
        all_data[p2].extend([v2]*(min_len-1))
        
        mc = MetropolisMCSimulator(energy_increments)
        acceptances[v1,v2] = mc.perform_simulation(min(100,len(energy_increments)), 20, 300)
        avg_energy[v1,v2] = numpy.mean(energy_increments)
        std_energy[v1,v2] = numpy.std(energy_increments)
        avg_rmsd[v1,v2] = (numpy.mean(rmsd_increments),numpy.std(rmsd_increments))
        avg_time[v1,v2] = (numpy.mean(raw_data["time_per_step"][1:min_len]),
                           numpy.std(raw_data["time_per_step"][1:min_len]))
        
        ## CAUTION: HARDCODED FILE ('out'). Must be extracted using the control file or experiment
        if experiment_details["prefix"] == "IC":
            relax_iters, times_conv = def_get_data_from_output_file(os.path.join(workspace, folder_name, "out"))
            if times_conv != 0:
                relax_iterations.append((v1,v2, numpy.mean(relax_iters[1:]), 
                                     numpy.std(relax_iters[1:]), 
                                     float(times_conv) / len(relax_iters[1:])))
    
    if experiment_details["prefix"] == "IC" and relax_iterations != []:
        save_relax_iterations(os.path.join(options.results_folder,os.path.basename(workspace),"relax.txt"),
                              relax_iterations)
           
    # Find limits for energy (has outliers)
    ener_low = min(all_data[ENERGY_LABEL])
    ener_high = max(sorted(all_data[ENERGY_LABEL])[:int(len(all_data[ENERGY_LABEL])*0.98)]) #98% should eliminate outlayers
    all_data[ENERGY_LABEL] = numpy.array(all_data[ENERGY_LABEL])
    all_data[ENERGY_LABEL][all_data[ENERGY_LABEL] > ener_high] = ener_high
    all_data[ENERGY_LABEL] = list(all_data[ENERGY_LABEL])
    
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
        g.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_u_rmsd.svg"))
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
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_vs_timespep.svg"))
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
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_rmsd_vs_tps_hue_acc.svg"))
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
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_energy_vs_tps_hue_acc.svg"))
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
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_rmsd_energy_vs_acc.svg"))
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
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_rmsd_vs_acc.svg"))
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
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_energy_vs_acc.svg"))
        plt.close()
        
        # Displacement (normalized) vs Frequency
            