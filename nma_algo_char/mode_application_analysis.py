"""
Created on Nov 19, 2015

@author: victor
"""
import numpy
import os
from min_bias.min_bias_calc import calculate_rmsds
from optparse import OptionParser
from nma_algo_char.common import load_control_json, pair_parameter_values,\
    parameter_value_to_string, create_directory, LineCounter, LineParser
from collections import defaultdict
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt

def load_data(data_folder, e_before, e_after,  coords_before, coords_after, step_time):
    data = {}
    # energies
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
        # in CC the main picker will otput the message once, then the NMA IC pickers. We need to eliminate the first ones.
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
        
        # skip first frame (usually an outlayer)
        print "DBG", min_len-1, len(modes), len(process_modes(experiment_details["prefix"], modes, 10))
        all_data["inc_U"].extend(process_energy_differences(raw_data)[1:])
        all_data["rmsd"].extend(process_after_perturb_rmsd(raw_data)[1:])
        all_data["mode"].extend(process_modes(experiment_details["prefix"], modes, 10)[1:min_len])
        all_data["time_per_step"].extend(raw_data["time_per_step"][1:min_len])
        all_data[p1].extend([v1]*(min_len-1))
        all_data[p2].extend([v2]*(min_len-1))
        
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
           
    db = DataFrame.from_dict(all_data, orient="index")
    db.transpose().to_csv(os.path.join(options.results_folder,os.path.basename(workspace),"data.csv"))
    db.columns = db.columns.get_level_values(0)
    
    if options.do_plots:
        import seaborn as sns
        g = sns.FacetGrid(db.transpose(), col=p1, row=p2, hue="mode", 
                          sharex = True, sharey = True, margin_titles=True,
                          legend_out = True)
        g.map(plt.scatter,  "rmsd", "inc_U")
        g.add_legend()
#         plt.show()
        g.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_u_rmsd.svg"))
        plt.close()
        
        # Do rmsd/energy time
        db = db.T
        db["rmsd_inc_U"] = db.T.loc["rmsd"] / db.T.loc["inc_U"]
        sns.lmplot("time_per_step", "rmsd_inc_U", db, hue= p1, fit_reg = False, scatter_kws={"alpha":0.6})
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_u_rmsd_vs_time_hue_%s.svg"%p1))
        sns.lmplot("time_per_step", "rmsd_inc_U", db, hue= p2, fit_reg = False, scatter_kws={"alpha":0.6})
        plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_u_rmsd_vs_time_hue_%s.svg"%p2))
        