"""
Created on Nov 19, 2015

@author: victor
"""
import numpy
import os
from min_bias.min_bias_calc import calculate_rmsds
from optparse import OptionParser
from ccvsic.common import load_control_json, pair_parameter_values,\
    parameter_value_to_string, create_directory
from collections import defaultdict
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
from anmichelpers.tools.tools import norm
from pyRMSD.RMSDCalculator import RMSDCalculator

def load_data(data_folder, e_before, e_after,  coords_before, coords_after):
    data = {}
    data["e_after"] = numpy.loadtxt(os.path.join(data_folder,e_after)).T[1]
    data["e_before"] = numpy.loadtxt(os.path.join(data_folder,e_before)).T[1]
    
    data["coords_after"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder, coords_after)),0,1) # -> delete first column (index=
    data["coords_before"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder, coords_before)),0,1)
    
    # skip first datum
    # trimm data as there can be an excess of 1
    min_len = min(len(data["e_after"]), len(data["e_before"]), len(data["coords_after"]), len(data["coords_before"]))
    data["e_after"] = data["e_after"][1:min_len]
    data["e_before"] = data["e_before"][1:min_len]
    data["coords_after"] = data["coords_after"][1:min_len]
    data["coords_before"] = data["coords_before"][1:min_len]
    
    return data, min_len-1

def load_ic_data(data_folder):
    return load_data(data_folder, "ener_mc_move_before.log", "ener_mc_move_after.log",  
                     "ca_mc_move_before.log", "ca_mc_move_after.log")
    
def load_cc_data(data_folder):
    return load_data(data_folder, "perturb_energy_before.log", "perturb_energy_after.log",  
                     "initial_cc.log", "after_anm_cc.log")

def process_energy_differences(data):
    return data["e_after"] - data["e_before"]

def process_after_perturb_rmsd(data):
    return calculate_rmsds(data["coords_before"], data["coords_after"])

def process_after_perturb_max_and_mean_disp(data):
    number_of_sets = len(data["coords_before"])
    num_coords = len(data["coords_before"][0])
    
    coordsets_before = numpy.array(data["coords_before"])
    coordsets_before = numpy.reshape(coordsets_before, (number_of_sets, num_coords/3, 3))
    coordsets_after = numpy.array(data["coords_after"])
    coordsets_after = numpy.reshape(coordsets_after, (number_of_sets, num_coords/3, 3))
    
    superimposed_translations = []
    for i in range(number_of_sets):
        coords = numpy.array([coordsets_before[i], coordsets_after[i]])
        
        calculator = RMSDCalculator(calculatorType = "QTRFIT_OMP_CALCULATOR",
                                    fittingCoordsets = coords)
        _, rearranged_coords = calculator.oneVsTheOthers(0, get_superposed_coordinates = True)
        superimposed_translations.append(rearranged_coords[1]-rearranged_coords[0])
    translations = numpy.array(superimposed_translations)
    norms = numpy.array([norm(t) for t in translations])
    return numpy.max(norms, axis = 1), numpy.mean(norms, axis = 1)

def def_get_data_from_output_file(out_path):
    """
    For NMA IC Output will have lines like:
    
    DBG: RELAX CLASH GRAD has converged
    
    DBG: RELAX CLASH GRAD - Iterations performed: 64

    DBG: - ANM Internal test - last: -14375.2036971242 current: -14377.42517911719 accepted: true
    """
    out_file = open(out_path)
    times_converged = 0
    step_is_accepted = []
    relax_iterations = []
    current_step = 0
    for line in out_file:
        if line[0:3]=="DBG":
            if line[0:24] == "DBG: - ANM Internal test":
                parts = line.split()
                if parts[11] == "true":
                    step_is_accepted.append(True)
                else:
                    step_is_accepted.append(False)
                current_step = current_step + 1
            if line[0:35] == "DBG: RELAX CLASH GRAD has converged":
                times_converged += 1
            if line[0:45] == "DBG: RELAX CLASH GRAD - Iterations performed:":
                parts = line.split(":")
                relax_iterations.append(int(parts[-1]))
    return numpy.array(step_is_accepted), relax_iterations, times_converged


def get_cc_time_from_log(log_path):
    times = []
    for line in open(log_path):
        if line[:7] == "ANM Ef:":
            parts = line.split()
            time = float(parts[8])
            times.append(time)
    return times[1:] # Always skip first structure

def get_ic_time_from_log(log_path):
    times = []
    for line in open(log_path):
        if line[:15] == "ANM_IC_Cycle t:":
            parts = line.split()
            time = float(parts[2][:-1])
            times.append(time)
    return times[1:] # Always skip first structure

def get_accepted_from_report(report_path):
    """
    This one is ok to get the accepted steps from the NMA CC report.
    """
    accepted_steps = []
    for line in open(report_path):
        if line[0] != "#":
            parts = line.split()
            accepted_steps.append(int(parts[1]))
    step_is_accepted = numpy.array([False]*(accepted_steps[-1]+1))
    step_is_accepted[accepted_steps] = True
    # The first one is the original structure (step 0) , 
    # skip the first one too (skipped when retrieving data)
    print "ACCEPTED", len(step_is_accepted[step_is_accepted==True]), len(step_is_accepted),accepted_steps[-1]+1
    return step_is_accepted[1:]

def pad_accepted(s, min_len):
    if len(s) == min_len:
        return s
    else:
        if len(s) > min_len:
            return s[:min_len]
        else:
            sn = numpy.resize(s, min_len)
            sn[min_len:] = False
            return sn

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--experiment", dest="experiment")
    parser.add_option("--results", dest="results_folder")
    parser.add_option("--workspace", dest="workspace")
    (options, args) = parser.parse_args()
    
    if not options.experiment:  
        parser.error('Experiment file not given')
    
    if not options.experiment:  
        parser.error('Results folder not given')
    
    experiment_details = load_control_json(options.experiment)
    if options.workspace is not None:
        workspace = os.path.normpath(options.workspace)
    else:
        workspace = os.path.normpath(experiment_details["workspace"])
    
    # create a place for the results
    create_directory(os.path.join(options.results_folder, os.path.basename(workspace)))
    all_data = defaultdict(list)
    relax_iterations = defaultdict(lambda: defaultdict(dict))
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        folder_name = "%s_%s_%s_%s_%s"%(experiment_details["prefix"],
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
        
        ## CAUTION: HARDCODED FOLDER ('info'). Must be extracted using the control file
        if experiment_details["prefix"] == "CC":
            raw_data, min_len = load_cc_data(os.path.join(workspace, folder_name,"info"))
            time_per_step = get_cc_time_from_log(os.path.join(workspace, folder_name,"log.txt"))
        
        if experiment_details["prefix"] == "IC":
            raw_data, min_len = load_ic_data(os.path.join(workspace, folder_name,"info"))
            time_per_step =  numpy.array(get_ic_time_from_log(os.path.join(workspace, folder_name,"log.txt")))/10. # the number of iterations
         
#         all_data["time_per_step"].append(({p1:v1,p2:v2}, time_per_step))
            
        all_data["inc_U"].extend(process_energy_differences(raw_data))
        all_data["rmsd"].extend(process_after_perturb_rmsd(raw_data))
        max_disp, avg_disp = process_after_perturb_max_and_mean_disp(raw_data)
        all_data["max_disp"].extend(max_disp.tolist())
        all_data["avg_disp"].extend(avg_disp.tolist())
        
        process_after_perturb_max_and_mean_disp(raw_data)
        all_data[p1].extend([v1]*min_len)
        all_data[p2].extend([v2]*min_len)
        
        ## CAUTION: HARDCODED FILE ('report'). Must be extracted using the control file
        if experiment_details["prefix"] == "CC":
            acc = get_accepted_from_report(os.path.join(workspace, folder_name, "report"))
            all_data["accepted"].extend(pad_accepted(acc, min_len))
            
        ## CAUTION: HARDCODED FILE ('out'). Must be extracted using the control file or experiment
        if experiment_details["prefix"] == "IC":
            acc, relax_iters, times_conv = def_get_data_from_output_file(os.path.join(workspace, folder_name, "out"))
            all_data["accepted"].extend( pad_accepted(acc, min_len))
            relax_iterations[(v1,v2)]["relax_avg"] = numpy.mean(relax_iters)
            relax_iterations[(v1,v2)]["relax_std"] = numpy.std(relax_iters)
            relax_iterations[(v1,v2)]["times_conv"] = float(times_conv) / len(relax_iters)
            
    db = DataFrame.from_dict(all_data, orient="index")
    db.transpose().to_csv(os.path.join(options.results_folder,os.path.basename(workspace),"data.csv"))
    db.columns = db.columns.get_level_values(0)
    
    import seaborn as sns
    g = sns.FacetGrid(db.transpose(), col=p1, row=p2, hue="accepted", 
                      sharex = True, sharey = True, margin_titles=True,
                      legend_out = True)
    g.map(plt.scatter,  "rmsd", "inc_U")
    g.add_legend()
    #plt.show()
    g.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_u_rmsd.svg"))
    plt.close()

    plt.figure()
    g = sns.FacetGrid(db.transpose(), col=p1, row=p2, hue="accepted", 
                      sharex = True, sharey = True, margin_titles=True,
                      legend_out = True)
    g.map(plt.scatter,  "max_disp", "inc_U")
    g.add_legend()
    #plt.show()
    g.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_u_max_disp.svg"))
    plt.close()
    
    plt.figure()
    g = sns.FacetGrid(db.transpose(), col=p1, row=p2, hue="accepted", 
                      sharex = True, sharey = True, margin_titles=True,
                      legend_out = True)
    g.map(plt.scatter,  "avg_disp", "inc_U")
    g.add_legend()
    #plt.show()
    g.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+"_u_avg_disp.svg"))
    plt.close()
    
    
    #Time analysis 
    time_related = ["rmsg", "iters", "steeringForce"]
    param = p1 if p1 in time_related else p2
    
#     keys = sorted(list(set(all_data[param])))
#     time_per_key = {}
#     for params, values in all_data["time_per_step"]:
#         time_per_key[params[param]] = numpy.mean(values)
#     print time_per_key
        