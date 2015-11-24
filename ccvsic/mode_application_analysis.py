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

def load_data(data_folder):
    data = {}
    data["e_after"] = numpy.loadtxt(os.path.join(data_folder,"perturb_energy_after.log")).T[1]
    data["e_before"] = numpy.loadtxt(os.path.join(data_folder,"perturb_energy_before.log")).T[1]
    
    data["coords_after"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder,"after_anm_cc.log")),0,1) # -> delete first column (index=
    data["coords_before"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder,"initial_cc.log")),0,1)
    
    # skip first datum
    # trimm data as there can be an excess of 1
    min_len = min(len(data["e_after"]), len(data["e_before"]), len(data["coords_after"]), len(data["coords_before"]))
    data["e_after"] = data["e_after"][1:min_len]
    data["e_before"] = data["e_before"][1:min_len]
    data["coords_after"] = data["coords_after"][1:min_len]
    data["coords_before"] = data["coords_before"][1:min_len]
    
    return data, min_len-1

def process_energy_differences(data):
    return data["e_after"] - data["e_before"]

def process_after_perturb_rmsd(data):
    return calculate_rmsds(data["coords_before"], data["coords_after"])


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

def get_accepted_from_report(report_path):
    """
    This one is ok to get the accepted steps from the NMA CC report.
    """
    accepted_steps = []
    for line in open(report_path):
        if line[0] != "#":
            parts = line.split()
            accepted_steps.append(int(parts[2]))
    step_is_accepted = numpy.array([False]*(accepted_steps[-1]+1))
    step_is_accepted[accepted_steps] = True
    # The first one is the original structure (step 0) , 
    # skip the first one too (skipped when retrieving data)
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
        raw_data, min_len = load_data(os.path.join(workspace, folder_name,"info"))
        all_data["inc_U"].extend(process_energy_differences(raw_data))
        all_data["rmsd"].extend(process_after_perturb_rmsd(raw_data))
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
    db.to_csv(os.path.join(options.results_folder,os.path.basename(workspace),"data.csv"))
    db.columns = db.columns.get_level_values(0)
    print db.index
    print db.columns
    
    import seaborn as sns
    g = sns.FacetGrid(db.transpose(), col=p1, row=p2, hue="accepted", margin_titles=True)
    g.map(plt.scatter,  "rmsd", "inc_U")
#     plt.show()
    print workspace
    print os.path.join(options.results_folder,os.path.basename(workspace)+".svg")
    plt.savefig(os.path.join(options.results_folder,os.path.basename(workspace)+".svg"))
        