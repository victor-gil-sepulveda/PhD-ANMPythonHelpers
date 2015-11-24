'''
Created on Nov 21, 2015

@author: victor
'''
import os
import copy
import ccvsic
from multiprocessing.pool import Pool
from ccvsic.common import set_parameter_value, pair_parameter_values,\
    wait_for_results_and_close, run_pele_in_folder, change_output_path_parameters,\
    prepare_workspace,load_control_json, parameter_value_to_string
from optparse import OptionParser
    
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-e", dest="exp_details_path")
    parser.add_option("--greasy", dest="greasy")
    parser.add_option("--mpi", action="store_true", default=False, dest="mpi")
    
    (options, args) = parser.parse_args()
    
    experiment_details = load_control_json(options.exp_details_path)
    control_file_template = load_control_json(experiment_details["control_file_template"])
    
    prepare_workspace(experiment_details)
    set_parameter_value("Initialization.Complex.files.path", control_file_template, 
                        os.path.join(experiment_details["workspace"],os.path.basename(experiment_details["initial_structure"])))
    
    if "license" in experiment_details:
        control_file_template["licenseDirectoryPath"] =  experiment_details["license"] 
    
    pool = Pool(experiment_details["number_of_processes"])
    results = []
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        folder_name = "%s_%s_%s_%s_%s"%(experiment_details["prefix"],
                                            experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                            experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
        control_dict = copy.deepcopy(control_file_template)
        folder = os.path.join(experiment_details["workspace"], folder_name)
        change_output_path_parameters(control_dict, experiment_details["common_changeable_paths"], folder)
        for p,v in [(p1,v1),(p2,v2)]:
            set_parameter_value(experiment_details["parameter_paths"][p], control_dict, v)
        
        #control_file_dict, folder, experiment_data, test = False, sleep_time = 10, return_command = False
        if options.greasy is None:
            results.append(pool.apply_async(run_pele_in_folder, (control_dict, folder, experiment_details, False, 5, False)))
        else:
            if not options.mpi:
                results.append("cd %s;"%(experiment_details["workspace"]) + run_pele_in_folder(control_dict, folder, experiment_details, True, 0, True))
            else:
                results.append("cd %s;"%(experiment_details["workspace"]) + 
                               "mpirun -np %s "%(experiment_details["number_of_processes"])+
                               run_pele_in_folder(control_dict, folder, experiment_details, True, 0, True))
    if options.greasy is None:
        wait_for_results_and_close(pool, results, 60)
    else:
        open(options.greasy, "w").write("\n".join(results)+"\n")
            
    
