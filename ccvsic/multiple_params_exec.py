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
    prepare_workspace,load_control_json
    
if __name__ == '__main__':
    data_path = ccvsic.__path__[0]
    test = os.path.join(data_path, "data", "anm_step_charachterization", "cc", "experiment.json")
    experiment_details = load_control_json(test)
    control_file_template = load_control_json(experiment_details["control_file_template"])
    
    prepare_workspace(experiment_details)
    set_parameter_value("Initialization.Complex.files.path", control_file_template, 
                        os.path.join(experiment_details["workspace"],os.path.basename(experiment_details["initial_structure"])))
    
    pool = Pool(experiment_details["number_of_processes"])
    results = []
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        folder_name = "%s_%s_%.2f_%s_%.2f"%(experiment_details["prefix"],
                                            experiment_details["parameter_abbv"][p1],v1,
                                            experiment_details["parameter_abbv"][p2],v2)
        control_dict = copy.deepcopy(control_file_template)
        folder = os.path.join(experiment_details["workspace"], folder_name)
        change_output_path_parameters(control_file_template, experiment_details["common_changeable_paths"], folder)
        for p,v in [(p1,v1),(p2,v2)]:
            set_parameter_value(experiment_details["parameter_paths"][p], control_file_template, v)
        
        #control_file_dict, folder, experiment_data, test = False, sleep_time = 10, return_command = False
        results.append(pool.apply_async(run_pele_in_folder, (control_dict, folder, experiment_details, False, 5, False)))
    
    wait_for_results_and_close(pool, results, 60)
    
