"""
Created on Nov 21, 2015

@author: victor
"""

import os
import json
import time
import sys
import re
import errno

def pair_parameter_values(parameter_keys, parameters):
    key1, key2 = parameter_keys[0], parameter_keys[1]
    for vals1 in parameters[key1]:
        for vals2 in parameters[key2]:
            yield (key1, vals1),(key2,vals2)

def set_parameter_value(key_description, param_dict, value):
    keys = key_description.split(".")
    tmp_dic = param_dict
    for i, k in enumerate(keys[:-1]):
        val = tmp_dic[k]
        if isinstance(val, list):
            # Change in all the items (all items must have same structure)
            sub_search = ".".join(keys[i+1:])
            for item in val:
                try:
                    set_parameter_value(sub_search, item, value)
                except:
                    pass
            return
        else: 
            tmp_dic = val
        
    tmp_dic[keys[-1]] = value

def prepare_workspace(base_data):
    create_directory(base_data["workspace"])
    os.system("cp -r %s %s"%(base_data["initial_structure"], base_data["workspace"]))
    if base_data["PELE_data"]["action"] == "COPY":
        os.system("cp -r %s %s"%(base_data["PELE_data"]["path"], base_data["workspace"]))
    if base_data["PELE_data"]["action"] == "LINK":
        os.system("ln -s %s %s"%(base_data["PELE_data"]["path"], base_data["workspace"]))
    if base_data["PELE_docs"]["action"] == "COPY":
        os.system("cp -r %s %s"%(base_data["PELE_docs"]["path"], base_data["workspace"]))
    if base_data["PELE_docs"]["action"] == "LINK":
        os.system("ln -s %s %s"%(base_data["PELE_docs"]["path"], base_data["workspace"]))

def run_pele_in_folder( control_file_dict, folder, experiment_data, test = False, sleep_time = 10, return_command = False):
    create_directory(folder)
    os.system("ln -s %s/Data %s"%(experiment_data["workspace"], folder ))
    os.system("ln -s %s/Documents %s"%(experiment_data["workspace"], folder ))
    create_directory(os.path.join(folder,"info"))
    control_file_path = os.path.join(folder,'control.json')
    out_file_path = os.path.join(folder,'out')
    with open(control_file_path, 'w') as outfile:
        json.dump(control_file_dict, outfile, indent=4)
    
    cmd = "%s %s > %s"%(experiment_data["PELE_exec_cmd"],control_file_path,out_file_path)
    
    print cmd
    sys.stdout.flush()

    if not test:
        os.system(cmd)
    else:
        time.sleep(sleep_time)
        
    if return_command:
        return cmd
    
def change_output_path_parameters(control_dict, params_dict, folder):
    for complete_key in params_dict:
        set_parameter_value(complete_key, control_dict, os.path.join(folder, params_dict[complete_key]))

def wait_for_results_and_close(pool, results, query_time):
    finished = 0
    while finished != len(results):
        finished = 0
        for result in results:
            if result.ready():
                finished += 1
        time.sleep(query_time)
    print "All processes have finished"
    pool.close()
    
def load_control_json(json_script):
    json_string = remove_comments(open(json_script).read())
    return convert_to_utf8(json.loads(json_string))

def remove_comments(string):
    """
    Removes /**/ and // comments from a string (used with the control script).
    From http://stackoverflow.com/questions/2319019/using-regex-to-remove-comments-from-source-files
    """
    string = re.sub(re.compile("/\*.*?\*/",re.DOTALL ) ,"" ,string) # remove all occurance streamed comments (/*COMMENT */) from string
    string = re.sub(re.compile("//.*?\n" ) ,"" ,string) # remove all occurance singleline comments (//COMMENT\n ) from string
    return string

def convert_to_utf8(my_input):
    """
    Recursively encodes all strings of an input dictionary as UTF-8. Useful to eliminate unicode strings.

    @param my_input: A dictionary object.

    @return: Encoded dictionary.
    """
    if isinstance(my_input, dict):
        return {convert_to_utf8(key): convert_to_utf8(value) for key, value in my_input.iteritems()}
    elif isinstance(my_input, list):
        return [convert_to_utf8(element) for element in my_input]
    elif isinstance(my_input, unicode):
        return my_input.encode('utf-8')
    else:
        return my_input
    
def create_directory(directory_path, ensure_writability = False):
    """
    Creates a directory (with subdirs) if it doesn't exist.
    
    @param directory_path: the path of the directory and subdirectories to be created. 
    """
    if ensure_writability:
        if not os.access(os.path.dirname(directory_path), os.W_OK):
            return False
    
    try:
        os.makedirs(directory_path)
        return True
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
    
    return False

