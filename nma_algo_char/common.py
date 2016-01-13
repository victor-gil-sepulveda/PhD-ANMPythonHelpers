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
import numbers
import numpy
from anmichelpers.tools.tools import norm
from pyRMSD.RMSDCalculator import RMSDCalculator
import random
from math import exp
from _collections import defaultdict
import matplotlib.pyplot as plt

class LineParser:
    def __init__(self, tag, position, conversion, split_token = None):
        self.data = []
        self.tag = tag
        self.tag_len = len(tag)
        self.position = position
        self.split_token = split_token
        self.conversion_function = conversion
        
    def parse_line(self, line):
        if line[0:self.tag_len] == self.tag:
            try:
                self.data.append(self.conversion_function(line.split(self.split_token)[self.position]))
            except Exception, e:
                print line
                raise e
class LineCounter:
    def __init__(self, tag):
        self.data = []
        self.tag = tag
        self.counter = 0
        
    def parse_line(self, line):
        if self.tag in line:
            self.counter = self.counter + 1
            
class MetropolisMCSimulator:
    BOLTZMANN_KCAL_MOL = 0.001987207
    
    def __init__(self, energy_increments):
        self.energy_increments = energy_increments

    def expected_probability(self, temperature):
        beta = 1.0/(MetropolisMCSimulator.BOLTZMANN_KCAL_MOL * temperature)
        experiment_acceptance = []
        for inc_u in self.energy_increments:
            experiment_acceptance.append(min(1.0,max(0.0,exp(-(beta*inc_u)))))
        return numpy.mean(experiment_acceptance), numpy.std(experiment_acceptance)
    
    def perform_simulation(self, number_of_samples, number_of_tries, temperature):
        beta = 1.0/(MetropolisMCSimulator.BOLTZMANN_KCAL_MOL * temperature)
        experiment_acceptance = []
        for _ in range(number_of_tries):
            energy_bootstrap = random.sample(self.energy_increments, number_of_samples)
            num_accepted = 0
            for energy_diff in energy_bootstrap: 
                if energy_diff <= 0:
                    num_accepted += 1
                else:
                    prob = exp(-(beta*energy_diff))
                    if prob > random.random():
                        num_accepted += 1
            experiment_acceptance.append(float(num_accepted)/number_of_samples)
        return numpy.mean(experiment_acceptance), numpy.std(experiment_acceptance)

def pair_parameter_values(parameter_keys, parameters):
    key1, key2 = parameter_keys[0], parameter_keys[1]
    for vals1 in parameters[key1]:
        for vals2 in parameters[key2]:
            yield ((key1, vals1),(key2,vals2))

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
        os.system("cp -r %s %s"%(base_data["PELE_docs"]["path"], 
                                           os.path.join(base_data["workspace"],"Documents")))
    if base_data["PELE_docs"]["action"] == "LINK":
        os.system("ln -s %s %s"%(base_data["PELE_docs"]["path"], 
                                           os.path.join(base_data["workspace"],"Documents")))

def run_pele_in_folder( control_file_dict, folder, experiment_data, test = False, sleep_time = 10, return_command = False):
    current_dir = os.getcwd()
    os.chdir(experiment_data["workspace"])
    create_directory(folder)
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
    
    os.chdir(current_dir)    
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

def parameter_value_to_string(val):
    if isinstance(val, basestring):
        return val
    else:
        if isinstance(val, numbers.Integral):
            return "%d"%val
        else:
            return "%.2f"%val
        
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

def get_values_by_hue(val_x, val_y, hue_vals):
    # assert all 3 have the same length
    if len(val_x) == len(val_y) and len(val_y)== len(hue_vals):
        vals_by_hue = defaultdict(lambda: {"x":[],"y":[]})
        for i in range(len(val_x)):
            vals_by_hue[hue_vals[i]]["x"].append(val_x[i])
            vals_by_hue[hue_vals[i]]["y"].append(val_y[i])
        return vals_by_hue
    else:
        print "Error::get_values_by_hue arays have not the same length (%d, %d, %d)"%(len(val_x), len(val_y), len(hue_vals))
        exit(-1)

def scatter_plot_by_hue(x, y, hue, colors):
    vals_by_hue = get_values_by_hue(x, y, hue)
    for i,hue_key in enumerate(vals_by_hue):
        plt.scatter(vals_by_hue[hue_key]["x"],vals_by_hue[hue_key]["y"], label = str(hue_key), color = colors[i], alpha = 0.6)

def prepare_subplots(row_len, col_len):
    if row_len > 1 or col_len > 1:
        f, axes = plt.subplots( col_len, row_len, sharey='row', sharex='col')
        f.subplots_adjust(hspace=0.4, wspace=0.3 )
        f.set_size_inches(12, 12, forward=True)
    else:
        f = plt.gcf()
        axes = {(0,0): plt.gca()}
    return f, axes
      