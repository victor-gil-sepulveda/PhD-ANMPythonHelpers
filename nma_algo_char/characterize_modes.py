'''
Created on Jan 28, 2016

@author: victor
'''
from optparse import OptionParser
from pyRMSD.utils.proteinReading import Reader
from anmichelpers.tools.tools import norm
import os.path
from anmichelpers.parsers.pronmd import ProdyNMDParser
from nma_algo_char.acceptance_and_rmsf_from_logs import load_single_proc_data
import numpy
from nma_algo_char.common import LineParser, MetropolisMCSimulator
from collections import Counter
import math
from anmichelpers.comparison.comparison import overlap
from nma_algo_char.mode_application_analysis import process_energy_differences

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--type", dest="sim_type")
    parser.add_option("-f",  dest="work_folder")
    parser.add_option("-o",  dest="open_pdb_file")
    parser.add_option("-c",  dest="closed_pdb_file")
    parser.add_option("-t", type="int",  dest="temperature")
    
    
    (options, args) = parser.parse_args()
    
    
    open_coords = Reader().readThisFile(options.open_pdb_file).gettingOnlyCAs().read()
    closed_coords = Reader().readThisFile(options.closed_pdb_file).gettingOnlyCAs().read()
   
    # Get eigenvectors, used modes and coordinates
    eigenfile = {"CC": "original_modes.1.nmd", "IC":"original_modes_cc.1.nmd"} 
    eigenvalues, eigenvectors, _ = ProdyNMDParser.read(os.path.join(options.work_folder,"info", eigenfile[options.sim_type]))
    raw_data, min_len = load_single_proc_data(options.sim_type, os.path.join(options.work_folder,"info"))
    mode_parser = LineParser("Picked", 2, int)
    for line in open(os.path.join(options.work_folder, "log.txt")):
        mode_parser.parse_line(line)
    raw_data["picked_mode"] =  numpy.array(mode_parser.data)[:min_len]
    frequencies = Counter(raw_data["picked_mode"])
    energy_increments = process_energy_differences(raw_data)[1:]
    mc = MetropolisMCSimulator(energy_increments)
    who_is_accepted = mc.who_is_accepted(options.temperature)
    acc_frequencies = Counter(raw_data["picked_mode"][who_is_accepted])
    
    # Calculate eigenfrequencies
    mode_info = {}
    for i,eigenvalue in enumerate(eigenvalues):
        mode_info[i] = {
                        "frequency": frequencies[i] / float(min_len),
                        "acc_frequency": acc_frequencies[i] / float(len(who_is_accepted)),
                        "eigen_frequency": math.sqrt(eigenvalue),
                        "avg_tr_overlap":[],
                        "acc_avg_tr_overlap":[]
                        }
    freq_sum = 0
    for i in range(len(mode_info.keys())):
        mode_info[i]["rel_frequency"] = mode_info[len(mode_info.keys())-1]["eigen_frequency"] / mode_info[i]["eigen_frequency"] 
        freq_sum += mode_info[i]["rel_frequency"]
    for i in range(len(mode_info.keys())):
        mode_info[i]["norm_frequency"] = mode_info[i]["rel_frequency"] / freq_sum
    
    
    # Calculate "real frequencies". I.e. the overlap of the chosen mode with the movement
    movements = raw_data["coords_after"]-raw_data["coords_before"]
    print len(movements), len(who_is_accepted)
    for i, movement in enumerate(movements):
        mode = raw_data["picked_mode"][i]
        mode_info[mode]["avg_tr_overlap"].append(overlap(eigenvectors[mode], movement))
        
    for i in who_is_accepted:
        movement = movements[i]
        mode = raw_data["picked_mode"][i]
        mode_info[mode]["acc_avg_tr_overlap"].append(overlap(eigenvectors[mode], movement))
    
    for mode in range(len(mode_info.keys())):
        mode_info[mode]["std_tr_overlap"] = numpy.std(mode_info[mode]["avg_tr_overlap"])
        mode_info[mode]["avg_tr_overlap"] = numpy.mean(mode_info[mode]["avg_tr_overlap"])
        mode_info[mode]["acc_std_tr_overlap"] = numpy.std(mode_info[mode]["acc_avg_tr_overlap"])
        mode_info[mode]["acc_avg_tr_overlap"] = numpy.mean(mode_info[mode]["acc_avg_tr_overlap"])
        mode_info[mode]["norm"] = norm(eigenvectors[mode])
    
    for mode in range(len(mode_info.keys())):
        print len(eigenvectors[mode]), len(closed_coords.flatten()-open_coords.flatten())
        mode_info[mode]["op_cl_overlap"] = overlap(eigenvectors[mode], closed_coords.flatten()-open_coords.flatten())
    
    
    ordered_tags = ["op_cl_overlap", "frequency", "acc_frequency", 
                    "eigen_frequency", "norm_frequency", "avg_tr_overlap",
                    "std_tr_overlap", "acc_avg_tr_overlap", "acc_std_tr_overlap", 
                    "norm"]
    
    print "Mode\t",
    for tag in ordered_tags:
        print tag,"\t",
    print
        
    for mode in sorted(mode_info.keys()):
        print mode,"\t", 
        for tag in ordered_tags:
            print "%.3f\t"%mode_info[mode][tag],
        print