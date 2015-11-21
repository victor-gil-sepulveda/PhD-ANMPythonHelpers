"""
Created on Nov 19, 2015

@author: victor
"""
import sys
import numpy
import os
from min_bias.min_bias_calc import calculate_rmsds

def load_data(data_folder):
    data = {}
    data["e_after"] = numpy.loadtxt(os.path.join(data_folder,"perturb_energy_after.log")).T[1]
    data["e_before"] = numpy.loadtxt(os.path.join(data_folder,"perturb_energy_before.log")).T[1][0:len(data["e_after"])] # -> trimmed as there can be an excess of 1
    
    data["coords_after"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder,"after_anm_cc.log")),0,1) # -> delete first column
    data["coords_before"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder,"initial_cc.log")),0,1)[0:len(data["coords_after"])]
    return data

def process_energy_differences(data):
    return data["e_after"] - data["e_before"]

def process_after_perturb_rmsd(data):
    return calculate_rmsds(data["coords_before"], data["coords_after"])

if __name__ == '__main__':
    data_folder = sys.argv[1]
    
    data = load_data(data_folder)
    
    print numpy.mean(process_energy_differences(data))
    print numpy.mean(process_after_perturb_rmsd(data))