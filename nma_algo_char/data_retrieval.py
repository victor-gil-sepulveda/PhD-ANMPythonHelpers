import numpy
from pyRMSD.RMSDCalculator import RMSDCalculator
import os.path
import math
from anmichelpers.parsers.pronmd import ProdyNMDParser

def load_energy(data_folder, energy_file):
    return numpy.loadtxt(os.path.join(data_folder,energy_file)).T[1]

def load_data(data_folder, e_before, e_after, coords_before, coords_after, 
              step_time, modes_selected, 
              skip_first = 1, max_samples = numpy.inf):
    data = {}

    # energies
    data["e_after"] = numpy.loadtxt(os.path.join(data_folder,e_after)).T[1]
    data["e_before"] = numpy.loadtxt(os.path.join(data_folder,e_before)).T[1]
    
    # coordinates
    data["coords_after"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder, coords_after)),0,1) # -> delete first column (the number of step)
    data["coords_before"] = numpy.delete(numpy.loadtxt(os.path.join(data_folder, coords_before)),0,1)
    
    # time per step
    try:
        data["time_per_step"] = numpy.loadtxt(os.path.join(data_folder, step_time)).T[1]
    except IOError:
        data["time_per_step"] = [0]*len(data["e_before"])
        
    # chosen mode at each step
    if modes_selected is not None:
        data["modes"] = numpy.loadtxt(os.path.join(data_folder, modes_selected)).T[1]
    else:
        data["modes"] = [0]*len(data["time_per_step"])
    
    # trim data as there can be an excess of 1
    min_len = min(max_samples+skip_first,
                  len(data["e_after"]), 
                  len(data["e_before"]), 
                  len(data["coords_after"]), 
                  len(data["coords_before"]),
                  len(data["time_per_step"]),
                  len(data["modes"])
                  )
    
    for key in data:
        data[key] = data[key][skip_first:min_len]
    
    return data, len(data["e_after"])

def load_ic_data(data_folder, skip_first = 1, max_samples = numpy.inf):
    return load_data(data_folder, 
                     "ener_mc_move_before.log", 
                     "ener_mc_move_after.log",  
                     "ca_mc_move_before.log", 
                     "ca_mc_move_after.log", 
                     "ic_movement_time.log", 
                     "ic_movement_mode.log",
                     skip_first, 
                     max_samples)

def load_ic_data_from_proc(data_folder, processor, skip_first = 1, max_samples = numpy.inf):
    return load_data(data_folder, 
                     "ener_mc_move_before.p%d.log"%processor, 
                     "ener_mc_move_after.p%d.log"%processor,  
                     "ca_mc_move_before.p%d.log"%processor, 
                     "ca_mc_move_after.p%d.log"%processor, 
                     "ic_movement_time.p%d.log"%processor, 
                     "ic_movement_mode.p%d.log"%processor,
                     skip_first, 
                     max_samples)
    

def load_cc_data_from_proc(data_folder, processor, skip_first = 1, full_pele_energy = False, max_samples = numpy.inf):
    if full_pele_energy:
        final_energy_file = "final_energy.p%d.log"%processor
    else:
        final_energy_file = "perturb_energy_after.p%d.log"%processor
    return load_data(data_folder, 
                     "perturb_energy_before.p%d.log"%processor, 
                     final_energy_file,  
                     "initial_cc.p%d.log"%processor, 
                     "after_anm_cc.p%d.log"%processor, 
                     "step_time.p%d.log"%processor,
                     None,
                     skip_first, 
                     max_samples)

def load_cc_data(data_folder, skip_first = 1, full_pele_energy = False, max_samples = numpy.inf):
    if full_pele_energy:
        final_energy_file = "final_energy.log"
    else:
        final_energy_file = "perturb_energy_after.log"
    return load_data(data_folder, 
                     "perturb_energy_before.log", 
                     final_energy_file,  
                     "initial_cc.log", 
                     "after_anm_cc.log", 
                     "step_time.log",
                     None,
                     skip_first,
                     max_samples)


def process_energy_differences(data):
    return data["e_after"] - data["e_before"]

def process_after_perturb_rmsd(data):
    return calculate_rmsds(data["coords_before"], data["coords_after"])


def calculate_rmsds(coord_set_ref, coord_set):
    number_of_steps = min(coord_set_ref.shape[0],coord_set.shape[0])
    number_of_coords = coord_set_ref.shape[1]
    
    print "Calculating RMSDs"
    print "\t- Number of steps:",number_of_steps
    print "\t- Number of coordinates:",number_of_coords
    print "\t- Shapes (ref/other):", coord_set_ref.shape, coord_set.shape
    
    rmsds = []
    for i in range(number_of_steps):
        ref = numpy.reshape(coord_set_ref[i], (number_of_coords/3, 3))
        conf = numpy.reshape(coord_set[i], (number_of_coords/3, 3))
        
        coords = numpy.array([ref, conf])
        
        calculator = RMSDCalculator(calculatorType = "QTRFIT_OMP_CALCULATOR",
                                    fittingCoordsets = coords)
        rmsds.append(calculator.oneVsFollowing(0)[0])
    
    return numpy.array(rmsds)

def get_mode_frequencies(modes, nmd_file):
    # returns a map from mode index to frequency
    evalues, _, _ =  ProdyNMDParser.read(nmd_file)
    frequencies = numpy.sqrt(evalues)/ (2*math.pi)
    mode_to_freq = dict(zip(range(len(frequencies)), frequencies))
    return numpy.array([mode_to_freq[m] for m in modes])
