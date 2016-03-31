import numpy
from pyRMSD.RMSDCalculator import RMSDCalculator
import os.path
import math
from anmichelpers.parsers.pronmd import ProdyNMDParser

def load_data_from_multiple_processors(sim_type, 
                                       num_procs, 
                                       data_folder, 
                                       skip_first,  
                                       max_samples = numpy.inf, 
                                       full_pele_energy = False):
    
    # processors from 1 to num_procs-1 have data
    all_raw_data = {}
    min_lens = []
    for proc in range(1,num_procs):
        if sim_type == "CC":
            raw_data, min_len = load_cc_data_from_proc(data_folder,
                                                        proc,
                                                        skip_first,
                                                        max_samples,
                                                        full_pele_energy)
        if sim_type == "IC":
            raw_data, min_len = load_ic_data_from_proc(data_folder, 
                                                        proc,
                                                        skip_first,
                                                        max_samples)
        min_lens.append(min_len)
        
        for key in raw_data:
            if key in all_raw_data:
                all_raw_data[key].extend(raw_data[key])
            else:
                all_raw_data[key] = list(raw_data[key])
    
    for key in all_raw_data:
        all_raw_data[key] = numpy.array(all_raw_data[key])
        print len(all_raw_data[key])
    
    return all_raw_data, min_lens


def load_single_proc_data(sim_type, 
                          data_folder, 
                          skip_first, 
                          max_samples = numpy.inf,
                          full_pele_energy = True):
    ## CAUTION: HARDCODED FOLDER ('info'). Must be extracted using the control file
    if sim_type == "CC":
        raw_data, min_len = load_cc_data(data_folder, skip_first, max_samples= max_samples, full_pele_energy=full_pele_energy)
    
    if sim_type == "IC":
        raw_data, min_len = load_ic_data(data_folder, skip_first, max_samples=max_samples)
    
    return raw_data, min_len


def add_coords_data(data, data_folder, filename, num_procs, min_len):
    # name = "coords_after_anm"
    data[filename] = []
    
    if num_procs > 1:
        for i, proc in enumerate(range(1,num_procs)):
            loaded_data = numpy.delete(numpy.loadtxt(os.path.join(data_folder, "%s.p%d.log"%(filename, proc))),0,1) 
            data[filename].extend(loaded_data[:min_len[i]]) 
    else:
        loaded_data = numpy.delete(numpy.loadtxt(os.path.join(data_folder, "%s.log"%filename )),0,1)
        data[filename].extend(loaded_data[:min_len])
   
    data[filename] = numpy.array(data[filename])
    return data

def add_time_like_data(data, data_folder, filename, num_procs, min_len, has_step_column = True):
    # name = "coords_after_anm"
    data[filename] = []
    
    if num_procs > 1:
        for i, proc in enumerate(range(1,num_procs)):
            if has_step_column:
                loaded_data = numpy.loadtxt(os.path.join(data_folder, "%s.p%d.log"%(filename, proc))).T[1]
            else:
                loaded_data = numpy.loadtxt(os.path.join(data_folder, "%s.log"%filename))
            data[filename].extend(loaded_data[:min_len[i]]) 
    else:
        if has_step_column:
            loaded_data = numpy.loadtxt(os.path.join(data_folder, "%s.p%d.log"%(filename, proc))).T[1]
        else:
            loaded_data = numpy.loadtxt(os.path.join(data_folder, "%s.log"%filename))
        data[filename].extend(loaded_data[:min_len])
    return data

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
        time_data = numpy.loadtxt(os.path.join(data_folder, step_time))
        if len(time_data.shape) > 1: # then it has the step column 
            data["time_per_step"] = numpy.loadtxt(os.path.join(data_folder, step_time)).T[1]
        else:
            data["time_per_step"] = numpy.loadtxt(os.path.join(data_folder, step_time))
            
    except IOError:
        data["time_per_step"] = [0]*len(data["e_before"])
        
    # chosen mode at each step
    if modes_selected is not None:
        try:
            data["modes"] = numpy.loadtxt(os.path.join(data_folder, modes_selected)).T[1]
        except IOError:
            data["modes"] = [0]*len(data["time_per_step"])
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
                     '''"step_time.log",''' "ic_movement_time.log",
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
    
# TODO change full_pele_energy by full_pele_step
def load_cc_data_from_proc(data_folder, processor, skip_first = 1, max_samples = numpy.inf, full_pele_energy = False):
    if full_pele_energy:
        final_energy_file = "final_energy.p%d.log"%processor
        final_coords_file = "after_minimization_cc.p%d.log"%processor
        final_time_file = "pele_step_time.p%d.log"%processor
    else:
        final_energy_file = "perturb_energy_after.p%d.log"%processor
        final_coords_file = "after_anm_cc.p%d.log"%processor
        final_time_file = "step_time.p%d.log"%processor
    
    return load_data(data_folder, 
                     "perturb_energy_before.p%d.log"%processor, 
                     final_energy_file,  
                     "initial_cc.p%d.log"%processor, 
                     final_coords_file, 
                     final_time_file,
                     None,
                     skip_first, 
                     max_samples)

def load_cc_data(data_folder, skip_first = 1, max_samples = numpy.inf, full_pele_energy = False):
    if full_pele_energy:
        final_energy_file = "final_energy.log"
        final_coords_file = "after_minimization_cc.log"
        final_time_file = "pele_step_time.log"
    else:
        final_energy_file = "perturb_energy_after.log"
        final_coords_file = "after_anm_cc.log"
        final_time_file = "step_time.log"
        
    return load_data(data_folder, 
                     "perturb_energy_before.log", 
                     final_energy_file,  
                     "initial_cc.log", 
                     final_coords_file, 
                     final_time_file,
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
