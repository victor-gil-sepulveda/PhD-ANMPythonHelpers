'''
Created on Feb 21, 2016

@author: victor
'''
from optparse import OptionParser
from nma_algo_char.common import load_control_json, pair_parameter_values,\
    parameter_value_to_string, MetropolisMCSimulator
import os.path
from nma_algo_char.data_retrieval import load_cc_data, load_ic_data,\
    process_energy_differences, process_after_perturb_rmsd
import numpy

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-e", dest="experiment")
    parser.add_option("-w", dest="workspace")
    parser.add_option("-f", dest="folder")
    
    (options, args) = parser.parse_args()
    
    if not options.experiment:  
        parser.error('Experiment file not given')
    
    
    experiment_details = load_control_json(options.experiment)
    if options.workspace is not None:
        workspace = os.path.normpath(options.workspace)
    else:
        workspace = os.path.normpath(experiment_details["workspace"])
    
    print "RESULT\tMean\tSigma\tRMSD\t\tAcceptance\t\t"
    for (p1,v1),(p2,v2) in pair_parameter_values(experiment_details["check"], experiment_details["parameter_values"]):
        
        folder_name = "%s_%s_%s_%s_%s"%(experiment_details["prefix"],
                                        experiment_details["parameter_abbv"][p1], parameter_value_to_string(v1),
                                        experiment_details["parameter_abbv"][p2], parameter_value_to_string(v2))
        dataAvailable = True
        try:
            if experiment_details["prefix"] == "CC":
                raw_data, data_len = load_cc_data(os.path.join(workspace, 
                                                              folder_name,
                                                              options.folder),
                                                 full_pele_energy = True,
                                                 skip_first = 5)
            
            if experiment_details["prefix"] == "IC":
                raw_data, data_len = load_ic_data(os.path.join(workspace, 
                                                              folder_name,
                                                              options.folder),
                                                 skip_first = 5)
        except:
            print "Data could not be loaded. Skipping %s = %.3f  %s = %.3f "%(p1,v1,p2,v2)
            dataAvailable = False
         
        if dataAvailable:    
            energy_increments = process_energy_differences(raw_data)
            rmsd_increments = process_after_perturb_rmsd(raw_data)
            acc_mean_and_avg = MetropolisMCSimulator(energy_increments).perform_simulation(
                                                        min(100,len(energy_increments)), 20, 300)
            
            mean = (float(v1)+float(v2))/2.
            sigma = (float(v1)-float(v2))/4.
            
            print "RESULT\t%.3f\t(%.3f)\t%.3f\t(%.3f)\t%.3f\t(%.3f)"%(mean, sigma, 
                            numpy.mean(rmsd_increments), numpy.std(rmsd_increments),
                            acc_mean_and_avg[0],acc_mean_and_avg[1]), v1,v2
        

