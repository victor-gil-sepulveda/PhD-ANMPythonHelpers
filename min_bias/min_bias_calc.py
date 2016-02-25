'''
Created on 14/9/2015

@author: victor
'''
import sys
import numpy
from nma_algo_char.data_retrieval import calculate_rmsds

if __name__ == '__main__':
    current_coords = numpy.delete(numpy.loadtxt(sys.argv[1]),0,1)
    proposal_coords = numpy.delete(numpy.loadtxt(sys.argv[2]),0,1)
    anm_final_coords = numpy.delete(numpy.loadtxt(sys.argv[3]),0,1)
    minim_final_coords = numpy.delete(numpy.loadtxt(sys.argv[4]),0,1)
    
    executions = {
                  "curr_vs_prop":(current_coords, proposal_coords),
                  "curr_vs_anmf":(current_coords, anm_final_coords),
                  "curr_vs_minim":(current_coords, minim_final_coords),
                  "prop_vs_anmf":(proposal_coords, anm_final_coords),
                  "prop_vs_minim":(proposal_coords, minim_final_coords),
                  "anmf_vs_minim":(anm_final_coords, minim_final_coords)
    }
    
    for execution_label in executions:
        print "Calculating ",execution_label
        reference = executions[execution_label][0]
        conf = executions[execution_label][1]
        numpy.savetxt(execution_label, calculate_rmsds(reference, conf))
    
    
    
    
    