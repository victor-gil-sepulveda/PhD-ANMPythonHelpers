'''
Created on 14/9/2015

@author: victor
'''
import sys
from pyRMSD.RMSDCalculator import RMSDCalculator
import numpy


def calculate_rmsds(coord_set_ref, coord_set):
    number_of_steps = min(coord_set_ref.shape[0],coord_set.shape[0])
    number_of_coords = coord_set_ref.shape[1]
    
    print coord_set_ref.shape, coord_set.shape
    
    rmsds = []
    for i in range(number_of_steps):
        ref = numpy.reshape(coord_set_ref[i], (number_of_coords/3, 3))
        conf = numpy.reshape(coord_set[i], (number_of_coords/3, 3))
        
        coords = numpy.array([ref, conf])
        
        calculator = RMSDCalculator(calculatorType = "QTRFIT_OMP_CALCULATOR",
                                    fittingCoordsets = coords)
        rmsds.append(calculator.oneVsFollowing(0)[0])
    
    return numpy.array(rmsds)

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
    
    
    
    
    