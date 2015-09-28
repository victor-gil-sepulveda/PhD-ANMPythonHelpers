'''
Created on 28/9/2015

@author: victor
'''
import unittest
import trajectory_comparison.data
import numpy
import os
import matplotlib.pyplot as plt
import math
from trajectory_comparison.tools import define_tableau20_cm
from trajectory_comparison.angvar_tools import get_side_window_mean,\
    get_maximums, remove_consecutive_values, calculate_distribution_range,\
    delete_overlapping_distributions, calculate_range_population_percent,\
    circular_range, circular_subtraction

class TestAngularStats(unittest.TestCase):
    not_working_observations = [23,24,27,28,29,30,46,47]
    
    ## JUNTAR MAXIMOS QUE NO ESTEN SEPARADOS POR MUCHO (ej. 10 casillas = 1 rad)

    def test_histograms(self):
        data_folder = trajectory_comparison.data.__path__[0]
        print "data folder", data_folder
        
        all_angles = numpy.loadtxt(os.path.join(data_folder, "pro_noh_md.pdb.rot.ang"))
        
        print "all_angles.shape", all_angles.shape
        print "all_angles[0].shape", all_angles[0].shape
        
        
        for observation in range(0,100): #len(all_angles)):
            
            print "OBSERVATION", observation
            hist_bin_size = 0.1
            hist, _ = numpy.histogram(all_angles.T[observation], 
                                           bins=int((math.pi*2)/ hist_bin_size),
                                           range = (-math.pi, math.pi))
            
            ## Get_maximums with window mean
            window_size = 3
#             maximums = get_maximums(hist, window_size)
            maximums = get_maximums(hist)
            print "All maximums:", maximums
            
            ## Remove consecutive values
            maximums = remove_consecutive_values(maximums)
            print "Maximums after sequence shrinking:",maximums 
            
            # Slope descent to calculate the distributions for each maximum
            distribution_ranges = {}
            for i, maximum in enumerate(maximums):
                distribution_ranges[maximum] = calculate_distribution_range(maximum, hist, window_size)
            print "Distribution ranges:", distribution_ranges
            
            delete_overlapping_distributions(distribution_ranges, hist)
            
            print "Overlaps deleted:", distribution_ranges
            
            # Calculate population for each range
            populations = {}
            total_population = hist.sum()
            for maximum in distribution_ranges:
                populations[maximum] = calculate_range_population_percent(total_population, hist, distribution_ranges[maximum])
            
            print "Population percents:", populations
            
            # If the population is less than 20%, aggregate to the closest maximum with more than 20% population
            big_pops = {}
            for pop in populations:
                if populations[pop] >= 0.2:
                    big_pops[pop] = [pop]
            print "big pops", big_pops
            
            if len(populations) > 1:
                for pop in populations:
                    if populations[pop] < 0.2:
                        subtractions = []
                        for big_pop in big_pops:
                            subtractions.append((circular_subtraction(pop,big_pop,len(hist)) ,pop, big_pop))
                        subtractions.sort()
                        big_pops[subtractions[0][2]].append(pop)
                        print subtractions
                print big_pops
            
            plt.subplot(211)
            plt.hist(all_angles.T[observation], 
                     bins = int((math.pi*2)/ hist_bin_size), 
                     range = (-math.pi, math.pi), 
                     normed=True)
              
            plt.subplot(212)
            colors = define_tableau20_cm()
            for i, population_index in enumerate(big_pops):
                pop_hist = numpy.zeros(len(hist))
                for maximum in big_pops[population_index]:
                    max_pop_range = circular_range(distribution_ranges[maximum][0],distribution_ranges[maximum][1],len(hist))
                    for j in max_pop_range:
                        pop_hist[j] = hist[j]
                
                plt.bar(range(len(hist)), pop_hist, width=1, color= colors[i], align='center')
            plt.show()
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()