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
from trajectory_comparison.angvar_tools import get_maximums, remove_consecutive_values, calculate_distribution_range,\
    delete_overlapping_distributions, calculate_range_population_percent,\
    merge_small_with_big,\
    merge_big_with_big_by_closeness, get_ordered_distributions,\
    get_angular_ranges, get_distribution_angles, shift_angular_distribution_to_0_2pi

class TestAngularStats(unittest.TestCase):
    
    ## JUNTAR MAXIMOS QUE NO ESTEN SEPARADOS POR MUCHO (ej. 10 casillas = 1 rad)

    def test_histograms(self):
        data_folder = trajectory_comparison.data.__path__[0]
        all_angles = numpy.loadtxt(os.path.join(data_folder, "pro_noh_md.pdb.rot.ang"))
        
        [32] # obs partida
        for observation in [32]: #range(3,100): #len(all_angles)):
            
            print "OBSERVATION", observation
            hist_bin_size = 0.1
            hist, _ = numpy.histogram(all_angles.T[observation], 
                                           bins=int((math.pi*2)/ hist_bin_size),
                                           range = (-math.pi, math.pi))
            
            ## Get_maximums with window mean
            maximums = get_maximums(hist)
            print "All maximums:", maximums
            
            ## Remove consecutive values
            maximums = remove_consecutive_values(maximums)
            print "Maximums after sequence shrinking:",maximums 
            
            # Slope descent to calculate the distributions for each maximum
            JUMP_SIZE = 15
            distribution_ranges = {}
            for i, maximum in enumerate(maximums):
                distribution_ranges[maximum] = calculate_distribution_range(maximum, hist, JUMP_SIZE)
            print "Distribution ranges:", distribution_ranges
            
            delete_overlapping_distributions(distribution_ranges, hist)
            
            print "Overlaps have been deleted:", distribution_ranges
            
            # Calculate population for each distribution
            distribution_population = {}
            total_population = hist.sum()
            for maximum in distribution_ranges:
                distribution_population[maximum] = calculate_range_population_percent(total_population, hist, 
                                                                                      distribution_ranges[maximum])
            print "Distribution population percents:", distribution_population
            
            # If the population is less than 20%, aggregate to the closest maximum with more than 20% population
            big_distributions = {}
            for maximum in distribution_population:
                if distribution_population[maximum] >= 0.2:
                    big_distributions[maximum] = [maximum]
            print "'Big' distributions:", big_distributions
            
            # Merge small distributions with the big ones
            merge_small_with_big(distribution_population,big_distributions, hist)
            print "'Big' distributions with merged small:", big_distributions
            
            
            # Merge big "close" distributions along with their "small" distributions
            CLOSENESS_THRESHOLD = 10 # for a 0.1 bin means 1rad
            merged_big_distributions = merge_big_with_big_by_closeness(big_distributions, distribution_population, CLOSENESS_THRESHOLD)
            print "Merged big distributions", merged_big_distributions


            # plot!                        
            plt.subplot(211)
            plt.hist(all_angles.T[observation], 
                     bins = int((math.pi*2)/ hist_bin_size), 
                     range = (-math.pi, math.pi), 
                     normed=True)
            
            plt.subplot(212)
            distributions = get_ordered_distributions(merged_big_distributions, distribution_ranges, hist)
            colors = define_tableau20_cm()
            for i, dist in enumerate(distributions):
                plt.bar(range(len(hist)), dist, width=1, color = colors[i], align='center', alpha = 0.5)
            plt.show()
            plt.clf()
            
            for distribution in distributions:
                print distribution
                index_range, angle_range = get_angular_ranges(distribution, -math.pi, hist_bin_size)
                print index_range, angle_range
                angles = get_distribution_angles(all_angles.T[observation], angle_range)
                shifted_angles, shift = shift_angular_distribution_to_0_2pi(angles, distribution, index_range, hist_bin_size)
                plt.hist(shifted_angles, 
                     bins = int((math.pi*3)/ hist_bin_size), 
                     range = (-math.pi, 2*math.pi), 
                     normed=True)
                
                plt.show()
                plt.clf()
            
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()