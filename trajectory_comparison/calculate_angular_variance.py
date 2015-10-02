'''
Created on 25/9/2015

@author: victor
'''
import sys
import numpy
import math
import trajectory_comparison.angvar_tools as ang_tools
from pca.dpca import circular_mean_and_variance, circular_variance


def calculate_distributions(observation, hist_bin_size = 0.1, jump_size = 15,  closeness_threshold = 10):
    
    hist, _ = numpy.histogram(observation, 
                                   bins=int((math.pi*2)/ hist_bin_size),
                                   range = (-math.pi, math.pi))
    
    ## Get_maximums with window mean
    maximums = ang_tools.get_maximums(hist)
    
    ## Remove consecutive values
    maximums = ang_tools.remove_consecutive_values(maximums)
    
    # Slope descent to calculate the distributions for each maximum
    distribution_ranges = {}
    for maximum in maximums:
        distribution_ranges[maximum] = ang_tools.calculate_distribution_range(maximum, hist, jump_size)
    
    ang_tools.delete_overlapping_distributions(distribution_ranges, hist)
    
    
    # Calculate population for each distribution
    distribution_population = {}
    total_population = hist.sum()
    for maximum in distribution_ranges:
        distribution_population[maximum] = ang_tools.calculate_range_population_percent(total_population, hist, 
                                                                              distribution_ranges[maximum])
    
    # If the population is less than 20%, aggregate to the closest maximum with more than 20% population
    big_distributions = {}
    for maximum in distribution_population:
        if distribution_population[maximum] >= 0.2:
            big_distributions[maximum] = [maximum]
    
    # Merge small distributions with the big ones
    ang_tools.merge_small_with_big(distribution_population,big_distributions, hist)
    
    
    # Merge big "close" distributions along with their "small" distributions
    merged_big_distributions = ang_tools.merge_big_with_big_by_closeness(big_distributions, distribution_population, closeness_threshold)

    # "Mount" the real distributions and order them by population (bigger go first)
    distributions = ang_tools.get_ordered_distributions(merged_big_distributions, distribution_ranges, hist)
    
    return distributions

if __name__ == '__main__':
    
    all_angles = numpy.loadtxt(sys.argv[1])
    
    HIST_BIN_SIZE = 0.1
    JUMP_SIZE = 15
    CLOSENESS_THRESHOLD = 10
    
    all_variances = []
    max_num_cols = 0
    for i, torsion_observations in enumerate(all_angles.T[3:]):
        try:
            variances = []
            distributions = calculate_distributions(torsion_observations,HIST_BIN_SIZE,JUMP_SIZE,CLOSENESS_THRESHOLD)
            for distribution in distributions:
                index_range, angle_range = ang_tools.get_angular_ranges(distribution, -math.pi, HIST_BIN_SIZE)
                angles = ang_tools.get_distribution_angles(torsion_observations, angle_range)
                shifted_angles, shift = ang_tools.shift_angular_distribution_to_0_2pi(angles, distribution, index_range, HIST_BIN_SIZE)
                ang_variance = circular_variance(shifted_angles)
                variances.append(ang_variance)
            all_variances.append(variances)
            max_num_cols = max(len(variances),max_num_cols)
        except:
            all_variances.append([0])

    padded_variances = numpy.zeros((len(all_angles.T), max_num_cols+1))
    for i, dist_vars in enumerate(all_variances):
        padded_variances[i][0] = i
        for j, var in enumerate(dist_vars):
            padded_variances[i][j+1] = var
    numpy.savetxt("%s.var"%sys.argv[1], padded_variances, fmt="%.5f")
    
    