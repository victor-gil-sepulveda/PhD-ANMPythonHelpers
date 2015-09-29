import numpy
import math

def get_centered_window_mean(array, position, win_side):
    """
    [***x***] here position is x, and win_side = 3 (includes x)
    """
    tmp = []
    for i in range(position-win_side, position +win_side+1):
        tmp.append(array[i % len(array)])
    return numpy.mean(tmp)

def get_side_window_mean(array, position, direction, win_size):
    """
    [***x] here position is x, direction is -1 and win_size = 3.
    Includes x
    """
    tmp = [array[position]]
    for i in range(position+direction, position + direction + direction*win_size, direction):
        tmp.append(array[i % len(array)]) 
    return numpy.mean(tmp)

def get_maximums(hist):
    ## get_maximums
    maximums = []
    for i in range(len(hist)):
        i_previous = (i-1) % len(hist)
        i_next = (i+1) % len(hist)
        if hist[i] > 0 and hist[i_previous] <= hist[i] and hist[i_next] <= hist[i]:
            maximums.append(i)
    return maximums

def remove_consecutive_values(maximums):
    """
    Removes the consecutive series in an array(leaving only the first in the series) and returns it.
    Ex. in 1,3,4,5,7,8,9 it would remove 4,5,8,9 
    """
    array_len = len(maximums)
    to_delete = []
    for i in range(array_len):
        next_i = (i+1) % array_len
        if maximums[i]+1 == maximums[next_i]:   
            to_delete.append(True)
        else:
            to_delete.append(False)
    #print to_delete
    no_reps = []
    for i in range(array_len):
        if not to_delete[i]:
            no_reps.append(maximums[i])
    return no_reps


def go_down_the_hill(starting_pos, direction, hist, allowed_jump):
    i = starting_pos
    next_i = (i+direction) % len(hist)
    jumped_once = False
    while hist[next_i] > 0 and hist[next_i]-hist[i] <= allowed_jump:
        if hist[next_i]-hist[i] > 0:
            if jumped_once:
                i = i-direction # we do not want to allow two jumps, it may mean we are going uphill
                                # we have to undo last step
                break
            else:
                jumped_once = True
        else:
            jumped_once = False

        i = next_i
        next_i = (i+direction) % len(hist)
    return i

def calculate_distribution_range(maximum_index, hist, allowed_jump):
    """
    Calculates where the angular distribution starts and ends.
    Can only perform ONE jump.
    """
    # Go down the hill until values decrease
    
    left_limit = go_down_the_hill(maximum_index, -1, hist, allowed_jump)
    right_limit = go_down_the_hill(maximum_index, 1, hist, allowed_jump)
   
    return  (left_limit, right_limit)

def circular_range(ini, end, length):
    """
    Works as the "range" function, but in a circular environment.
    Ex. range 5, 2 given an array length of 7 would return 5,6,0,1,2
    """
    if ini > end:
        return range(ini, length)+range(0, end+1)
    else:
        return range(ini, end+1)

def delete_overlapping_distributions(distribution_ranges, hist):
    # Merge those that completely overlap (delete the smaller one)
    delete_minor_overlaps = []
    maximums = distribution_ranges.keys()
    for i, i_maximum in enumerate(maximums):
        i_range = set(circular_range(distribution_ranges[i_maximum][0],distribution_ranges[i_maximum][1],len(hist)))
        for j_maximum in maximums[i+1:]:
            j_range = set(circular_range(distribution_ranges[j_maximum][0],distribution_ranges[j_maximum][1],len(hist)))
            if i_range.issubset(j_range):
                delete_minor_overlaps.append(i_maximum)
            elif j_range.issubset(i_range):
                delete_minor_overlaps.append(j_maximum)
    # Delete overlaps
    for to_del in set(delete_minor_overlaps):
        del distribution_ranges[to_del]
        
def calculate_range_population_percent(total_population, hist, p_range):
    max_range = circular_range(p_range[0], p_range[1],len(hist))
    population = 0.
    for i in max_range:
        population += hist[i]
    return float(population)/total_population

def circular_subtraction(a,b,length):
    # shortest path substraction
    return min(abs(a-b), abs(a+length-b))

def merge_small_with_big(distribution_population,big_distributions, hist):
    """
    Inplace merging
    """
    for maximum in distribution_population:
        if distribution_population[maximum] < 0.2:
            subtractions = []
            for big_distribution_maximum in big_distributions:
                subtractions.append((circular_subtraction(maximum, big_distribution_maximum, 
                                                          len(hist)) ,
                                     maximum, big_distribution_maximum))
            subtractions.sort()
            big_distributions[subtractions[0][2]].append(maximum)

def merge_big_with_big_by_closeness(big_distributions, distribution_population, CLOSENESS_THRESHOLD):
    to_merge = []
    # sort maximums by REAL normalized total population
    new_dist_populations = []
    for maximum in big_distributions:
        new_population = 0
        for other_max in big_distributions[maximum]:
            new_population += distribution_population[other_max]
        new_dist_populations.append((new_population,maximum))
    new_dist_populations.sort(reverse=True)
    
    for i, (_, i_maximum) in enumerate(new_dist_populations):
        for _, j_maximum in new_dist_populations[i+1:]:
            if abs(i_maximum-j_maximum) < CLOSENESS_THRESHOLD: 
                to_merge.append([abs(i_maximum-j_maximum), i_maximum, j_maximum, True])
    
    to_merge.sort() # small distances first
    print "Merging units:", to_merge
    
    merged_big_distributions = {}
    mergeable = {}
    for maximum in big_distributions:
        mergeable[maximum] = True
        
    for _, maximum in new_dist_populations:
        if mergeable[maximum]:
            merged_big_distributions[maximum] = big_distributions[maximum]
            print "Adding", maximum
            for _, i, j, active in to_merge:
                if active and maximum == i: # there's a smaller one to merge
                    merged_big_distributions[maximum].extend(big_distributions[j])
                    print "Merging ", i , "with", j
                    # Change j by i if it is in first position, remove the merging if it is 
                    # in second position
                    for merging_unit in to_merge:
                        if merging_unit[1] == j:
                            merging_unit[1] = i
                            
                        if merging_unit[2] == j:
                            merging_unit[3] = False
                    mergeable[j] = False
    return merged_big_distributions

def get_ordered_distributions(merged_big_distributions, distribution_ranges, hist):
    tmp_distributions = []
    for maximum in merged_big_distributions:
        base_dist_hist = numpy.zeros(len(hist))
        for i in merged_big_distributions[maximum]:
            dist_range = circular_range(distribution_ranges[i][0],distribution_ranges[i][1],len(hist))
            for j in dist_range:
                base_dist_hist[j] = hist[j]
        tmp_distributions.append((base_dist_hist.sum(),base_dist_hist))
    tmp_distributions.sort(reverse = True)
    
    distributions = []
    for i in range(len(tmp_distributions)):
        distributions.append(tmp_distributions[i][1])
    return distributions

def find_first_zero(distribution, starting_position, direction):
    index = starting_position
    while distribution[index] != 0:
        index = (index+direction) % len(distribution)
    return index

def find_first_non_zero(distribution, starting_position, direction):
    index = starting_position
    while distribution[index] == 0:
        index = (index+direction) % len(distribution)
    return index
    
def to_0_2pi_range(angles):
    """
    Each angle is in -pi, pi
    """
    new_angles = []
    for angle in angles:
        if angle < 0:
            new_angles.append(angle+2*math.pi)
        else:
            new_angles.append(angle)

def get_angular_ranges(distribution, min_val, bin_size):
    
    if distribution[0] != 0:
        print "*"
        # If starts with non-zero, may be a dist that crosses pi
        # search from the middle (0rad)
        mid = len(distribution) / 2
        print "DIST MID", distribution[mid]
        right_index =  find_first_non_zero(distribution, 
                                           find_first_zero(distribution,mid,1),
                                           -1)
        left_index = find_first_non_zero(distribution, 
                                         right_index+1, 
                                         1)
        return (left_index, right_index), (min_val+(left_index*bin_size), min_val+(right_index*bin_size))
    else:
        print "**"
        left_index = find_first_non_zero(distribution, 0, 1)
        right_index =  find_first_non_zero(distribution, 0, -1)
        return (left_index,right_index), (min_val+(left_index*bin_size), min_val+(right_index*bin_size))
    
def get_distribution_angles(full_observation, angle_range):
    distribution_angles = []
    if  angle_range[1] < angle_range[0]:#crosses_pi:
        for angle in full_observation:
            if angle <= angle_range[1] or angle >= angle_range[0]:
                distribution_angles.append(angle)   
    else:
        for angle in full_observation:
            if angle >= angle_range[0] and angle <= angle_range[1]:
                distribution_angles.append(angle)
    return distribution_angles

def shift_angular_distribution_to_0_2pi(distribution_angles, distribution, dist_limits, bin_size):
    #maximum_index = distribution.index(max(distribution))
    shift = -math.pi + float(dist_limits[0])*bin_size
    print "SHIFT", shift
    shifted_angles = []
    for angle in distribution_angles:
        if angle >= 0:
            shifted_angles.append(angle - shift)
        else:
            shifted_angles.append(angle - shift + 2*math.pi)
    return shifted_angles, shift
    