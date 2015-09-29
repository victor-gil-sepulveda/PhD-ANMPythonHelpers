import numpy

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
#         print i,
#     print "window side", position, tmp
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

# def get_maximums(hist, window_size):
#     maximums = []
#     for i in range(len(hist)):
#         
#         left_side_mean = get_side_window_mean(hist, i, -1, window_size)
#         right_side_mean = get_side_window_mean(hist, i, 1, window_size)
# #         left_side_mean = get_centered_window_mean(hist, i-1,  window_size)
# #         right_side_mean = get_centered_window_mean(hist, i+1, window_size)
# #         if i == 34:
# #             print hist[i], "left", left_side_mean, "right",  right_side_mean 
# #             exit()
#         if hist[i] > 0 and  left_side_mean <= hist[i] and  right_side_mean <= hist[i]:
#             maximums.append(i)
#     return maximums

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


def calculate_distribution_range(maximum_index, hist, window_size):
    """
    Calculates where the angular distribution starts and ends.
    """
    # Go down the hill until values increase
    i = maximum_index
#     while get_side_window_mean(hist, i, -1, window_size) <= hist[i]:
    next_i = (i-1) % len(hist)
    while get_side_window_mean(hist, i, -1, window_size) <= hist[i]:
#     while get_centered_window_mean(hist, i-1, window_size) <= hist[i]:
#         if maximum_index == 31:
#             print "left slope",i, hist[i], get_side_window_mean(hist, i, 1, window_size)
        i = next_i
        next_i = (i-1) % len(hist)
    left_limit = i
#     if maximum_index == 31:
#         print "left slope",i, hist[i], get_side_window_mean(hist, i, 1, window_size)
       
    
    i = maximum_index
#     while get_side_window_mean(hist, i, 1, window_size) <= hist[i]:
    next_i = (i+1) % len(hist)
    while hist[next_i] <= hist[i]:
#         if maximum_index == 31:
#             print "right slope",i,hist[i], get_side_window_mean(hist, i, 1, window_size)
#     while get_centered_window_mean(hist, i+1, window_size) <= hist[i]:
        i = next_i
        next_i = (i+1) % len(hist)
#     if maximum_index == 31:
#         print "right slope",i,hist[i], get_side_window_mean(hist, i, 1, window_size)
    right_limit = i
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
        print i_maximum, i_range
        for j_maximum in maximums[i+1:]:
            j_range = set(circular_range(distribution_ranges[j_maximum][0],distribution_ranges[j_maximum][1],len(hist)))
            if i_range.issubset(j_range):
                delete_minor_overlaps.append(i_maximum)
#                 print i_range, "overlaps with", j_range
            elif j_range.issubset(i_range):
                delete_minor_overlaps.append(j_maximum)
#                 print j_range, "overlaps with", i_range
    # Delete overlaps
    for to_del in set(delete_minor_overlaps):
        del distribution_ranges[to_del]
        
def calculate_range_population_percent(total_population, hist, p_range):
    max_range = circular_range(p_range[0], p_range[1],len(hist))
    population = 0
    for i in max_range:
        population += hist[i]
    return float(population)/total_population

def circular_subtraction(a,b,length):
    # shortest path substraction
    return min(abs(a-b), abs(a+length-b))