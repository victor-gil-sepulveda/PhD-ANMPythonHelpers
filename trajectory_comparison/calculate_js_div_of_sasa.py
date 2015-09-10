"""
This script does two things. First compares the SASA distributions of a given set of trajectories, then 
it compares that distributions with a reference distribution.
Some parts of the code can be uncommented to get, if needed, the plot of the distributions in the 
second case (can be a way of detecting if they overlap visually).
Input files can be manually changed in order to analyze another distributions of values. 
"""

#----------INPUT FILES ----------------

# input_files = {
# 
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/1432_0_25.sasa":{"T":1432, "Ca_dist":0.25},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/1432_0_66.sasa":{"T":1432, "Ca_dist":0.66},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/1432_1_08.sasa":{"T":1432, "Ca_dist":1.08},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/1432_1_5.sasa":{"T":1432, "Ca_dist":1.5},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/2000_0_25.sasa":{"T":2000, "Ca_dist":0.25},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/2000_0_66.sasa":{"T":2000, "Ca_dist":0.66},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/2000_1_08.sasa":{"T":2000, "Ca_dist":1.08},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/2000_1_5.sasa":{"T":2000, "Ca_dist":1.5},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/300_0_25.sasa":{"T":300, "Ca_dist":0.25},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/300_0_66.sasa":{"T":300, "Ca_dist":0.66},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/300_1_08.sasa":{"T":300, "Ca_dist":1.08},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/300_1_5.sasa":{"T":300, "Ca_dist":1.5},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/866_0_25.sasa":{"T":866, "Ca_dist":0.25},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/866_0_66.sasa":{"T":866, "Ca_dist":0.66},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/866_1_08.sasa":{"T":866, "Ca_dist":1.08},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/866_1_5.sasa":{"T":866, "Ca_dist":1.5}
# }
# 
# ordered_input_files = [
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/300_0_25.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/300_0_66.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/300_1_08.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/300_1_5.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/866_0_25.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/866_0_66.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/866_1_08.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/866_1_5.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/1432_0_25.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/1432_0_66.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/1432_1_08.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/1432_1_5.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/2000_0_25.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/2000_0_66.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/2000_1_08.sasa",
#     "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/results/2000_1_5.sasa",
# ]
# 
# # Ordered labels
# Ts = [300, 866, 1432, 2000]
# Cas = [0.25, 0.66, 1.08, 1.5]

#--------------------------------------
# 
# input_files = {
# # "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/TestANMCCBest10modes/all.pdb.sasa":{"T":2000, "Ca_dist":1.5},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_sdmin/all.pdb.sasa":{"T":300, "Ca_dist":2.},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_500_sdmin/all.pdb.sasa":{"T":500, "Ca_dist":2.},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_700_sdmin/all.pdb.sasa":{"T":700, "Ca_dist":2.},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Shaw/pro_noh_md.pdb.residue_8_to_266.sasa":{"T":"MD", "Ca_dist":2.}
# 
# }
# 
# ordered_input_files = [
# # "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/TestANMCCBest10modes/all.pdb.sasa",
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_sdmin/all.pdb.sasa",
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_500_sdmin/all.pdb.sasa",
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_700_sdmin/all.pdb.sasa"
# 
# ]
# 
# # Ordered labels
# Ts = [300, 500, 700]
# Cas = [2.,2.,2.]
#--------------------------------------

input_files = {
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/TestANMCCBest10modes/all.pdb.residue_8_to_266.sasa":{"T":2000, "Ca_dist":1.5},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_sdmin/all.pdb.residue_8_to_266.sasa":{"T":300, "Ca_dist":0.2},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_nosdmin/all.pdb.residue_8_to_266.sasa":{"T":"300_no", "Ca_dist":0.2},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_sdmin_1.0/all.pdb.residue_8_to_266.sasa":{"T":300, "Ca_dist":0.1},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_500_sdmin/all.pdb.residue_8_to_266.sasa":{"T":500, "Ca_dist":0.2},
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_700_sdmin/all.pdb.residue_8_to_266.sasa":{"T":700, "Ca_dist":0.2},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Shaw/pro_noh_md.pdb.residue_8_to_266.sasa":{"T":"MDShaw", "Ca_dist":""},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Suwipa/mdgb1.pdb.residue_8_to_266.sasa":{"T":"MDSuw", "Ca_dist":"igb1"}
}

ordered_input_files = [
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/TestANMCCBest10modes/all.pdb.residue_8_to_266.sasa",
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_sdmin/all.pdb.residue_8_to_266.sasa",
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_nosdmin/all.pdb.residue_8_to_266.sasa",
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_sdmin_1.0/all.pdb.residue_8_to_266.sasa",
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_500_sdmin/all.pdb.residue_8_to_266.sasa",
# "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_700_sdmin/all.pdb.sasa",
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Shaw/pro_noh_md.pdb.residue_8_to_266.sasa",
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Suwipa/mdgb1.pdb.residue_8_to_266.sasa"
]

# Ordered labels
Ts = [2000, 300, "300_no", 300, 500, "MDSuw", "MDShaw"]
Cas = [1.5, 0.2, 0.2, 0.1, 0.2, "igb1", "",]
#--------------------------------------

import numpy 
import math
from scipy.stats import entropy
from numpy.linalg import norm
import matplotlib.pyplot as plt

def JSD(P, Q):
    """
    Calculates the Jensen-Shannon divergence as a metric (sq_root)
    See: http://www.researchgate.net/publication/3084774_A_new_metric_for_probability_distributions
    """
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return math.sqrt(0.5 * (entropy(_P, _M) + entropy(_Q, _M)))

def smoothed(distribution,small_value = 1.0e-8):
    """
    Applies a smoothing process to the distribution.
    See http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
    for an explanation about the problem and the solution.
     
    @param distribution: distribution to be smoothed
    @param small_value: value to be set to those bins with 0 probability
     
    @return: The smoothed distribution.
    """
    total_number_of_samples = len(distribution)
    samples_in_distrib = numpy.count_nonzero(distribution)
    pc = small_value * (total_number_of_samples - samples_in_distrib) / samples_in_distrib
    smoothed_distrib = numpy.empty(len(distribution))
    for i in range(len(distribution)):
        if distribution[i] == 0:
            smoothed_distrib[i] = small_value
        else:
            smoothed_distrib[i] = distribution[i] - pc
    return numpy.array(smoothed_distrib)

def show_table(data, x_label, y_label):
    """
    Shows a "heatmap"-like plot, adds the labels and the color bar legend.
    """
    _, ax = plt.subplots()
    plt.subplots_adjust(top = 0.85, left = 0.135)
    heatmap = ax.pcolormesh(data, cmap=plt.cm.Blues)
    # put the major ticks at the middle of each cell
    ax.set_xticks(numpy.arange(data.shape[0])+1, minor=False)
    ax.set_yticks(numpy.arange(data.shape[1])+0.5, minor=False)
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.set_xticklabels(x_label, minor=False, rotation=45)
    ax.set_yticklabels(y_label, minor=False)
    plt.colorbar(heatmap)
    plt.show()

def load_values(files):
    """
    Loads the sasa values from the files in the files dictionary. Returns the maximum and minimum values.
    Values are loaded into the "files" structure.
    """
    min_val = float("inf")
    max_val =  0.0
    for filename in input_files:
        input_files[filename]["values"] = numpy.loadtxt(filename)
        max_val = max(numpy.max(input_files[filename]["values"]), max_val)
        min_val = min(numpy.min(input_files[filename]["values"]), min_val) 

    return max_val,  min_val

def calculate_distributions(input_files, max_val,  min_val):
    """
    Calculates the histogram (distribution) for the values inside the "files" structure.
    """
    for filename in input_files:
        input_files[filename]["distrib"] = smoothed(numpy.histogram(input_files[filename]["values"], 
                                                            bins = 100, range = (min_val, max_val), 
                                                            normed=True)[0])
#         print filename, input_files[filename]["T"], input_files[filename]["Ca_dist"]
#         plt.hist(input_files[filename]["values"], bins = 100, range = (min_val, max_val))
#         plt.show()

def test_tables_correctly_plotted():
    """
    Test snippet to visually prove that tables are correctly plotted.
    """
    table = [ [1,2,3],
              [4,5,6],
              [7,8,9]]
              
    x_labels = ["a","b","c"]
    y_labels = ["d","e","f"]

    show_table(numpy.array(table), x_labels, y_labels)

#------------------------------
# All CC vs All CC 
#------------------------------
max_val,  min_val = load_values(input_files)

calculate_distributions(input_files, max_val,  min_val )

# generate table labels
labels = []
for filename in ordered_input_files:
    labels.append(str((input_files[filename]["T"],input_files[filename]["Ca_dist"])))

# calculate table
table = []
for i in range(len(ordered_input_files)):
    row = []
    for j in range(len(ordered_input_files)):
        jsd = JSD(input_files[ordered_input_files[i]]["distrib"], 
                           input_files[ordered_input_files[j]]["distrib"])
        row.append(jsd)
    table.append(row)
    
data = numpy.array(table)
show_table(data, labels, labels)


#------------------------------
# All CC vs Reference (MD)
#------------------------------
<<<<<<< HEAD
reference_filename = "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Shaw/pro_noh_md.pdb.residue_8_to_266.sasa"
input_files[reference_filename] = {"T":0,"Ca_dist":0.}
=======
reference_filename = "/home/victor/Desktop/prots/pro_noh_md.pdb.sasa"
input_files[reference_filename] = {"T":5000,"Ca_dist":1}
>>>>>>> branch 'master' of https://github.com/victor-gil-sepulveda/PhD-ANMPythonHelpers.git

max_val,  min_val = load_values(input_files)

calculate_distributions(input_files, max_val,  min_val )
<<<<<<< HEAD
=======

# Ordered labels
Ts = [3000, 5000]
Cas = [1]
>>>>>>> branch 'master' of https://github.com/victor-gil-sepulveda/PhD-ANMPythonHelpers.git

# Reorder files indexing by T and Ca dist
files_per_CaT = {}
for filename in ordered_input_files:
    T = input_files[filename]["T"]
    Ca = input_files[filename]["Ca_dist"]
    files_per_CaT[(T,Ca)] = filename
    print (T,Ca)

<<<<<<< HEAD
# # Calculate table
# table = []
# for T in Ts:
#     row = []
#     for ca_dist in Cas:
#         jsd = JSD(input_files[files_per_CaT[(T,ca_dist)]]["distrib"], 
#                            input_files[reference_filename ]["distrib"])
#         plt.hist(input_files[files_per_CaT[(T,ca_dist)]]["values"], bins = 100, range = (min_val, max_val), normed=True)
#         plt.hist(input_files[reference_filename]["values"], bins = 100, range = (min_val, max_val), normed=True)
#           
#         plt.show()
#         row.append(jsd)
#     table.append(row)
# 
# data = numpy.array(table)
# show_table(data, Cas, Ts)
=======
# Calculate table
table = []
for T in Ts:
    row = []
    for ca_dist in Cas:
        jsd = JSD(input_files[files_per_CaT[(T,ca_dist)]]["distrib"], 
                           input_files[reference_filename ]["distrib"])
        plt.hist(input_files[files_per_CaT[(T,ca_dist)]]["values"], bins = 100, range = (min_val, max_val), normed=True)
        plt.hist(input_files[reference_filename]["values"], bins = 100, range = (min_val, max_val), normed=True)
         
        plt.show()
        row.append(jsd)
    table.append(row)
>>>>>>> branch 'master' of https://github.com/victor-gil-sepulveda/PhD-ANMPythonHelpers.git

# Calculate table
table = []
for i in range(len(Cas)):
    ca_dist  = Cas[i]
    T = Ts[i]
    print files_per_CaT[(T,ca_dist)]
    jsd = JSD(input_files[files_per_CaT[(T,ca_dist)]]["distrib"], 
                       input_files[reference_filename ]["distrib"])
    print ca_dist, T, jsd
    
    plt.hist(input_files[files_per_CaT[(T,ca_dist)]]["values"], bins = 100, range = (min_val, max_val), normed=True)
    plt.hist(input_files[reference_filename]["values"], bins = 100, range = (min_val, max_val), normed=True)
      
    plt.show()


