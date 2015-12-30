'''
Created on Dec 29, 2015

@author: victor
'''
import cPickle as pickle
import os
import seaborn as sns
from optparse import OptionParser
import numpy
import matplotlib.pyplot as plt

ENERGY_LABEL = "$\Delta$ U"
RMSD_LABEL = "RMSD"

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("--wp1", dest="wp1")
    parser.add_option("--wp2", dest="wp2")
    parser.add_option("--results", dest="results")
    
    (options, args) = parser.parse_args()
    
    workspaces = [options.wp1, options.wp2]
    
    sns.set_style("whitegrid")
    colors = sns.color_palette("hls", 2)
    for c, wp in enumerate(workspaces):
        l = os.path.basename(os.path.normpath(wp))
        all_data = pickle.load(open(os.path.join(wp,"all_data.pickle")))
        all_data["cost"] = numpy.array(all_data["time_per_step"]*numpy.array(all_data[ENERGY_LABEL]))
        ps = plt.scatter(all_data[RMSD_LABEL], all_data["cost"], c = colors[c],
                    label = l, alpha=0.7)
    plt.legend()
    plt.xlabel("RMSD")
    plt.ylabel("Cost")
    plt.show()
#     plt.savefig("test.svg")
    
        
        
        