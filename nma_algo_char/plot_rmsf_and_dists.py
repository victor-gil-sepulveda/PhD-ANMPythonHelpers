'''
Created on Feb 28, 2016

@author: victor
'''
from optparse import OptionParser
from nma_algo_char.common import create_directory
import numpy
import matplotlib.pyplot as plt
import seaborn as sns
import os.path
from trajectory_comparison.compare_two_rmsfs import rms

sns.set_style("whitegrid")

def rmsf_rescale(reference, rmsf):
    scaled = None
    results = []
    for scale in numpy.arange(1, 5, 0.1):
        scaled = rmsf * scale
        results.append((rms(reference, scaled), scale))
#     data = numpy.array(results)
#     plt.plot(data.T[1], data.T[0])
#     plt.show()
#     plt.close()
    scale = min(results)[1]
    return rmsf*scale, scale

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--results", dest="results_folder")
    parser.add_option("-d", default = False, action="store_true", dest="plot_distances")
    
    (options, args) = parser.parse_args()
    
    files = ["md","cc","ic_m1"]
    file_legend = {"md":"MD","cc":"CC","ic_m1": "IC (m1)"}
    
    create_directory(options.results_folder)
    
    # RMSF plot
    for filename in files:
        rmsf = numpy.loadtxt("%s.rmsf"%filename)
        if len(rmsf.shape)>1:
            rmsf = rmsf.T[1]
        plt.plot(rmsf, label = file_legend[filename])
    plt.legend()
    plt.ylabel("Root mean square fluctuation ($\AA$)")
    plt.xlabel("Residue number")
#     plt.xticks(range(1,len(rmsf),20), [str(s) for s in range(259,533+1,20)], rotation=45)
    plt.xticks(range(1,len(rmsf),10), [str(s) for s in range(1,len(rmsf),10)], rotation=45)
    plt.savefig(os.path.join(options.results_folder,"rmsf.svg"))
    plt.close()
    
    # RMSF rescaling without first residues
#     reference = numpy.loadtxt("md.rmsf").T[1][6:-6]
    reference = numpy.loadtxt("md.rmsf").T[1][:-6]
    plt.plot(reference, label = file_legend["md"])
    for filename in ["cc","ic_m1"]:
        rmsf = numpy.loadtxt("%s.rmsf"%filename)
        print len(rmsf), len(numpy.loadtxt("md.rmsf").T[1])
        if len(rmsf.shape)>1:
            rmsf = rmsf.T[1]
#         rescaled, scale = rmsf_rescale(reference, rmsf[6:-5])
        rescaled, scale = rmsf_rescale(reference, rmsf[:-3])
        plt.plot(rescaled, label = "%s (x %.2f)"%(file_legend[filename], scale))
#     plt.xticks(range(1,len(rmsf),20), [str(s) for s in range(259+6,533+1,20)], rotation=45)
    plt.xticks(range(1,len(rmsf),10), [str(s) for s in range(1,len(rmsf), 10)], rotation=45)
    plt.ylabel("Root mean square fluctuation ($\AA$)")
    plt.xlabel("Residue number")
    plt.legend()
    plt.savefig(os.path.join(options.results_folder,"rescaled_rmsf.svg"))
    plt.close()
    
    if options.plot_distances:
        min_dist = numpy.inf
        max_dist = 0
        for filename in files:
            dists = numpy.loadtxt("%s.dists"%filename)
            if len(dists.shape)>1:
                dists = dists.T[1]
            min_dist = min(min_dist, min(dists))
            max_dist = max(max_dist, max(dists))
            
        for filename in files:
            dists = numpy.loadtxt("%s.dists"%filename)
            if len(dists.shape)>1:
                dists = dists.T[1]
            plt.hist(dists, bins = 100, range =(min_dist, max_dist), normed=True, alpha = 0.7, label = file_legend[filename])
        plt.xlabel("Distance")
        plt.legend()
        plt.savefig(os.path.join(options.results_folder,"domain_distances.svg"))
    
    
    