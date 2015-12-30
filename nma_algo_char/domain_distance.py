"""
Created on Dec 19, 2015

@author: victor
"""
import glob
import os
from prody.proteins.pdbfile import parsePDB
from anmichelpers.tools.tools import norm
import matplotlib.pyplot as plt
import seaborn as sns
import cPickle as pickle
from optparse import OptionParser
from nma_algo_char.common import create_directory
import numpy


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-d", dest="data")
    parser.add_option("--results", dest="results")
    
    (options, args) = parser.parse_args()
    if not options.results:  
        parser.error('Results folder not given')
    else:
        create_directory(options.results)
        
    folders = glob.glob("CC_*")
    workspaces = ["cc_closed_rmsg_dispf","cc_open_rmsg_dispf"]
    workspace_titles = {"cc_closed_rmsg_dispf":"Closed",
                        "cc_open_rmsg_dispf":"Open"}
    workspace_colors = {"cc_closed_rmsg_dispf":"g",
                        "cc_open_rmsg_dispf":"b"}
    if options.data is None:    
        distances = {}
        for workspace in workspaces:
            distances[workspace] = {}
            for folder in folders:
                T = int(folder.split("_")[1])
                inner_folders = glob.glob(os.path.join(folder,workspace,"*rmsg*"))
                distances[workspace][T] = {}
                for inner_folder in inner_folders:
                    name = os.path.basename(inner_folder)
                    prefix,p1,v1,p2,v2 = name.split("_")
                    
                    try:
                        traj_file = os.path.join(inner_folder,"trajectory.pdb")
                        pdb = parsePDB(traj_file, subset='ca')
                        res_coords = pdb.select("resid 277 or resid 387").getCoordsets()
                        ds = []
                        for cys_coords, leu_coords in res_coords:
                            ds.append(norm(leu_coords-cys_coords))
                        distances[workspace][T][v1,v2] = ds
                    except IOError:
                        pass
        pickle.dump(distances, open(os.path.join(options.results, "domain_dist_data"),"w"))
    else:
        distances = pickle.load(open(options.data))
    
#     for wp in distances:
#         for T in sorted(distances[wp].keys()):
#             for v1,v2 in distances[wp][T].keys():
#                 plt.plot(distances[wp][T][v1,v2], label = "%d"%(T))    
#         plt.title(workspace_titles[wp])
#         plt.show()
    
    for ws in distances:
        max_samples = 0
        for T in distances[ws]:
            for key in distances[ws][T]:
                max_samples = max(max_samples, len(distances[ws][T][key]))
        all_samples = []
        for T in distances[ws]:
            for key in distances[ws][T]:
                distances[ws][T][key] = numpy.array(distances[ws][T][key])
                distances[ws][T][key].resize(max_samples)
                all_samples.append(distances[ws][T][key])
        all_samples = numpy.array(all_samples)
        
        # Change 0s by avg
        for i in range(len(all_samples.T)):
            avg_nz = [] 
            for n in all_samples.T[i]:
                if n != 0: avg_nz.append(n)
            avg_nz = numpy.mean(avg_nz)
            for j,n in enumerate(all_samples.T[i]):
                if n == 0: all_samples.T[i][j] = avg_nz
        
        sns.tsplot(all_samples,ci=[68, 95, 100], c= workspace_colors[ws], err_style="unit_traces")
    plt.show()
            
    