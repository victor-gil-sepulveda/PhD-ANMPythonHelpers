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
import matplotlib.patches as mpatches

def prepare_subplots(row_len, col_len):
        if row_len > 1 or col_len > 1:
            f, axes = plt.subplots( col_len, row_len, sharey='row', sharex='col')
            f.subplots_adjust(hspace=0.4, wspace=0.3 )
            f.set_size_inches(12, 12, forward=True)
        else:
            f = plt.gcf()
            axes = {(0,0): plt.gca()}
        return f, axes
    
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
    max_confs = 1500
    
#     folders = glob.glob("IC_*")
#     workspaces = ["ic_open_rmsg_dispf"]
#     workspace_titles = {"ic_open_rmsg_dispf":"Open"}
#     workspace_colors = {"ic_open_rmsg_dispf":"b"}
#     max_confs = 300
    
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
                        if max_confs is not None:
                            res_coords = res_coords[0:max_confs]
                            
                        ds = []
                        for cys_coords, leu_coords in res_coords:
                            ds.append(norm(leu_coords-cys_coords))
                        distances[workspace][T][v1,v2] = ds
                    except IOError:
                        print "Trajectory not read", traj_file
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
    
    temperatures = sorted(distances[distances.keys()[0]].keys())
    
    collapsed_values = {}
    for ws in workspaces:
        collapsed_values[ws] = {}
        
        # Substitute missing values with averages
        for T in temperatures:
            max_samples = 0
            for key in distances[ws][T]:
                max_samples = max(max_samples, len(distances[ws][T][key]))

            for key in distances[ws][T]:
                distances[ws][T][key] = numpy.array(distances[ws][T][key])
                distances[ws][T][key].resize(max_samples)
            
            for i in range(max_samples):
                samples = []
                for key in distances[ws][T]:
                    value = distances[ws][T][key][i]
                    if value != 0: samples.append(value)
                avg_value = numpy.mean(samples)
                for key in distances[ws][T]:
                    if distances[ws][T][key][i] == 0: distances[ws][T][key][i] = avg_value
            
            collapsed_values[ws][T] = []
            for key in distances[ws][T]:
                collapsed_values[ws][T].append(distances[ws][T][key])
            collapsed_values[ws][T] = numpy.array(collapsed_values[ws][T])
    
    sns.set_style("whitegrid")
    row_len = 3
    col_len = 3
    f, axes = prepare_subplots(row_len, col_len)
    extra_temperatures = temperatures
    extra_temperatures.extend([300,300])
    for i,T in enumerate(extra_temperatures):
        ax = axes[i/row_len, i%row_len]
        handles = []
        labels = []
        for ws in workspaces:
            sns.tsplot(data = collapsed_values[ws][T],  ci=[68, 95, 100], 
                   c= workspace_colors[ws], 
                   ax = ax)
            labels.append(ws)
            handles.append(mpatches.Patch(color=workspace_colors[ws], label=ws))
        if i%row_len == 0:
            ax.set_ylabel("Distance($\AA$)")
        if i/row_len == 2:
            ax.set_xlabel("Steps (accepted)")
        ax.set_title("T = %d"%T)
            
        ax.legend(handles, labels)
    plt.show()
            
    