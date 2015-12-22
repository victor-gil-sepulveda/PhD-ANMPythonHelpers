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


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-d", dest="data")
    parser.add_option("--results", dest="results")
    
    (options, args) = parser.parse_args()
    if not options.results:  
        parser.error('Results folder not given')
        
    folders = glob.glob("CC_*")
    workspaces = ["cc_closed_rmsg_dispf","cc_open_rmsg_dispf"]
    workspace_titles = {"cc_closed_rmsg_dispf":"Open",
                        "cc_open_rmsg_dispf":"Closed"}
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
        distances = pickle.load(options.data)
    
    for wp in distances:
        for T in sorted(distances[wp].keys()):
            for v1,v2 in distances[wp][T].keys():
                plt.plot(distances[wp][T][v1,v2], label = "%d"%(T))    
        plt.title(workspace_titles[wp])
        plt.show()
            
    