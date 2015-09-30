"""
Created on 30/9/2015

@author: victor
"""
import sys
from trajectory_comparison.T_Disp_super_batch_analysis import get_folders_for_analysis
import os
import glob
import numpy

def get_num_models(merged_pdb):
    models = 0
    handler = open(merged_pdb,"r")
    for line in handler:
        if "MODEL" == line[0:5]:
            models += 1
    handler.close()
    return models

if __name__ == '__main__':
    folders = get_folders_for_analysis(sys.argv[1])
    base_path = sys.argv[2]

    results = {}
    expected_data = ["rgyr.jsd",
                     "sasa.jsd",
                     "rms_rmsfs",
                     "acc",
                     "models_per_h_node"]
    
    ordered_data = ["T","disp","it"]
    ordered_data.extend(expected_data)
    
    num_processors = int(sys.argv[3])
    num_hours = int(sys.argv[4])
    for folder, data in folders:
        path = os.path.join(sys.argv[2], folder)
        print "Summarizing folder: ", path
        
        key = (int(data[0]), data[1], data[2])
        results[key] = {"T":data[0],"disp":data[1],"it":data[2]}
        
        for ext in expected_data:
            files = glob.glob(os.path.join(path, "*.%s"%ext))
            if len(files) != 1:
                print "PROBLEM in %s finding files with extension %s. Num files: %d"%(path, ext, len(files))
            else:
                results[key][ext] = "%.3f"%numpy.loadtxt(files[0])
        try:
            merged_pdb = glob.glob(os.path.join(path, "*.pdb"))[0]
            acc_steps = get_num_models(merged_pdb)
            total_steps = acc_steps / float(results[key]["acc"])
            results[key]["models_per_h_node"] = "%.3f"%(total_steps / (num_processors*num_hours))
        except:
            pass
        
    all_ordered_keys = sorted(results.keys())
    for key in all_ordered_keys:
        for data_type in ordered_data:
            try:
                print "%6s  "%results[key][data_type],
            except KeyError:
                print "%6s  "%"---",
        print
        
        