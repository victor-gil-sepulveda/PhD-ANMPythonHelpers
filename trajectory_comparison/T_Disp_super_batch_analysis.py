'''
Created on 14/9/2015

@author: victor
'''
import sys
import tools
from analysis_pipeline import execute_pipeline
import os

def get_folder_name(t,disp, it = ""):
    if it == "":
        return "%s_%s"%(t,disp)
    else:
        return "%s_%s"%(get_folder_name(t, disp),it)


def get_folders_for_analysis(temp_and_disp_file):
    ts_and_disp_file = temp_and_disp_file
    lines = open(ts_and_disp_file,"r").readlines()
    temperatures = lines[0].split()
    displacements = lines[1].split()
    use_iterations = True
    try:
        iterations = lines[2].split()
    except:
        use_iterations = False
    folders = []
    for T in temperatures:
        for disp in displacements:
            if use_iterations:
                for it in iterations:
                    folders.append( (get_folder_name(T, disp, it),(T,disp,it)))
            else:
                folders.append((get_folder_name( T, disp),(T,disp,-1)))
    return folders

if __name__ == '__main__':
    folders = get_folders_for_analysis(sys.argv[1])
    
    template_parameters_file = sys.argv[2]
    parameters =  tools.load_control_json(template_parameters_file)
    
    BASE_PATH = "Results"
    for folder, _ in folders:
        print "Executing pipeline for folder: ", folder
        
        # Change script
        parameters["workspace"] = os.path.join(BASE_PATH, folder)
        parameters["merge"]["working_directory"] = os.path.join(os.getcwd(),folder)
        parameters["merge"]["merged_file"] = "%s.merged.pdb"%folder
        parameters["analyze"]["acceptance"]["working_directory"] = folder
        
        # Execute
        execute_pipeline(parameters)
    
    
    
