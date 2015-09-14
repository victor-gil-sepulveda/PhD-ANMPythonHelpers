"""
Created on Sep 11, 2015

@author: victor
"""

import sys
import tools
import os
        
MERGING_COMMAND = "python %(script_location)s -d %(working_directory)s -p %(trajectory_prefix)s -o %(merged_file)s"

def execute_pipeline(parameters):
    # Retrieve workspace
    tools.create_directory(parameters["workspace"])
    
    # Merge files
    parameters["merge"]["merged_file"] = os.path.join(parameters["workspace"], parameters["merge"]["merged_file"])
    if parameters["merge"]["do"]:
        command = MERGING_COMMAND%parameters["merge"]
        print "Issuing >> ", command
        os.system(command)
        if parameters["merge"]["delete_parts"]:
            os.system("rm %s*"%(os.path.join(parameters["merge"]["working_directory"], parameters["merge"]["trajectory_prefix"])))
        if parameters["merge"]["compress_logs"]:
            os.system("tar -zcvf %(folder)s/logs.tar.gz %(folder)s/log_* --remove-files"%({"folder":os.path.join(parameters["merge"]["working_directory"])}))
    
    # Perform all analyses
    if parameters["analyze"]["do"]:
        args = ["-i", parameters["merge"]["merged_file"]]
        
        if "sasa_rgyr" in parameters["analyze"]:
            if parameters["analyze"]["sasa_rgyr"]["sasa"]:
                args.append("--sasa-vmd")
            if parameters["analyze"]["sasa_rgyr"]["rgyr"]:
                args.append("--rgyr-vmd")
            if parameters["analyze"]["sasa_rgyr"]["selection"] != "":
                args.append("--vmd-sel")
                args.append('"%s"'%parameters["analyze"]["sasa_rgyr"]["selection"])
                
        if parameters["analyze"]["rmsf"]:
            args.append("--rmsf")
            
        if "acceptance" in parameters["analyze"]:
            args.append("--report")
            args.append(parameters["analyze"]["acceptance"]["report_glob"])
            args.append("--report-dir")
            args.append(parameters["analyze"]["acceptance"]["working_directory"])
        
        command = "python %s %s"%(parameters["analyze"]["script_location"], 
                           " ".join(args))
        print "Issuing >> ", command
        os.system(command)
    
    if parameters["compare"]["do"]:
        try:
            selection = parameters["analyze"]["sasa_rgyr"]["selection"]
        except:
            selection = "all"
            
        base_path_for_sasa_rgyr = "%s.%s"%(parameters["merge"]["merged_file"],
                              selection.replace(" ","_"))
        if "sasa" in parameters["compare"]:
            parameters["compare"]["sasa"]["input"] = "%s.sasa"%(base_path_for_sasa_rgyr)
            command = "python %(script_location)s %(reference)s %(input)s save"%(parameters["compare"]["sasa"])
            print "Issuing >> ", command
            os.system(command)
            
        if "rgyr" in parameters["compare"]:
            parameters["compare"]["rgyr"]["input"] = "%s.rgyr"%(base_path_for_sasa_rgyr)
            command = "python %(script_location)s %(reference)s %(input)s save"%(parameters["compare"]["rgyr"])
            print "Issuing >> ", command
            os.system(command)
        
        if "rmsf" in parameters["compare"]:
            parameters["compare"]["rmsf"]["input"] = "%s.rmsf"%(parameters["merge"]["merged_file"])
            command = "python %(script_location)s %(reference)s %(input)s save"%(parameters["compare"]["rmsf"])
            print "Issuing >> ", command
            os.system(command)

if __name__ == '__main__':
    
    # Get parameters
    parameters = tools.load_control_json(sys.argv[1])
    
    # Do all the stuff
    execute_pipeline(parameters)
    
    
            
            