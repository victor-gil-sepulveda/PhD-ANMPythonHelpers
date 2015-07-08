"""
Calculates the acceptances for all the pdbs in a folder and leaves the result in the 
same folder the original file was.  
"""

import os
rootdir = '.'

for subdir, dirs, files in os.walk(rootdir):
    if not "results" in subdir:
        command = "python /home/victor/git/PhD-ANMPythonHelpers/trajectory_full_analysis.py --report %s/run_report_* -i %s"%(subdir,subdir)
        os.system(command)
        
