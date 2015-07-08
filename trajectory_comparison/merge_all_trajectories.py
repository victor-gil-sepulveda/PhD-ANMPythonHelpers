"""
Uses the trajectory merger to merge all small trajectories inside a folder. Does it for all 
the available folders.
"""

import os
rootdir = '.'

for subdir, dirs, files in os.walk(rootdir):
    if not "results" in subdir:
        print subdir
        command = "python /home/victor/git/PhD-PDBMerger/PDBTrajectoryMerger.py -p run_trajectory_ -d %s -o results/%s.pdb"%(subdir,subdir)
        os.system(command)
        
