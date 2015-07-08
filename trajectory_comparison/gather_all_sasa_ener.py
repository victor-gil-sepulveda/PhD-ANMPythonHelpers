"""
Simply walks the directories to gather all SASA files and merge them in one file that will be
stored somewhere else.
"""

import os
rootdir = '.'

for subdir, dirs, files in os.walk(rootdir):
    if not "results" in subdir:
        for filename in files:
            if ".ener" in filename:
                os.system("cat %s >> results/%s.ener"%(os.path.join(subdir,filename), subdir))
            if ".sasa" in filename:
                os.system("cat %s >> results/%s.sasa"%(os.path.join(subdir,filename), subdir))
        
