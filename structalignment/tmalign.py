"""
Created on 14/05/2015

@author: vgil
"""
import os
import glob

class TMAlign:
    def __init__(self, exe_name = "TMalign"):
        self.exe_name = exe_name
    
    def run(self, prot_reference, prot_2):
        """
        Executes TM-align and returns the mapping of residue ids between the template and 
        the other protein.
        Creates the files "TM.sup*" and removes them.
        """
        cmd_str = "./%s %s %s -o TM.sup"%(self.exe_name, prot_reference, prot_2)
        os.system(cmd_str)
        select_lines = filter(lambda x: "select" in x and len(x.split())>2,[line[0:-1] for line in open("TM.sup").readlines()])
        mapping = {}
        for select_line in select_lines:
            parts = select_line.split()
            res1 = int(parts[1].split(":")[0])
            res2 = int(parts[2].split(":")[0])
            mapping[res1] = res2
            
        filelist = glob.glob("*.sup*")
        for f in filelist:
            os.remove(f)
        
        return mapping
        
#TMAlign().run("9WVG.initial.pdb", "3E7C.pdb")