'''
Created on 14/05/2015

@author: vgil
'''
import os

class MICAN(object):

    def __init__(self, exe_name = "mican_linux_64"):
        self.exe_name = exe_name
        
    def run(self, prot_reference, prot_2):
        cmd_str = "./%s %s %s -a alignment "%(self.exe_name, prot_reference, prot_2)
        os.system(cmd_str)
        select_lines = filter(lambda x: x[0] != "#" and len(x.split()) == 6,[line[0:-1] for line in open("alignment").readlines()])
        mapping = {}
        for select_line in select_lines:
            parts = select_line.split()
            res1 = parts[0]
            res2 = parts[3]
            
            if res1 == ".":
                print "No correspondence in reference for res %s of %s"%(res2, prot_2)
            
            elif res2 == ".":
                print "No correspondence in %s for res %s of reference"%(prot_2, res1)
            else:
                res1 = int(parts[0])
                res2 = int(parts[3])
                mapping[res1]=res2
        os.remove("alignment")
        return mapping

#MICAN().run("9WVG.initial.pdb", "1M2Z.pdb")

