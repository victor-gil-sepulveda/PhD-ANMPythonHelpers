'''
Created on 3/3/2016

@author: victor
'''
import os.path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy
import sys

sns.set_style("whitegrid")

if __name__ == '__main__':
    labels = {
              "cc":"CC",
              "ic":"IC",
              "md":"MD"
              }
    
    files = {
             "src": {
                     "cc":"merged.pdb.residue_8_to_266",
                     "ic":"merged.pdb.residue_8_to_266",
                     "md":"md_no_het.pdb.residue_8_to_266"
                     },
             "ubi" : {"cc":"merged.pdb.residue_1_to_73",
                      "ic":"merged.pdb.residue_1_to_73",
                      "md":"opls_skip10.pdb.residue_1_to_73"
                      }
             }
    
    extension = sys.argv[1]
    for prot_folder in files:
        for sim_type in files[prot_folder]:
            filename = os.path.join(prot_folder, sim_type, files[prot_folder][sim_type])
            sasa = numpy.loadtxt("%s.%s"%(filename,extension))
            plt.hist(sasa, 100, label = labels[sim_type], normed = True, alpha = 0.75)
        plt.title(prot_folder)
        plt.legend()
#         plt.show()
        plt.savefig("%s_%s.svg"%(prot_folder,extension))