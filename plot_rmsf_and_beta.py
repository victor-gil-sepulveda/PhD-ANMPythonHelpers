"""
Created on 10/2/2015

@author: victor
"""
from anmichelpers.tools.tools import load_all_pdbs_ca, get_all_betas, normalize,\
    normalize_zero_one
from anmichelpers.tools.measure import rmsf
import matplotlib.pyplot as plt

if __name__ == '__main__':
    data, sources = load_all_pdbs_ca([{
                                   "source":"/home/victor/Escritorio/Christoph/structures/*.pdb", 
                                   "base_selection":"name CA"   
                                   }])
    betas = get_all_betas(sources)
    
    rmsf_array = rmsf(data.structure_ensemble)
    
    plt.plot(normalize_zero_one(betas.mean(axis = 0)))
    plt.plot(normalize_zero_one(rmsf_array))
    plt.show()
    