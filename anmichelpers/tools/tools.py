"""
Created on 9/2/2015

@author: victor
"""

import  numpy
import math
import prody
try:
    from pyproct.data.handler.sourceGenerator import SourceGenerator
    from pyproct.data.handler.protein.proteinEnsembleDataLoader import ProteinEnsembleDataLoader
except:
    print "[WARNING] pyProCT was not found. Some functions cannot be used"
    
def norm(v):
    """
    Numpy compliant norm implementation.
    
    @param v: The vector used to calculate the norm.
    
    @return: A norm or an array of norms.
    """
    if len(v.shape) == 1:
        return numpy.sqrt(numpy.dot(v,v))
    elif len(v.shape) == 2:
        norms = []
        for i in range(len(v)):
            norms.append(norm(v[i]))
        return norms
    else:
        return None

def frec_from_eigenvalue(e_val):
    """
    Calculates the proportional frequency of a given eigenvalue (if it comes from a 
    vibrational study).
    
    @param e_val: The eigenvalue.
    
    @return: The computed frequency (no units).
    """
    return e_val / (2*math.pi)

def ensure_modes_layout(modes):
    """
    If the layout of the modes is flat, it converts it to a (M,N,3) layout.
    
    @param modes: [In/Out] A numpy array containing all the modes. 
    
    @return: The same numpy array with a (M,N,3) layout or (N,3) 
    """
    if len(modes.shape) == 3:
        return modes
    elif len(modes.shape) == 2:
        number_of_modes = len(modes)
        number_of_nodes = modes.shape[1] / 3
        return numpy.reshape(modes, (number_of_modes, number_of_nodes, 3))
    else:
        raise ValueError("The array has an unexpected size")
    
def load_all_pdbs_ca(pdb_list):
    """
    Loads a list of pdbs in pyproct format (this includes the use of globs and 'base_selection'.
    
    @param pdb_list: A list of pdbs in pyproct format.
    
    @return: The pyproct data object and the list of sources (prody pdb structure -> data.structure_ensemble
    source from pyproct source -> s.source["source"] )
    """
    class MockParams:
        def __init__(self):
            pass
        def get_value(self,a,b):
            return ""
    sources = SourceGenerator(pdb_list).source_list
    loader = ProteinEnsembleDataLoader(MockParams())
    for source in sources:
        loader.load(source)
    # Retrieve the data object
    data = loader.close()
    return data, sources

def get_all_betas(sources):
    """
    Loads CA temperature factors from a list of pyproct sources.
    
    @return: A matrix with all the beta factors.
    """
    betas = []
    for s in sources:
        pdb = prody.parsePDB(s.source["source"]).select("name CA")
        betas.append(pdb.getBetas())
    
    betas = numpy.array(betas)
    
    mean_betas = betas.mean(axis = 0)
    
    for beta_array in betas:
        for i in range(len(beta_array)):
            if  beta_array[i] == 0.0:
                print "Zero beta value @ %d; exchanging with mean."%i
                beta_array[i] = mean_betas[i]
        
    return betas

def normalize(v):
    max_val = max(v)
    return v / abs(max_val)

def normalize_zero_one(v):
    max_val = max(v)
    min_val = min(v)
    val_range = max_val - min_val
    return (v - min_val) / val_range

def is_int(this_str):
    try:
        int(this_str)
        return True
    except ValueError:
        return False

