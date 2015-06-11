"""
Created on 10/2/2015

@author: victor
"""
from anmichelpers.tools.tools import load_all_pdbs_ca, normalize_zero_one
from anmichelpers.tools.measure import rmsf
import matplotlib.pyplot as plt
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-s", dest="start_res", type = "int")
    parser.add_option("-p", dest="padding", action="store_true", default=False)
    
    parser.add_option("--normalize", dest="normalize", action = "store_true", default=False)
    
    (options, args) = parser.parse_args()
    
    for path in options.input.split(","):
        print "processing", path
        data, sources = load_all_pdbs_ca([{
                                       "source":path, 
                                       "base_selection":"name CA"   
                                       }])
        rmsf_array = rmsf(data.structure_ensemble)
        
        num_residues =  data.structure_ensemble.numResidues()
        if options.start_res is not None:
            res_range = range(options.start_res, options.start_res+num_residues)
            xlimits = [options.start_res, options.start_res+num_residues]
        else:
            res_range = range(num_residues)
            xlimits = [0, num_residues]
           
        if options.normalize:
#            plt.plot(normalize_zero_one(rmsf_array), )
            plt.plot(res_range, normalize_zero_one(rmsf_array), 'o-')
            ylimits = [0, 1]
        else:
            plt.plot(res_range, rmsf_array)
            
        if options.padding:
            if options.normalize:
                ylimits = [ylimits[0]-0.1, ylimits[1]+0.1]
            xlimits = [xlimits[0]-2, xlimits[1]+2]
            
        
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.xlabel('Residues')
    plt.ylabel('RMSF ($\AA$)')
    plt.show()
    