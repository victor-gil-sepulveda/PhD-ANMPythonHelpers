'''
Created on 19/03/2015

@author: user
'''
import numpy
from optparse import OptionParser
import matplotlib.pyplot as plt
import anmichelpers.tools.measure as measure
from pyRMSD.RMSDCalculator import RMSDCalculator

def load_file(path):
    try:
        v = numpy.loadtxt(path)
        if len(v.shape) == 1:
            return numpy.array([v])
        else:
            return v
    except:
        print "[ERROR] Impossible to load %s"%path
        exit()

def needs_one_file(i1):
    assert i1 is not None, "An input file is needed"
    return load_file(i1)

def needs_two_files(i1, i2):
    needs_one_file(i1)
    needs_one_file(i2)
    return load_file(i1), load_file(i2)

def can_have_different_numberof_rows(v1, v2):
    #can have different number of rows
    if len(v1) != len(v2):
        min_c = min(len(v1),len(v2))
        v1 = v1[:min_c,:]
        v2 = v2[:min_c,:]
    return v1, v2

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--i1", dest="input1")
    parser.add_option("--i2", dest="input2")
    parser.add_option("-f","--from", type= "int", dest="_from")
    parser.add_option("-t","--to", type= "int", dest="to")
    parser.add_option("-p","--plot_type", default="normal", dest="plot_type")
    parser.add_option("-s","--skip_step", action= "store_true", dest="skip_step")
    
    (options, args) = parser.parse_args()
    
    plot_types = ["ccdist", "absdiff", "diff", "rmsd", "normal"]
    assert options.to >= options._from,  "[ERROR] 'from' value is bigger than 'to'. "
    assert options.plot_type in plot_types, "[ERROR] plot type ('-p','--plot_type') must be one of %s"%str(plot_types)

    result = None
    if options.plot_type == "diff" or options.plot_type == "absdiff":
        v1,v2 = needs_two_files(options.input1, options.input2)
        
        assert len(v1[0]) == len(v2[0]),"[ERROR] arrays must have the same length (%s vs %s)."
        
        v1,v2 = can_have_different_numberof_rows(v1, v2)
        
        result = v2-v1
        
        if options.plot_type == "absdiff":
            result = abs(result)
        
        if options.skip_step:
            # skip first number
            result = result[:,1:]
        
    if options.plot_type == "ccdist":
        v1,v2 = needs_two_files(options.input1, options.input2)
        
        assert len(v1[0]) == len(v2[0]),"[ERROR] arrays must have the same length (%s vs %s)."
        
        if len(v1) != len(v2):
            min_c = min(len(v1),len(v2))
            v1 = v1[:min_c,:]
            v2 = v2[:min_c,:]
            
        result = v2-v1
        tmp_result = []
        for r in result:
            tmp_result.append(measure.calculate_mode_magnitudes(r))
        result = numpy.array(tmp_result)
        if options.skip_step:
            # skip first number
            result = result[:,1:]
            
    if options.plot_type == "rmsd":
        v1,v2 = needs_two_files(options.input1, options.input2)
        v1,v2 = v1[:,1:], v2[:,1:]
        v1,v2 = can_have_different_numberof_rows(v1, v2)
        
        result = []
        for i in range(len(v1)):
            coordset1 = numpy.resize(v1[i], (len(v1[i])/3,3))
            coordset2 = numpy.resize(v2[i], (len(v2[i])/3,3))
            coordsets = numpy.array([coordset1, coordset2])
            calculator = RMSDCalculator("QCP_SERIAL_CALCULATOR", coordsets)
            result.append(calculator.pairwise(0,1))
    
    elif options.plot_type == "normal":
        result = needs_one_file(options.input1)

        if options.skip_step:
            # skip first number
            result = result[:,1:]
    
    if options.plot_type in ["ccdist", "absdiff", "diff", "normal"]:
        number_of_plots = options.to-options._from+1 
        subplot_shape = (number_of_plots,1)
        plt.title('----')
        for i in range(number_of_plots):
            ax = plt.subplot2grid(subplot_shape, (i,0))
            plt.plot(range(len(result[i])), result[i])
        plt.ylabel("value")
        plt.xlabel("id")
        plt.show()
    elif options.plot_type in ["rmsd"]:
        plt.plot(result)
        plt.ylabel("value")
        plt.xlabel("id")
        plt.show()
