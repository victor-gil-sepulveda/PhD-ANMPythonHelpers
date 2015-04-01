'''
Created on 19/03/2015

@author: user
'''
import numpy
from optparse import OptionParser
import matplotlib.pyplot as plt
import anmichelpers.tools.measure as measure
from xdg.Menu import tmp

def load_file(path):
    v = numpy.loadtxt(path)
    if len(v.shape) == 1:
        return numpy.array([v])
    else:
        return v

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--i1", dest="input1")
    parser.add_option("--i2", dest="input2")
    parser.add_option("-f","--from", type= "int", dest="_from")
    parser.add_option("-t","--to", type= "int", dest="to")
    parser.add_option("-p","--plot_type", default="normal", dest="plot_type")
    parser.add_option("-s","--skip_step", action= "store_true", dest="skip_step")
    
    (options, args) = parser.parse_args()
    
    plot_types = ["ccdist", "absdiff", "diff", "normal"]
    assert options.to >= options._from,  "[ERROR] 'from' value is bigger than 'to'. "
    assert options.plot_type in plot_types, "[ERROR] plot type ('-p','--plot_type') must be one of %s"%str(plot_types)

    result = None
    if options.plot_type == "diff" or options.plot_type == "absdiff":
        v1 = load_file(options.input1)
        v2 = load_file(options.input2)
        assert len(v1[0]) == len(v2[0]),"[ERROR] arrays must have the same length (%s vs %s)."
        #can have different number of rows
        if len(v1) != len(v2):
            min_c = min(len(v1),len(v2))
            v1 = v1[:min_c,:]
            v2 = v2[:min_c,:]
        result = v2-v1
        if options.plot_type == "absdiff":
            result = abs(result)
        
    if options.plot_type == "ccdist":
        v1 = load_file(options.input1)
        v2 = load_file(options.input2)
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
        
    elif options.plot_type == "normal":
        result = load_file(options.input1)

    if options.skip_step:
        # skip first number
        result = result[:,1:]
    
    number_of_plots = options.to-options._from+1 
    subplot_shape = (number_of_plots,1)
    
    plt.title('----')
    for i in range(number_of_plots):
        ax = plt.subplot2grid(subplot_shape, (i,0))
        plt.plot(range(len(result[i])), result[i])
    plt.ylabel("value")
    plt.xlabel("id")
    plt.show()
    
