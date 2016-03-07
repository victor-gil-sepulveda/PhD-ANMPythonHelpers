'''
Created on 20/1/2016

@author: victor
'''

import seaborn as sns
import sys
import numpy
import matplotlib.pyplot as plt

if __name__ == '__main__':
    plots = []
    for i in range(1,(len(sys.argv)/2)+1):
        plots.append((sys.argv[2*i],sys.argv[2*i-1]))
    
    sns.set_style("whitegrid")
    
    for label, data_path in plots:
        
        data = numpy.loadtxt(data_path)
        
        if len(data.shape) > 1:
            plt.hist(data.T[1], 50, label = label, normed=1, range = (0,150.), alpha = 0.7)
        else:
            # remove anything beyond 150 to improve visualization
            filtered_data = data[data<=150]
            print "FILTERED: %.3f %%"%(100. -(len(filtered_data)*100./ len(data)))
            plt.hist(filtered_data, 50, label = label, normed=1, range = (0,150.), alpha = 0.7)
#         sns.distplot(numpy.loadtxt(data).T[1], label = label)
    plt.xlabel("Time (s)")
    plt.title("Performance of IC vs CC methods")
    plt.legend()
    plt.show()


    
    