'''
Created on 20/1/2016

@author: victor
'''

import seaborn as sns
import sys
import numpy
import matplotlib.pyplot as plt

if __name__ == '__main__':
    times_1 = sys.argv[1]
    times_1_name = sys.argv[2]
    times_2 = sys.argv[3]
    times_2_name = sys.argv[4]
    
    sns.set_style("whitegrid")
    
    for label, data in [(times_1_name, times_1),(times_2_name, times_2)]:
        sns.distplot(numpy.loadtxt(data).T[1], label = label)
    plt.legend()
    plt.show()


    
    