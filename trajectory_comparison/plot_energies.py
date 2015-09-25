'''
Created on 21/9/2015

@author: victor
'''
import sys
import numpy
import matplotlib.pyplot as plt
from trajectory_comparison import tools

if __name__ == '__main__':
    energies_file = sys.argv[1]
    energy_handler = open(energies_file, "r")
    
    if len(sys.argv[1]) > 2:
        skip = int(sys.argv[2])
    else:
        skip = 0
    
    energies = []
    traj_energy = []
    for line in energy_handler:
        ml = line.strip('\n')
        if ml != "###":
            traj_energy.append(float(ml))
        else:
            energies.append(traj_energy[skip:-1])
            traj_energy = []
    energies.append(traj_energy)
    energy_handler.close()
    
    colors = tools.define_tableau20_cm()
    
    for i in range(len(energies)):
        ax = plt.plot(energies[i], color=colors[i%len(colors)])
        plt.ticklabel_format(useOffset=False)
    plt.savefig("%s.svg"%sys.argv[1])