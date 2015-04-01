'''
Created on 26/03/2015

@author: user
'''

import prody
import sys
import matplotlib.pyplot as plt
from anmichelpers.tools.tools import norm

pdb = prody.parsePDB(sys.argv[1])

A = pdb.select("resid 25 and name CA").getCoordsets()
B =  pdb.select("resid 106 and name CA").getCoordsets()

dists = [norm(coords[0]) for coords in  B-A]

plt.plot(dists)
plt.ylabel('Distance')
plt.xlabel('Step (Accepted)')

plt.show()