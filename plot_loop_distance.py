'''
Created on 26/03/2015

@author: user
'''

import prody
import sys
import matplotlib.pyplot as plt
from anmichelpers.tools.tools import norm

pdb = prody.parsePDB(sys.argv[1])

if sys.argv[2] == "closed":
    A = pdb.select("resid 25 and name CA").getCoordsets()
    B =  pdb.select("resid 106 and name CA").getCoordsets()

elif sys.argv[2] == "open":
    A = pdb.select("resid 31 and name CA").getCoordsets()
    B =  pdb.select("resid 112 and name CA").getCoordsets()
else:
    print "You need to state if you'r plotting the 'open' or 'closed' conformations"
    exit(-1)

dists = [norm(coords[0]) for coords in  B-A]

plt.plot(dists)
plt.ylabel('Distance')
plt.xlabel('Step (Accepted)')

plt.show()