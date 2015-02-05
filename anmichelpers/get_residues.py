import sys
from prody import *

pdb = sys.argv[1]
structure = parsePDB(pdb)

calphas = parsePDB(pdb, subset='ca')
writePDB(pdb + '_ca.pdb', calphas)

hv = structure.getHierView()
a_chain = hv['A']

atoms = hv.getAtoms()

n = 0
indexes = []
for i, a in enumerate(atoms.iterAtoms()):
    if a.getName() == 'CA':
        # print str(n) + ' - ' + str(i) + ' - ' + a.getName()
        n += 1
        indexes.append(i)

print indexes
# a_residues = list(a_chain)

# for i, r in enumerate(a_residues):
#     print str(r.getResnum()) + ' - ' + r.getResname()
#     print r.getNames()
