import sys
from prody import *

pdb = sys.argv[1]
backbone = parsePDB(pdb, subset='bb')

writePDB(pdb + '_backbone.pdb', backbone)
