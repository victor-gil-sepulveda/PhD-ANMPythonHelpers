import sys
import prody 

pdb = sys.argv[1]
backbone = prody.parsePDB(pdb, subset='bb')

prody.writePDB(pdb + '_backbone.pdb', backbone)
