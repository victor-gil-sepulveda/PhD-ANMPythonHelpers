from anmichelpers.parsers.pronmd import ProdyNMDParser
from prody.dynamics.compare import calcCumulOverlap
from prody.dynamics.nma import NMA
from anmichelpers.comparison.comparison import cumulative_overlap, overlap
from anmichelpers.tools.tools import norm


eval1, evec1, _ = ProdyNMDParser.read("/home/victor/Escritorio/NMA_modes_charachterization/cc/CC_prot_ubi_start.pdb_cutt_9/info/normalized_modes.1.nmd")
eval2, evec2, _ = ProdyNMDParser.read("/home/victor/Escritorio/NMA_modes_charachterization/ic/IC_cutt_9_prot_ubi_start.pdb/info/normalized_modes_cc.1.nmd")


for i in range(len(evec2)-1):
    for j in range(i, len(evec2)):
        print overlap(evec2[i],evec2[j]/norm(evec2[j])),
    print

for i in range(len(evec1)-1):
    for j in range(i, len(evec1)):
        print overlap(evec1[i],evec1[j]),
    print
    

nma1 = NMA("cc")
nma1.setEigens(evec1.T, eval1)
nma2 = NMA("ic")
nma2.setEigens(evec2.T, eval2)
print calcCumulOverlap(nma1[1], nma2[:10])

print cumulative_overlap(evec1[1], evec2[0:10])
