"""
Plots the indicated precalculated rmsfs all together. From "rmsf_files" plots the 4 that are closer to
a reference rmsf. It will also plot the rmsfs in "control_rmsf_files" with a thicker line. 
"""

import numpy
import matplotlib.pyplot as plt
import math 

common = "/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_CC/Benchmark/"
rmsf_files = {
common+"results/300_0_25.pdb.rmsf" :{"T":300, "Ca_dist":0.25},
common+"results/300_0_66.pdb.rmsf" :{"T":300, "Ca_dist":0.66},
common+"results/300_1_08.pdb.rmsf" :{"T":300, "Ca_dist":1.08},
common+"results/300_1_5.pdb.rmsf"  :{"T":300, "Ca_dist":1.5},
common+"results/866_0_25.pdb.rmsf" :{"T":866, "Ca_dist":0.25},
common+"results/866_0_66.pdb.rmsf" :{"T":866, "Ca_dist":0.66},
common+"results/866_1_08.pdb.rmsf" :{"T":866, "Ca_dist":1.08},
common+"results/866_1_5.pdb.rmsf"  :{"T":866, "Ca_dist":1.5},
common+"results/1432_0_25.pdb.rmsf":{"T":1432, "Ca_dist":0.25},
common+"results/1432_0_66.pdb.rmsf":{"T":1432, "Ca_dist":0.66},
common+"results/1432_1_08.pdb.rmsf":{"T":1432, "Ca_dist":1.08},
common+"results/1432_1_5.pdb.rmsf" :{"T":1432, "Ca_dist":1.5},
common+"results/2000_0_25.pdb.rmsf":{"T":2000, "Ca_dist":0.25},
common+"results/2000_0_66.pdb.rmsf":{"T":2000, "Ca_dist":0.66},
common+"results/2000_1_08.pdb.rmsf":{"T":2000, "Ca_dist":1.08},
common+"results/2000_1_5.pdb.rmsf" :{"T":2000, "Ca_dist":1.5}
}

control_rmsf_files = {
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Suwipa/mdgb1.pdb.rmsf":{"T":"MD", "Ca_dist":"Suw1"},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Suwipa/mdgb2.pdb.rmsf":{"T":"MD", "Ca_dist":"Suw2"},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Suwipa/mdgb5.pdb.rmsf":{"T":"MD", "Ca_dist":"Suw5"},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_300_nosdmin/all.pdb.rmsf":{"T":"IC300", "Ca_dist":"no_sdmin"},
"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/ANM_IC/T_500_nosdmin/all.pdb.rmsf":{"T":"IC500", "Ca_dist":"no_sdmin"},

"/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Shaw/pro_noh_md.pdb.rmsf":{"T":"MD", "Ca_dist":"Shaw"}
}

REFERENCE_RMSF = numpy.loadtxt("/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/ANMIC/MD/Shaw/pro_noh_md.pdb.rmsf")[:-1]

def rms(one, two):
    return math.sqrt( numpy.dot(one-two, one-two)/len(one))

def get_best(number_of_best, reference_rmsf, tested_rmsfs_dic):
    rmsds = []
    for filename in tested_rmsfs_dic:
        rmsds.append((rms(reference_rmsf, numpy.loadtxt(filename)), filename))
    rmsds.sort()
    
    for rmsf_score in rmsds:
        print "%.2f"%rmsf_score[0], rmsf_score[1]
        
    return numpy.array([rmsd for _, rmsd in rmsds])
    

best = get_best(2,REFERENCE_RMSF, rmsf_files)[0:4]

labels = []
# for filename in best:
#     label ="%d_%.2f"%(rmsf_files[filename]["T"],rmsf_files[filename]["Ca_dist"])
#     rmsf = numpy.loadtxt(filename)
#     plt.plot(rmsf, label = label)

for filename in control_rmsf_files:
    label = "%s_%s"%(control_rmsf_files[filename]["T"],control_rmsf_files[filename]["Ca_dist"])
    rmsf = numpy.loadtxt(filename)
    lines = plt.plot(rmsf, label = label)
    plt.setp(lines, linewidth=2)

plt.legend(loc=9, ncol = 4)

plt.show()
