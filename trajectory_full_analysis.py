"""
Created on 29/06/2015

@author: user
@brief: Performs some automatic analysis over the trajectories specified in a 
file.
"""

import os
import glob
import numpy
from optparse import OptionParser

mdtraj_accessible = True

try:
    import mdtraj
except:
    mdtraj_accessible = False

anmhelpers_accessible = True
try:    
    from anmichelpers.tools.tools import load_all_pdbs_ca
    from anmichelpers.tools.measure import rmsf
except:
    anmhelpers_accessible = False

def analyze_trajectory(traj_path, do_sasa, do_rgyr, do_rmsf, report_pattern):
    #Shrake, A; Rupley, JA. (1973) J Mol Biol 79 (2): 351--71.
    
    trajectory = None
    data = None
    
    if do_sasa:
        print "Calculating SASA ..."
        if trajectory is None:
            trajectory = mdtraj.load(traj_path)
        sasa = mdtraj.shrake_rupley(trajectory, mode = 'residue').sum(axis=1)
        numpy.savetxt(traj_path+'.sasa', sasa)
    
    if do_rgyr:
        print "Calculating Radius of Gyration ..."
        if trajectory is None:
            trajectory = mdtraj.load(traj_path)
        rgyr = mdtraj.compute_rg(trajectory)
        numpy.savetxt(traj_path+'.rgyr', rgyr)
    
    if trajectory is not None:
        del trajectory
    
    if do_rmsf:
        print "Calculating RMSF ..." 
        data, _ = load_all_pdbs_ca([{
                                       "source": traj_path, 
                                       "base_selection":"name CA"   
                                    }])
        rmsf_array = rmsf(data.structure_ensemble)
        numpy.savetxt(traj_path+'.rmsf', rmsf_array)

    if data is not None:
        del data
    
    if report_pattern != "":
        print "Extracting acceptance and energies from report files with pattern %s ..."%report_pattern
        directory, _ = os.path.split(traj_path)
        
        files = glob.glob(os.path.join(directory, report_pattern))
        
        assert len(files)!=0, "No report file with pattern %s found inside %s"%(report_pattern,
                                                                                 directory)
        all_accepted = []
        all_total = []
        all_energies = []
        for report_file in files:
            total, accepted, energies = process_report_file(report_file)
            all_total.append(total)
            all_accepted.append(accepted)
            all_energies.extend(energies)
        total = numpy.sum(all_total)
        accepted = numpy.sum(all_accepted)
        acceptance = accepted / total
        numpy.savetxt(traj_path+'.acc', [acceptance], fmt = "%.4f ")
        numpy.savetxt(traj_path+'.ener', all_energies)

def process_report_file(report_file):
    handler = open(report_file)
    lines = handler.readlines()
    header = lines[0].split()
    total_steps_index = header.index("Step")
    accepted_steps_index = header.index("AcceptedSteps")
    energy_index = header.index("Energy")
    handler.close()
    
    values = numpy.loadtxt(report_file,comments = "#").T
    
    return values[total_steps_index][-1], values[accepted_steps_index][-1], values[energy_index]
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="traj_path")
    parser.add_option("--sasa", dest="do_sasa", action="store_true", default=False)
    parser.add_option("--rgyr", dest="do_rgyr", action="store_true", default=False)
    parser.add_option("--rmsf", dest="do_rmsf", action="store_true", default=False)
    parser.add_option("--report", dest="report", default="")
    
    (options, args) = parser.parse_args()

    assert  options.traj_path is not None, "It is mandatory to choose a valid trajectory file (-i option)"
    
    if  (options.do_sasa or options.do_rgyr) and not mdtraj_accessible: 
        print "It was not possible to load the 'mdtraj' module. Using the '--sasa' or '--rgyr' options is not possible."
        exit()
        
    if options.do_rmsf and not anmhelpers_accessible:
        print "It was not possible to load the 'anmichelpers' module. Using the '--rmsf' option is not possible."
        exit()
    
    analyze_trajectory(options.traj_path, 
                       options.do_sasa,
                       options.do_rgyr,
                       options.do_rmsf,
                       options.report)

#import matplotlib.pyplot as plt
#    plt.hist(sasa)
#    plt.show()