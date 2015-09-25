"""
Created on 29/06/2015

@author: vgil
@brief: Performs some automatic analysis over the trajectories. The kind of analysis to be done are 
specificied using the command line (rmsf, sasa, radius of gyration or a report analysis if it is a 
PELE result)
"""

import os
import glob
import numpy
from optparse import OptionParser
from anmichelpers.tools.measure import ca_rmsf

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


VMD_MEASURE_SCRIPT = """# from http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/7502.html
mol load pdb  %s
set outfile [open %s w]
set nf [molinfo top get numframes]
set all [atomselect top "%s"]
for {set i 0} {$i<$nf} {incr i} {
    $all frame $i
    $all update
    set sasa [measure %s $all]
    puts $outfile "$sasa" 
}
close $outfile
quit
"""

VMD_MEASURE_SASA_AND_RGYR_SCRIPT = """mol load pdb  %s
set outfile_sasa [open %s w]
set outfile_rgyr [open %s w]
set nf [molinfo top get numframes]
set all [atomselect top "%s"]
for {set i 0} {$i<$nf} {incr i} {
    $all frame $i
    $all update
    set sasa [measure sasa 1.4 $all]
    set rgyr [measure rgyr $all]
    puts $outfile_sasa "$sasa" 
    puts $outfile_rgyr "$rgyr" 
    
}
close $outfile_sasa
close $outfile_rgyr
quit
"""

def execute_vmd_script(vmd_code):
    open("tmp_vmd_script","w").write(vmd_code)
    os.system("vmd -dispdev none -e tmp_vmd_script > vmd_out")

def get_vmd_array():
    a = numpy.loadtxt("tmp_vmd_out")
    os.system("rm tmp_vmd_script tmp_vmd_out vmd_out")
    return a

def calculate_sasa_with_vmd(pdb_file, outfile, vmd_selection):
    """
    It executes a vmd script in order to calculate trajectory' SASA (Angstrom^2). 
    'mdtraj' calculates it in nm^2.
    """
    vmd_code = VMD_MEASURE_SCRIPT%(pdb_file, "tmp_vmd_out", vmd_selection, "sasa 1.4")
    
    execute_vmd_script(vmd_code)
    sasa = get_vmd_array()
    
    sasa = sasa /100.
    numpy.savetxt(outfile, sasa, "%.4f")

def calculate_rgyr_with_vmd(pdb_file, outfile, vmd_selection):
    """
    http://www.ks.uiuc.edu/Research/vmd/vmd-1.9/ug/node135.html
    """
    vmd_code = VMD_MEASURE_SCRIPT%(pdb_file, "tmp_vmd_out", vmd_selection, "rgyr")
    
    execute_vmd_script(vmd_code)
    rgyr = get_vmd_array()
    
    numpy.savetxt(outfile, rgyr, "%.4f")
    
def calculate_sasa_and_rgyr_with_vmd(pdb_file, sasa_outfile, rgyr_outfile, vmd_selection):
    vmd_code = VMD_MEASURE_SASA_AND_RGYR_SCRIPT%(pdb_file, "tmp_vmd_out_sasa", "tmp_vmd_out_rgyr", vmd_selection)
    execute_vmd_script(vmd_code)
    
    sasa = numpy.loadtxt("tmp_vmd_out_sasa")/100.
    rgyr = numpy.loadtxt("tmp_vmd_out_rgyr")
    
    #os.system("rm tmp_vmd_script tmp_vmd_out_sasa tmp_vmd_out_rgyr vmd_out")
    
    numpy.savetxt(sasa_outfile, sasa, "%.4f")
    numpy.savetxt(rgyr_outfile, rgyr, "%.4f")
    
# def calculate_torsional_variation(pdb_file, outfile):
    

def analyze_trajectory(traj_path, 
                       do_sasa, do_sasa_vmd, 
                       do_rgyr, do_rgyr_vmd, vmd_selection,
                       do_rmsf, 
                       report_pattern,
                       report_dir,
                       disp_logfile):
    
    trajectory = None
    data = None
    
    if do_sasa:
        print "Calculating SASA ..."
        if trajectory is None:
            trajectory = mdtraj.load(traj_path)
        #Shrake, A; Rupley, JA. (1973) J Mol Biol 79 (2): 351--71.
        sasa = mdtraj.shrake_rupley(trajectory, mode = 'residue').sum(axis=1)
        numpy.savetxt(traj_path+'.sasa', sasa)

    if do_sasa_vmd and do_rgyr_vmd:
        calculate_sasa_and_rgyr_with_vmd(traj_path, 
                                '%s.%s.sasa'%(traj_path, vmd_selection.replace(" ","_")),
                                '%s.%s.rgyr'%(traj_path, vmd_selection.replace(" ","_")),
                                vmd_selection)
    
    if do_sasa_vmd and not do_rgyr_vmd:
        calculate_sasa_with_vmd(traj_path, 
                                '%s.%s.sasa'%(traj_path, vmd_selection.replace(" ","_")), 
                                vmd_selection)
        
    if do_rgyr_vmd and not do_sasa_vmd:
        calculate_rgyr_with_vmd(traj_path, 
                                '%s.%s.rgyr'%(traj_path, vmd_selection.replace(" ","_")), 
                                vmd_selection)
    
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
#         data, _ = load_all_pdbs_ca([{
#                                        "source": traj_path, 
#                                        "base_selection":"name CA"   
#                                     }])
#         rmsf_array = rmsf(data.structure_ensemble)
        import prody
        
        rmsf_array = ca_rmsf(prody.proteins.pdbfile.parsePDB(traj_path, subset='ca'))
        numpy.savetxt(traj_path+'.rmsf', rmsf_array)

    if data is not None:
        del data
    
    if report_pattern != "":
        print "Extracting acceptance and energies from report files with pattern %s inside %s ..."%(report_pattern, 
                                                                                                    report_dir)
        
        files = glob.glob(os.path.join(report_dir, report_pattern))
        
        assert len(files)!=0, "No report file with pattern %s found inside %s"%(report_pattern,
                                                                                 report_dir)
        all_accepted = []
        all_total = []
        all_energies = []
        for report_file in files:
            total, accepted, energies = process_report_file(report_file)
            all_total.append(total)
            all_accepted.append(accepted)
            all_energies.append(list(energies))
        total = numpy.sum(all_total)
        accepted = numpy.sum(all_accepted)
        acceptance = accepted / total
        numpy.savetxt(traj_path+'.acc', [acceptance], fmt = "%.4f ")
        
        energy_handler = open(traj_path+'.ener', "w")
        for i in range(len(all_energies)):
            for j in range(len(all_energies[i])):
                energy_handler.write("%f\n"%all_energies[i][j])
            energy_handler.write("###\n")
            
    
    if disp_logfile != "":
        handler = open(disp_logfile)
        fractions = []
        for line in handler:
            if line[0:4] == "DBG:":
                if "iterations performed" in line:
                    parts = line.split()
                    fractions.append(float(parts[1])/float(parts[3]))
        numpy.savetxt(traj_path+'.frac', [numpy.mean(fractions)], fmt = "%.4f ")
        handler.close()

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
    parser.add_option("--sasa-vmd", dest="do_sasa_vmd", action="store_true", default=False)
    parser.add_option("--vmd-sel", dest="vmd_sel", default="all")
    parser.add_option("--rgyr", dest="do_rgyr", action="store_true", default=False)
    parser.add_option("--rgyr-vmd", dest="do_rgyr_vmd", action="store_true", default=False)
    parser.add_option("--rmsf", dest="do_rmsf", action="store_true", default=False)
    parser.add_option("--report", dest="report", default="")
    parser.add_option("--report-dir", dest="report_dir", default="")
    parser.add_option("--disp", dest="disp", default="")
    
    (options, args) = parser.parse_args()

    assert options.traj_path is not None, "It is mandatory to choose a valid trajectory file (-i option)"
    
    if  (options.do_sasa or options.do_rgyr) and not mdtraj_accessible: 
        print "It was not possible to load the 'mdtraj' module. Using the '--sasa' or '--rgyr' options is not possible."
        exit()
        
    if options.do_rmsf and not anmhelpers_accessible:
        print "It was not possible to load the 'anmichelpers' module. Using the '--rmsf' option is not possible."
        exit()
    
    _, ext = os.path.splitext(options.traj_path)
    if ext == ".txt":
        print "Reading from ", options.traj_path
        for line in open(options.traj_path).readlines():
            parts = line.rstrip('\r\n').split()
            print "Processing ", parts[0], " ..."
            if len(parts) == 2:
                analyze_trajectory(parts[0], 
                           options.do_sasa,
                           options.do_sasa_vmd,
                           options.do_rgyr,
                           options.do_rgyr_vmd,
                           options.vmd_sel,
                           options.do_rmsf,
                           parts[1],
                           parts[2],
                           options.disp)
            else:
                analyze_trajectory(parts[0], 
                           options.do_sasa,
                           options.do_sasa_vmd,
                           options.do_rgyr,
                           options.do_rgyr_vmd,
                           options.vmd_sel,
                           options.do_rmsf,
                           options.report,
                           options.report_dir,
                           options.disp)
    else:
        analyze_trajectory(options.traj_path, 
                           options.do_sasa,
                           options.do_sasa_vmd,
                           options.do_rgyr,
                           options.do_rgyr_vmd,
                           options.vmd_sel,
                           options.do_rmsf,
                           options.report,
                           options.report_dir,
                           options.disp)
