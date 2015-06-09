'''
Created on 14/05/2015

@author: user
'''
import numpy
from optparse import OptionParser
from prody.proteins.pdbfile import parsePDB, writePDB
from structalignment.mican import MICAN
from structalignment.tmalign import TMAlign
from numpy.core.numeric import inf
from anmichelpers.tools.prody_amends import calcPsi, calcPhi

def get_residue_angles(residue, use_degrees):
    try:
        phi = calcPhi(residue, radian= not use_degrees)
    except Exception, e:
        print "Adding phi = inf in residue %s:%d because: %s"%(residue.getResname(), residue.getResnum(), e)
        phi = inf
    try:
        psi = calcPsi(residue, radian = not use_degrees)
    except Exception, e:
        print "Adding psi = inf in residue %s:%d because: %s"%(residue.getResname(), residue.getResnum(), e)
        psi = inf
    return phi, psi

def get_residue(pdb, res_id):
    residue_struct = pdb.select("resid %d"%res_id).copy()
    i = 0
    for residue in residue_struct.iterResidues():
        i+=1
    assert i == 1, "[ERROR] resid %d selection returns more than one residues (%d)"%(res_id,i)
    return residue

def get_reference_structure(path, selection_str, expected_num_res):
    reference = parsePDB(path)
    reference_structure = reference.select(selection_str)
    writePDB("ref.pdb.chunk", reference_structure)
    #The reference chunk must NOT have RESIDUE GAPS
    assert expected_num_res == reference_structure.getHierView().numResidues(),\
    "[ERROR] There are gaps in the reference structure inside this residue range.%d %d"%(options.final_res - options.initial_res +1, reference_structure.getHierView().numResidues())

def get_best_res_mapping(reference_path, target_path):
    res_mapping_1 = MICAN().run(reference_path, target_path)
    res_mapping_2 = TMAlign().run(reference_path, target_path)
    
    # We get the one that got more residues aligned
    (_, res_mapping) = max((len(res_mapping_1.keys()), res_mapping_1),
                                 (len(res_mapping_2.keys()), res_mapping_2))
    return res_mapping

def get_means_for_non_inf(all_angles):
    """
    Calculates the arithmetic average of all the non inf values 
    """
    means = []
    for res_angles in all_angles.T:
        mean_candidates = []
        for angle in res_angles:
            if angle != inf:
                mean_candidates.append(angle)
        if mean_candidates == []:
            means.append(0.0)
        else:
            means.append(numpy.mean(mean_candidates))
    return means

def get_angles_and_sequence_for_pdb_with_mapping(pdb, res_mapping, pdb_path, options):
    """
    template pdb used to generate res_mapping has no gaps
    'res_mapping' must contain sequential ids. This means that [(2,3),(5,6)] is a correct map, 
    but [(2,3), (5,2)] is not.
    """
    
    num_res = options.final_res - options.initial_res +1
    
    sequence = ["GAP"]*num_res
    angles = [inf]*2*num_res
    used_residues = []
    for r_id in range(options.initial_res, options.final_res+1):
        if r_id in res_mapping:
            target_id = res_mapping[r_id]
            index = r_id - options.initial_res
            offset = index*2
            
            # Residue copy to get information and assert uniqueness
            residue = get_residue(pdb, target_id)               
            resname = residue.getResname()
    
            # Modification and query of the real residue
            real_residue = pdb.getHierView().getResidue(residue.getChid(), 
                                                        residue.getResnum())
            real_residue.setResnum(target_id)
            phi, psi = get_residue_angles(real_residue, options.use_degrees)
            angles[offset] = phi
            angles[offset+1] = psi
            sequence[index] = resname
            used_residues.append(target_id)
    
    return angles, sequence, used_residues


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--ref", dest="reference")
    parser.add_option("--db", dest="db_list")
    parser.add_option("-o", dest="output")
    parser.add_option("-f", dest="initial_res", type="int")
    parser.add_option("-t", dest="final_res", type="int")
    parser.add_option("--degrees", dest="use_degrees", action="store_true", default=False)
    parser.add_option("--geopca", dest="geopca_format", action="store_true", default=False)
    (options, args) = parser.parse_args()
    
    sequences = []
    all_angles = []
    
    reference_structure = get_reference_structure(options.reference, 
                                        "resid "+" ".join([str(i) for i in range(options.initial_res, options.final_res+1)]),
                                        options.final_res - options.initial_res +1)
    
    for pdb_path in open(options.db_list).readlines():
        sequence = []
        pdb_path = pdb_path.rstrip('\r\n')
        print pdb_path
        
        res_mapping = get_best_res_mapping(options.reference, pdb_path)
        
        pdb = parsePDB(pdb_path)
        angles, sequence, res_selection = get_angles_and_sequence_for_pdb_with_mapping(pdb, res_mapping, pdb_path, options)
        sequences.append(sequence)
        all_angles.append(angles)

        # Extract the chunk        
        res_sel_string ="resid "+" ".join( [str(r_id) for r_id in res_selection])
        selected_residues_structure = pdb.select(res_sel_string)
        writePDB(pdb_path+".chunk", selected_residues_structure)
    
    all_angles = numpy.array(all_angles)
    
#    # Eliminate first column and change inf by mean values    
#    angle_means = get_means_for_non_inf(numpy.array(all_angles)) # <-- must be calculated with circular stats!
#    all_new_angles = []
#    for i in range(1, len(all_angles.T)):
#        angles = all_angles.T[i]
#        new_angles = []
#        for j in range(len(angles)):
#            if angles[j] == inf:
#                new_angles.append( angle_means[i])
#            else:
#                new_angles.append( angles[j])
#        all_new_angles.append(new_angles)
#    
#    all_angles = numpy.array(all_new_angles).T
    
    # Write the sequences
    open("%s.seq"%(options.output),"w").write("\n".join([" ".join(sequence) for sequence in sequences]))
    
    if options.geopca_format:
        open("%s.angles"%(options.output),"w").write("\n".join([",".join([str(a) for a in angles]) for angles in numpy.array(all_angles).T]))
    else:
        
        numpy.savetxt("%s.angles"%(options.output), all_angles, fmt="%.4f")
    
    
    
                
                 
