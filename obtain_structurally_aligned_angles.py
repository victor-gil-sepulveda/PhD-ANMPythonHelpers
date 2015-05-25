'''
Created on 14/05/2015

@author: user
'''
from optparse import OptionParser
from prody.proteins.pdbfile import parsePDB, writePDB
from structalignment.mican import MICAN
from structalignment.tmalign import TMAlign
from prody.measure.measure import calcPhi, calcPsi
from scipy import stats
from numpy.core.numeric import inf
import numpy

def get_residue_angles(residue):
    try:
        phi = calcPhi(residue, radian=False)
    except:
        phi = 0
    
    try:
        psi = calcPsi(residue, radian=False)
    except:
        psi = 0
    
    return phi, psi

def add_angles(phi, psi, angles):
    angles.extend([phi,psi])

def get_residue(pdb, res_id):
    residue_struct = pdb.select("resid %d"%res_id).copy()
    i = 0
    for residue in residue_struct.iterResidues():
        i+=1
    assert i == 1, "[ERROR] resid %d selection returns more than one residues (%d)"%(res_id,i)
    return residue

def add_single_residue_angles(residue, angles):
    phi, psi = get_residue_angles(residue)
    add_angles(phi, psi, angles)

def get_reference_structure(path, selection_str, expected_num_res):
    reference = parsePDB(path)
    reference_structure = reference.select(selection_str)
    writePDB("ref.pdb.chunk", reference_structure)
    #The reference chunk must NOT have RESIDUE GAPS
    assert expected_num_res == reference_structure.getHierView().numResidues(),\
    "[ERROR] There are gaps in the reference structure inside this residue range.%d %d"%(options.to_res - options.from_res +1, reference_structure.getHierView().numResidues())

def get_best_res_mapping(reference_path, target_path):
    res_mapping_1 = MICAN().run(reference_path, target_path)
    res_mapping_2 = TMAlign().run(reference_path, target_path)
    
    # We get the one that got more residues aligned
    (_, res_mapping) = max((len(res_mapping_1.keys()), res_mapping_1),
                                 (len(res_mapping_2.keys()), res_mapping_2))
    return res_mapping

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--ref", dest="reference")
    parser.add_option("--db", dest="db_list")
    parser.add_option("-o", dest="output")
    parser.add_option("-f", dest="from_res", type="int")
    parser.add_option("-t", dest="to_res", type="int")
    (options, args) = parser.parse_args()
    
    sequences = []
    all_angles = []
    
    reference_structure = get_reference_structure(options.reference, 
                                        "resid "+" ".join([str(i) for i in range(options.from_res, options.to_res+1)]),
                                        options.to_res - options.from_res +1)
    
    for pdb_path in open(options.db_list).readlines():
        sequence = []
        pdb_path = pdb_path.rstrip('\r\n')
        
        res_mapping = get_best_res_mapping(options.reference, pdb_path)
        
        pdb = parsePDB(pdb_path)
        res_selection = []
        angles = []
        
        for r_id in range(options.from_res, options.to_res+1):
            try:
                residue = get_residue(pdb, res_mapping[r_id])               
                resname = residue.getResname()
                # Get the real residue (not a copy)
                real_residue = pdb.getHierView().getResidue(residue.getChid(), 
                                                            residue.getResnum())
                real_residue.setResnum(r_id)
                res_selection.append(r_id)
                add_single_residue_angles(real_residue, angles)

            except KeyError:
                print "Correspondence not found for residue %d of template in pdb %s"%(r_id, pdb_path)
                # Only reason to not to have a residue mapped is that it is not present
                # in the target pdb file (as it is mandatory that there are  no gaps
                # in the reference structure).
                resname = "GAP"
                angles.extend([inf,inf])
                
            sequence.append(resname)
        sequences.append(sequence)
        res_sel_string ="resid "+" ".join( [str(r_id) for r_id in res_selection])
        selected_residues_structure = pdb.select(res_sel_string)
        writePDB(pdb_path+".chunk", selected_residues_structure)
        all_angles.append(angles)
        
    # Save the angles    
    angle_modes = stats.mode(angles)
    for i in range(len(sequences)):
        sequence = sequences[i]
        for j in range(len(sequence)):
            res = sequence[j]
            if res == "GAP":
                all_angles[i][j*2] = angle_modes[j*2]
                all_angles[i][(j*2)+1] = angle_modes[(j*2)+1]
    
    open("%s.angles"%(options.output),"w").write("\n".join([",".join([str(a) for a in angles]) for angles in numpy.array(all_angles).T]))
    # Write the sequences
    open("%s.seq"%(options.output),"w").write("\n".join([" ".join(sequence) for sequence in sequences]))
    
                
                 
