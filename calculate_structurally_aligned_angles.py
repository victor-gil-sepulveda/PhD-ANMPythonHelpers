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

def get_residue_angles(residue):
    try:
        phi = calcPhi(residue, radian=True)
    except:
        phi = 0
    
    try:
        psi = calcPsi(residue, radian=True)
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

def add_single_residue_angles(pdb, res_id, angles):
    residue_struct = pdb.select("resid %d"%res_id).copy()
    i = 0
    for residue in residue_struct.iterResidues():
        phi, psi = get_residue_angles(residue)
        add_angles(phi, psi, angles)
        i+=1
    assert i == 1, "[ERROR] resid %d selection returns more than one residues (%d)"%(res_id,i)
    return residue

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
    
    reference = parsePDB(options.reference)
    reference_structure = reference.select("resid "+" ".join([str(i) for i in range(options.from_res, options.to_res+1)]))
    writePDB("ref.pdb.chunk", reference_structure)
    #The reference chunk must NOT have RESIDUE GAPS
    assert options.to_res - options.from_res +1 == reference_structure.getHierView().numResidues(),\
    "[ERROR] There are gaps in the reference structure inside this residue range.%d %d"%(options.to_res - options.from_res +1, reference_structure.getHierView().numResidues())
    
    for pdb_path in open(options.db_list).readlines():
        sequence = []
        pdb_path = pdb_path.rstrip('\r\n')
        res_mapping_1 = MICAN().run(options.reference, pdb_path)
        res_mapping_2 = TMAlign().run(options.reference, pdb_path)
        
        # We get the one that got more residues aligned
        (key_len, res_mapping) = max((len(res_mapping_1.keys()), res_mapping_1),
                                     (len(res_mapping_2.keys()), res_mapping_2))
        
        pdb = parsePDB(pdb_path)
        res_selection = []
        angles = []
        
        for r_id in range(options.from_res, options.to_res+1):
            try:
                res_selection.append(res_mapping[r_id])
                residue = get_residue(pdb, res_mapping[r_id])               
                resname = residue.getResname()
                # Get the real residue (not a copy)
                real_residue = pdb.getHierView().getResidue(residue.getChain(), residue.getResnum())
                residue.setResnum(r_id)
            except KeyError:
                print "Correspondence not found for residue %d of template in pdb %s"%(r_id, pdb_path)
                # Only reason to not to have a residue mapped is that it is not present
                # in the target pdb file (as it is mandatory that there are  no gaps
                # in the reference structure).
                resname = "GAP"
                # angles.extend([inf,inf])
                
            sequence.append(resname)
        sequences.append(sequence)
        res_sel_string ="resid "+" ".join( [str(r_id) for r_id in res_selection])
        selected_residues_structure = pdb.select(res_sel_string)
        writePDB(pdb_path+".chunk",selected_residues_structure)
        all_angles.append(angles)
        
    # Save the angles    
#    angle_modes = stats.mode(angles)
#    for i in range(len(sequences)):
#        sequence = sequences[i]
#        for j in range(len(sequence)):
#            res = sequence[j]
#            if res == "GAP":
#                all_angles[i][j*2] = angle_modes[j*2]
#                all_angles[i][(j*2)+1] = angle_modes[(j*2)+1]
#    open("%s.angles"%(options.output),"w").write("\n".join([",".join([str(a) for a in angles]) for angles in all_angles]))
    
    # Write the sequences
    
                
                 
