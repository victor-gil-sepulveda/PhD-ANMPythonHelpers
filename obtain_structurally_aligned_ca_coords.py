"""
Created on 22/05/2015

@author: vgil
"""
from optparse import OptionParser
from prody.proteins.pdbfile import parsePDB
from obtain_structurally_aligned_angles import get_reference_structure,\
    get_best_res_mapping
import numpy

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--ref", dest="reference")
    parser.add_option("--db", dest="db_list")
    parser.add_option("-o", dest="output")
    parser.add_option("-f", dest="from_res", type="int")
    parser.add_option("-t", dest="to_res", type="int")
    (options, args) = parser.parse_args()
    
    sequences = []
    all_ca_coords = []
    
    reference_structure = get_reference_structure(options.reference, 
                                        "name CA and resid "+" ".join([str(i) for i in range(options.from_res, options.to_res+1)]),
                                        options.to_res - options.from_res +1)
    
    for pdb_path in open(options.db_list).readlines():
        sequence = []
        pdb_path = pdb_path.rstrip('\r\n')
        
        res_mapping = get_best_res_mapping(options.reference, pdb_path)
        
        pdb = parsePDB(pdb_path)
        res_selection = []
        coords = []
        for r_id in range(options.from_res, options.to_res+1):
            try:
                res_id = res_mapping[r_id]
                res_selection.append(res_id)
                
                residue_ca = pdb.select("name CA and resid %d"%res_id).copy()
                residue_coordsets = residue_ca.getCoordsets()
                resname = str(residue_ca.getResnames()[0])
                assert residue_coordsets.shape[1] == 1 , \
                        "[ERROR] resid %d selection returns more than one residues (%d)"%(res_id,
                                                                                  residue_coordsets.shape[1])
                coords.extend(residue_coordsets[0][0].tolist())
            except KeyError:
                print "Correspondence not found for residue %d of template in pdb %s"%(r_id, pdb_path)
                # Only reason to not to have a residue mapped is that it is not present
                # in the target pdb file (as it is mandatory that there are  no gaps
                # in the reference structure).
                resname = "GAP"
                coords.extend([0,0,0])
                
            sequence.append(resname)
        sequences.append(sequence)
        all_ca_coords.append(coords)
    
    numpy.savetxt("%s.ca_coords"%(options.output),all_ca_coords)
    open("%s.seq"%(options.output),"w").write("\n".join([" ".join(sequence) for sequence in sequences]))
    