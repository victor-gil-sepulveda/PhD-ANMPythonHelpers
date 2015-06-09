"""
Created on 22/05/2015

@author: vgil
"""
from optparse import OptionParser
from prody.proteins.pdbfile import parsePDB, writePDB
from obtain_structurally_aligned_angles import get_reference_structure,\
    get_best_res_mapping
import numpy
from numpy.core.numeric import inf
import os

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--ref", dest="reference")
    parser.add_option("--db", dest="db_list")
    parser.add_option("-o", dest="output")
    parser.add_option("-f", dest="from_res", type="int")
    parser.add_option("-t", dest="to_res", type="int")
    parser.add_option("--chunks", dest="chunk_path", type="string")
    
    (options, args) = parser.parse_args()
    
    sequences = []
    all_ca_coords = []
    reference_structure = get_reference_structure(options.reference, 
                                        "name CA and resid "+" ".join([str(i) for i in range(options.from_res, options.to_res+1)]),
                                        options.to_res - options.from_res +1)
    
    for pdb_path in open(options.db_list).readlines():
        pdb_path = pdb_path.rstrip('\r\n')
        
        res_mapping = get_best_res_mapping(options.reference, pdb_path)
        pdb = parsePDB(pdb_path)

        num_res = options.to_res - options.from_res +1
        sequence = ["GAP"]*num_res
        res_selection = []
        coords = [inf,inf,inf]*num_res
        
        index = 0
        center = numpy.array([0., 0., 0.])
        for r_id in range(options.from_res, options.to_res+1):
            if r_id in res_mapping:
                t_id = res_mapping[r_id]
                res_selection.append(t_id)
                residue_ca = pdb.select("name CA and resid %d"%t_id).copy()
                residue_coordsets = residue_ca.getCoordsets()
                sequence [index] = str(residue_ca.getResnames()[0])
                assert residue_coordsets.shape[1] == 1 , \
                        "[ERROR] resid %d selection returns more than one residues (%d)"%(t_id,residue_coordsets.shape[1])
                coords_list = residue_coordsets[0][0].tolist()                                                                  
                offset = index*3
                coords[offset] = coords_list[0]
                coords[offset+1] = coords_list[1]
                coords[offset+2] = coords_list[2]
                center += coords_list
            else:
                print "Correspondence not found for residue %d of template in pdb %s"%(r_id, pdb_path)
            index = index + 1
        center = center / len(res_selection)
        # Center coordinates
        for i in range(num_res):
            offset = i*3
            if sequence[i]!= "GAP":
                coords[offset] -= center[0]
                coords[offset+1] -= center[1]
                coords[offset+2] -= center[2]
        
        sequences.append(sequence)
        all_ca_coords.append(coords)
        
        if options.chunk_path is not None:
            pdb_name = os.path.basename(pdb_path)
            # Extract the chunk
            res_sel_string ="resid "+" ".join( [str(residue_index) for residue_index in res_selection])
            selected_residues_structure = pdb.select(res_sel_string)
            writePDB(os.path.join(options.chunk_path, pdb_name+".chunk"), selected_residues_structure)
    
    all_ca_coords = numpy.array(all_ca_coords)
    
    numpy.savetxt("%s.coords"%(options.output), all_ca_coords, fmt = "%.5f")
    open("%s.seq"%(options.output),"w").write("\n".join([" ".join(sequence) for sequence in sequences]))
    