'''
Created on 14/05/2015

@author: user
'''
from optparse import OptionParser
from prody.proteins.pdbfile import parsePDB, writePDB
from structalignment.mican import MICAN
from structalignment.tmalign import TMAlign


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--ref", dest="reference")
    parser.add_option("--db", dest="db_list")
    parser.add_option("-o", dest="output")
    parser.add_option("-f", dest="from_res", type="int")
    parser.add_option("-t", dest="to_res", type="int")
    (options, args) = parser.parse_args()
    
    reference = parsePDB(options.reference)
    writePDB("ref.pdb.chunk", reference.select("resid "+" ".join([str(i) for i in range(options.from_res, options.to_res)])))
    
    for pdb_path in open(options.db_list).readlines():
        pdb_path = pdb_path.rstrip('\r\n')
        res_mapping_1 = MICAN().run(options.reference, pdb_path)
        res_mapping_2 = TMAlign().run(options.reference, pdb_path)
        
        # We get the one that got more residues aligned
        (key_len, res_mapping) = max((len(res_mapping_1.keys()), res_mapping_1),
                                     (len(res_mapping_2.keys()), res_mapping_2))
        
        pdb = parsePDB(pdb_path)
        res_selection = []
        for r_id in range(options.from_res, options.to_res+1):
            try:
                res_selection.append(res_mapping[r_id])
            except KeyError:
                print "SELSTR Correspondence not found for residue %d of template in pdb %s"%(r_id, pdb_path)

        res_sel_string ="resid "+" ".join( [str(r_id) for r_id in res_selection])
        print  "SELSTR", res_sel_string
        selected_residues_structure = pdb.select(res_sel_string)
        writePDB(pdb_path+".chunk",selected_residues_structure)
