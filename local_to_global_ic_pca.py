"""
Created on 25/05/2015

@author: vgil
"""
from optparse import OptionParser
from anmichelpers.parsers.pronmd import ProdyNMDParser
import numpy
from anmichelpers.writers.pronmd import ProdyNMDWriter

def fill_with_zeros(evecs, from_res, to_res, width):
    """
    In the IC case we have an initial residue with only one angle
    """
    real_from_res = from_res - 1
    real_to_res = to_res - 1
    real_width = (width*2)-1
     
    ending_left_pad_index = (real_from_res*2)-1
    start_right_pad_index = (real_to_res*2)+1
    left_padding = [0.]*ending_left_pad_index
    right_padding = [0.]*(real_width - start_right_pad_index)
    new_evecs = []
    for evec in evecs:
        new_evec = []
        new_evec.extend(left_padding)
        new_evec.extend(evec)
        new_evec.extend(right_padding)
        new_evecs.append(new_evec)
    return numpy.array(new_evecs)

def extract_evecs(evecs, from_res, to_res, sequence):
    """
    from and to are included
    """
    num_residues = to_res - from_res +1
    real_from_res = from_res - 1
    start_index = real_from_res * 2
    end_index = start_index + (num_residues*2) -1 
    # <- we neglect last phi angle because it
    # belongs to the other residue
    
    # we also neglect one phi angle per proline
    proline_phi_indices = []
    for i, r_name in enumerate(sequence):
        if i!=0 and r_name == "PRO" and i%2 == 1: 
            # first position does not have phi
            # and we have phis at all even positions
            proline_phi_indices.append(2*i-1)
            
    new_evecs = []
    for evec in evecs:
        new_evecs.append(numpy.delete(evec[start_index:end_index], proline_phi_indices))
    
    return numpy.array(new_evecs)

def extract_coords(coords, from_res, to_res):
    return coords[(from_res-1)*3:to_res*3]

def extract_resseq(seq, from_res, to_res):
    return extract_evecs([seq], from_res, to_res)[0]

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-g", dest="global_file")
    parser.add_option("-l", dest="local_file")
    parser.add_option("-f", dest="from_res", type = "int")
    parser.add_option("-t", dest="to_res", type = "int")
    parser.add_option("-m", dest="mode", type = "string")
    parser.add_option("-o", dest="output")
    (options, args) = parser.parse_args()

    # from_res starts by 1 (first residue is 1)
    assert options.from_res > 0 and options.to_res > 0 and options.from_res < options.to_res
    assert options.global_file is not None and options.local_file is not None and options.output is not None and options.mode is not None
    assert options.mode in ["LOCAL_TO_GLOBAL",
                            "GLOBAL_TO_LOCAL"]
    
    g_evals, g_evecs, g_header = ProdyNMDParser.read(options.global_file)
    l_evals, l_evecs, l_header = ProdyNMDParser.read(options.local_file)
    
    if "LOCAL_TO_GLOBAL" == options.mode :
        # we do not have all the information, that's why we have to rely in the
        # global data
        header = g_header
        evals = l_evals
        evecs = fill_with_zeros(l_evecs, 
                                options.from_res, 
                                options.to_res,
                                len(g_evecs[0])/3)
    else: #  "GLOBAL_TO_LOCAL"
        header = {}
        
        assert "resnames" in g_header, "The global file must define the 'resnames' section"
        
        header["resnames"] = extract_resseq(g_header["resnames"] , 
                                            options.from_res, 
                                            options.to_res)
        if "coordinates" in g_header:
            header["coordinates"] = extract_coords(g_header["coordinates"] , 
                                                options.from_res, 
                                                options.to_res)
        
        evals = g_evals
        evecs = extract_evecs(g_evecs, 
                              options.from_res, 
                              options.to_res,
                              header["resnames"]) # The extracted sequence
        
    header["title"] = options.output
    ProdyNMDWriter.write(options.output, evals, evecs, header)