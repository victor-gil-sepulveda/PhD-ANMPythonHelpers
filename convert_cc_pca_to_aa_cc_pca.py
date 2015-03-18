'''
Created on 18/03/2015

@author: vgil
'''
from anmichelpers.parsers.pronmd import ProdyNMDParser
import numpy
from optparse import OptionParser
from anmichelpers.writers.pronmd import ProdyNMDWriter
import os

class Atom:
    def __init__(self,name, resid, coords = [0.,0.,0.,], mode_v = [0.,0.,0.,]):
        self.name = name
        self.resid = resid
        self.coords = coords
        self.mode_v = mode_v
    
    def __eq__(self, other):
        return (self.name == other.name) and (self.resid == other.resid)
    
    def __hash__(self):
        return hash((self.name, self.resid))

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-o", dest="output")
    parser.add_option("-t", dest="mode_type")
    parser.add_option("-m", dest="filling_method")
    parser.add_option("-f", dest="atom_order_file")
    (options, args) = parser.parse_args()

    filling_methods = ["ZEROS","PROPAGATE_CA"]
    assert options.filling_method in filling_methods,\
     "[ERROR] The filling method can only be one of these %s"%(str(filling_methods))

    # Parse ordered atoms
    handler = open(options.atom_order_file)
    ordered_atoms = []
    lines = handler.readlines()
    for line in lines:
        name, resid_s, x_s, y_s, z_s  = line.split()
        ordered_atoms.append(Atom(name, int(resid_s), 
                                  [float(x_s),
                                   float(y_s),
                                   float(z_s)]))
    
    # Parse the nmd file
    input_atoms = []
    eigenvalues, eigenvectors, header = ProdyNMDParser.read(options.input)
    assert "resids" in header, '[ERROR] The "resids" tag is not present into the nmd input file'
    assert "atomnames" in header, '[ERROR] The "atomnames" tag is not present into the nmd input file'
    assert len(eigenvectors[0]) == len(header["resids"])*3, '[ERROR] "mode" and "resids" do not have compatible sizes'
    assert len(eigenvectors[0]) == len(header["atomnames"])*3, '[ERROR] "mode" and "atomnames" do not have compatible sizes'
#    coords_3t = numpy.resize(header["coordinates"], (len(header["coordinates"])/3, 3)
    
    for i in range(len(header["resids"])):
        input_atoms.append(Atom(header["atomnames"][i],
                                header["resids"][i]))
    
    # Create data for the new file
    name, extension = os.path.splitext(options.output)
    new_header = {
                  "name":name,
                  "type":options.mode_type,
                  "coordinates":[],
                  "resids":[],
                  "atomnames":[]
                  }
    
    # New header
    for atom in ordered_atoms:
        new_header["coordinates"].extend(atom.coords)
        new_header["atomnames"].append(atom.name)
        new_header["resids"].append(atom.resid)
    
    # New modes
    backbone_atoms = ["N","O","C","CA","OXT"]
    new_eigenvectors = []
    for mode in eigenvectors:
        # Get the mode correspondence with input atoms
        evect_3t = numpy.resize(mode, (len(mode)/3,3))
        mode_map = {}
        for i,atom in enumerate(input_atoms):
            mode_map[input_atoms[i]] = evect_3t[i]

        # Obtain new modes for the ordered atoms array
        new_mode = []
        for atom in ordered_atoms:
            try:
                old_atom_mode = mode_map[atom]
                new_mode.extend(old_atom_mode)
            except:
                print "Atom not found: ", atom.name, atom.resid
                # If not found we have to decide how to fill it    
                if options.filling_method == "PROPAGATE_CA":
                    ca_mode_v = mode_map[Atom("CA",atom.resid)].mode_v
                    new_mode.extend(ca_mode_v)
                elif options.filling_method == "ZEROS":
                    new_mode.extend([0.,0.,0.])
        
        new_eigenvectors.append(new_mode)

    # Write it
    ProdyNMDWriter.write(name, eigenvalues, new_eigenvectors, new_header)