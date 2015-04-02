'''
Created on 18/03/2015

@author: vgil
'''
from anmichelpers.parsers.pronmd import ProdyNMDParser
import numpy
from optparse import OptionParser
from anmichelpers.writers.pronmd import ProdyNMDWriter
import os
import math

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-a", dest="alternate_input")
    parser.add_option("-o", dest="output")
    parser.add_option("-t", dest="mode_type")
    parser.add_option("-m", dest="filling_method")
    parser.add_option("-f", dest="atom_order_file")
    (options, args) = parser.parse_args()

    filling_methods = ["ZEROS","PROPAGATE_CA", "PROPAGATE_CA_TURN"]
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
    
    for i in range(len(header["resids"])):
        if not "coordinates" in header:
            input_atoms.append(Atom(header["atomnames"][i],
                                    header["resids"][i]))
        else:
            offset = i*3
            input_atoms.append(Atom(header["atomnames"][i],
                                    header["resids"][i],
                                    coords = header["coordinates"][offset:offset+3])
                               )

    """
    Load the alternate input file. If an atom is not present in the input file, we can still
    look for it here.
    """
    alt_eigenvectors = None
    alt_atoms = None
    if options.alternate_input:
        _, alt_eigenvectors, alt_header = ProdyNMDParser.read(options.alternate_input)
        print alt_header.keys()
        alt_atoms = []
        for i in range(len(alt_header["resids"])):
            alt_atoms.append(Atom(alt_header["atomnames"][i],
                                    alt_header["resids"][i]))
    
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
    for i,mode in enumerate(eigenvectors):
        # Get the mode correspondence with input atoms
        mode_map = {}
        evect_3t = numpy.resize(mode, (len(mode)/3,3))
        coords_map = {}
        for j in range(len((input_atoms))):
            mode_map[input_atoms[j]] = evect_3t[j]
            coords_map[input_atoms[j]] = input_atoms[j].coords

        if alt_atoms is not None:
            alt_mode_map = {}
            alt_evect_3t = numpy.resize(alt_eigenvectors[i], (len(alt_eigenvectors[i])/3,3))
            for j in range(len((alt_atoms))):
                alt_mode_map[alt_atoms[j]] = alt_evect_3t[j]
            
        # Obtain new modes for the ordered atoms array
        new_mode = []
        for atom in ordered_atoms:
            try:
                old_atom_mode = mode_map[atom]
                new_mode.extend(old_atom_mode)
            except:
                try:
                    old_atom_mode = alt_mode_map[atom]
                    new_mode.extend(old_atom_mode)
                except:
                    print "Not found ", atom.name, atom.resid
                    # If not found we have to decide how to fill it    
                    if options.filling_method == "PROPAGATE_CA":
                        """
                        The variance of the atoms with unknown pca values
                        is the same as the CA atom (they follow a linear 
                        diferential displacement)
                        """
                        ca_mode_v = mode_map[Atom("CA",atom.resid)]
                        new_mode.extend(ca_mode_v)
                    elif options.filling_method == "ZEROS":
                        """
                        Unknown atoms will be filled with [0,0,0], indicating
                        that there is no variance for them.
                        """
                        new_mode.extend([0.,0.,0.])
                    elif options.filling_method == "PROPAGATE_CA_TURN":
                        """
                        Unknown BB atoms have the same axis of variance than
                        CA atoms, but as the linear difference will come from a 
                        rotation of the backbone torsions, the sidechain variation
                        will take this into account.
                        """ 
                        ca_atom = Atom("CA",atom.resid)
                        ca_mode_v = mode_map[ca_atom]
                        if atom.name in backbone_atoms:
                            new_mode.extend(ca_mode_v)
                        else:
                            n_atom = Atom("N", atom.resid)
                            #look for the atoms into the "ordered atoms" vector
                            ca_atom_o =  find(lambda a: a == ca_atom, ordered_atoms)
                            n_atom_o =  find(lambda a: a == n_atom, ordered_atoms)
                            # vector from N to CA 
                            n_c = n_atom_o.coords - ca_atom_o.coords
                            # vector from CA to Atom
                            c_a = atom.coords - ca_atom_o.coords
                            # cross
                            n_cxmode = numpy.cross( n_c, ca_mode_v)
                            # cross again
                            c_axn_cxmode = numpy.cross(c_a, n_cxmode)
                            # module is CA mode v module
                            new_mode.extend( c_axn_cxmode  * ( math.sqrt(numpy.dot(ca_mode_v, ca_mode_v)) / math.sqrt(numpy.dot(c_axn_cxmode,c_axn_cxmode))))
        
        new_eigenvectors.append(new_mode)
    new_eigenvectors = numpy.array(new_eigenvectors)
    
    # Write it
    ProdyNMDWriter.write(name, eigenvalues, new_eigenvectors, new_header)
    numpy.savetxt(name+".values", eigenvalues)
    numpy.savetxt(name+".vectors", new_eigenvectors)
    