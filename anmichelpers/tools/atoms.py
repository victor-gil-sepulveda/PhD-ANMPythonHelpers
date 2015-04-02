import numpy

class Atom:
    def __init__(self,name, resid, coords = [0.,0.,0.,], mode_v = []):
        self.name = name
        self.resid = resid
        self.coords = numpy.array(coords)
        self.mode_v = mode_v
    
    def __eq__(self, other):
        return (self.name == other.name) and (self.resid == other.resid)
    
    def __hash__(self):
        return hash((self.name, self.resid))

def atoms_from_header(header, eigenvectors):
    atoms = []
    
    for i in range(len(header["atomnames"])):
        if not "coordinates" in header:
            atoms.append(Atom(header["atomnames"][i],
                                    header["resids"][i]))
        else:
            offset = i*3
            atoms.append(Atom(header["atomnames"][i],
                                    header["resids"][i],
                                    coords = header["coordinates"][offset:offset+3])
                               )
    for mode in eigenvectors:
        mode_3t = numpy.resize(mode, (len(mode)/3,3))
        for i in range(len(atoms)):
            atoms[i].mode_v.append(mode_3t[i])
            
    return atoms

def get_CA_atoms( atoms ):
    ca_atoms = []
    for atom in atoms:
        if atom.name == "CA":
            ca_atoms.append(atom)
    return ca_atoms
