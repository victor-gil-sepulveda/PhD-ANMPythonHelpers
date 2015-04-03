import numpy

class Atom:
    def __init__(self,name, resid, coords = [0.,0.,0.,], mode_v = None):
        self.name = name
        self.resid = resid
        self.coords = numpy.array(coords)
        if mode_v is None:
            self.mode_v = []
        else:
            self.mode_v = mode_v
    
    def __eq__(self, other):
        return (self.name == other.name) and (self.resid == other.resid)
    
    def __hash__(self):
        return hash((self.name, self.resid))

def atoms_from_header(header, eigenvectors):
    atoms = []
    coords_v3 = numpy.resize(header["coordinates"], (len(header["coordinates"])/3, 3))
    
    for i in range(len(header["atomnames"])):
        if not "coordinates" in header:
            atoms.append(Atom(header["atomnames"][i],
                                    header["resids"][i]))
        else:
            atoms.append(Atom(header["atomnames"][i],
                                    header["resids"][i],
                                    coords = coords_v3[i])
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

def get_CA_modes(header,eigenvectors):
    # Get only the CAs
    atoms = atoms_from_header(header, eigenvectors)
    ca_atoms = get_CA_atoms(atoms)
    # remount eigenvalues and header
    new_header = {
                  "atomnames":[],
                  "resids":[],
                  "coordinates":[]
                  }
    new_modes = []
    for ca in ca_atoms:
        new_header["atomnames"].append("CA")
        new_header["resids"].append(ca.resid)
        new_header["coordinates"].extend(ca.coords)
        for i in range(len(ca.mode_v)):
            try:
                new_modes[i].extend(ca.mode_v[i])
            except IndexError:
                new_modes.append(list(ca.mode_v[i]))
    
    new_evecs = numpy.array(new_modes)
    return new_header, new_evecs

