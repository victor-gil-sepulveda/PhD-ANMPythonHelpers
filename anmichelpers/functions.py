from prody import *
import math

def read_file_cc(path):
    file = open(path, "r")
    lines = []
    for line in file:
        line = line.split()
        lines.append(line[1:])
    file.close()

    return lines

def read_file_ci(path):
    file = open(path, "r")
    lines = []
    for line in file:
        line = line.split()
        lines.append(line)
    file.close()

    lines = lines[2:]

    return lines

    
def parse_imods(path, type):
    lines = read_file_ci(path)
    values = []
    vectors = []

    vector = []
    prev_line = '****'
    for l in range(len(lines)):
        line = lines[l]
        if l < len(lines) - 1:
            next_line = lines[l+1]
        if l > 0:
            prev_line = lines[l-1]

        if len(line) == 2 and prev_line[0] == '****':
            # inserta el eigenvalue
            values.append(float(line[1]))
        elif len(line) != 1:
            # transforma los valores a floats
            line = [float(x) for x in line]
            # inserta la linea actual al eigenvector
            vector += line

        elif line[0] == '****':
            # nuevo value, se inserta el eigevector completo a la lista
            if type == 'cc':
                vector = zip(*[iter(vector)]*3)
            vectors.append(vector)
            vector = []

    # inserta el ultimo eigenvector
    if type == 'cc':
        vector = zip(*[iter(vector)]*3)
    vectors.append(vector)

    return vectors[1:], values

def get_ca_positions(pdb):
    protein = parsePDB(pdb)
    hv = protein.getHierView()
    atoms = hv.getAtoms()

    indexes = []
    for i, a in enumerate(atoms.iterAtoms()):
        if a.getName() == 'CA':
            indexes.append(i)

    return indexes

def get_ca_vectors(vectors, mode, ca_positions):
    eigenvectors_ca = []
    for i in ca_positions:
        eigenvectors_ca.append(vectors[mode][i])
        # print i, vectors[mode][i]

    return eigenvectors_ca

def norma(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
