import sys
from functions import *

file_values = sys.argv[1]
file_vectors = sys.argv[2]
vectors = parse_eigenvectors_cc(file_vectors)
values = parse_eigenvalues_cc(file_values)
values = values[:20]

# print jp cc formato alba
for i in range(len(values)):
    print values[i],
    for j in range(len(vectors[i])):
       print vectors[i][j][0], vectors[i][j][1], vectors[i][j][2],
    print
