import sys
from functions import *

file = sys.argv[1]
vectors, values = parse_imods(file, 'cc')

# print jp imods formato alba
for i in range(len(values)):
    print values[i],
    for j in range(len(vectors[i])):
       print vectors[i][j][0], vectors[i][j][1], vectors[i][j][2],
    print
