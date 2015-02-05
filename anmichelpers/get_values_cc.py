import sys
from functions import *

file = sys.argv[1]
values = parse_eigenvalues_cc(file)

print values
