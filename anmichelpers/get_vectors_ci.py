import sys
from functions import *

file = sys.argv[1]
vectors, values = parse_imods(file, 'ci')

print vectors
print values