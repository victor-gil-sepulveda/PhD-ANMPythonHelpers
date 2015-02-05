import numpy as np
import sys
import math
from functions import *
import matplotlib.pyplot as plt

file = sys.argv[1]
freqs_cc = []
values_cc = parse_eigenvalues_cc(file)
for value in values_cc:
    f = value / (2*math.pi)
    freqs_cc.append(f)

file = sys.argv[2]
freqs_ci = []
vectors_ci, values_ci = parse_ci(file)
for value in values_ci:
    f = value / (2*math.pi)
    freqs_ci.append(f)


plt.title('freqs_ci')
plt.bar(range(0, 20), freqs_ci)
plt.show()

plt.title('freqs_cc')
plt.bar(range(0, 20), freqs_cc[:20])
plt.show()

width = 0.4 # the width of the bars
ran1 = range(0, 20)
ran2 = [i+width for i in ran1]

fig, ax = plt.subplots()
rects1 = ax.bar(ran1, freqs_cc[:20], width, color='g')
rects2 = ax.bar(ran2, freqs_ci, width, color='r')

plt.title('freqs')
plt.show()
