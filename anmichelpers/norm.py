import numpy as np
from functions import *

vo1 = [0.014375, -0.022887, 0.013381]
vo1 = np.array(vo1)/norma(vo1)

vo2 = [0.012958, -0.020079, 0.012157]
vo2 = np.array(vo2)/norma(vo2)

vo3 = [0.004634, -0.012146, 0.013471]
vo3 = np.array(vo3)/norma(vo3)

print

vt1 = [0.001764, -0.002802, 0.001654]
vt1 = np.array(vt1)/norma(vt1)

vt2 = [0.001519, -0.002451, 0.001564]
vt2 = np.array(vt2)/norma(vt2)

vt3 = [0.001404, -0.002416, 0.001632]
vt3 = np.array(vt3)/norma(vt3)

print np.dot(vo1, vt1)
print np.dot(vo2, vt2)
print np.dot(vo3, vt3)