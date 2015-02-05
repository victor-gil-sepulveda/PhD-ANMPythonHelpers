import sys
import numpy as np
from functions import *
import matplotlib.pyplot as plt

protein = sys.argv[1]
# modo escogido
mode_anm = int(sys.argv[2])
mode_imods = int(sys.argv[3])

file_anm_cc = "files/" + protein + "_cc.slwevs"
file_imods_cc = "files/" + protein + "_imode_cart.evec"
protein_pdb = "files/" + protein + "_imode_model.pdb"

eigenvectors_anm_cc = parse_eigenvectors_cc(file_anm_cc)
eigenvectors_imods_cc, eigenvalues_imods_cc = parse_imods(file_imods_cc, 'cc')

# print len(eigenvectors_anm_cc[0]) # 214
# print len(eigenvectors_imods_cc[0]) # 1656

ca_positions = get_ca_positions(protein_pdb)

# contiene las tripletas correspondientes a las posiciones de los Ca
eigenvectors_ca = get_ca_vectors(eigenvectors_imods_cc, mode_imods, ca_positions)

norm_a = []
norm_b = []
dot_products = []
for i in range(len(ca_positions)):
    a = np.array(eigenvectors_anm_cc[mode_anm][i])/norma(eigenvectors_anm_cc[mode_anm][i])
    norm_a.append(norma(eigenvectors_anm_cc[mode_anm][i]))

    b = np.array(eigenvectors_ca[i])/norma(eigenvectors_ca[i])
    norm_b.append(norma(eigenvectors_ca[i]))
    dot_products.append(np.dot(a, b))

max_a = max(norm_a)
norm_a = [i/max_a for i in norm_a]

max_b = max(norm_b)
norm_b = [i/max_b for i in norm_b]

print np.mean(dot_products)

plt.plot(dot_products)
plt.plot(norm_a)
plt.plot(norm_b)
plt.ylabel('Dot products')
plt.xlabel('Residue number')
plt.plot(dot_products, color="blue", label="Dot products")
plt.plot(norm_a, color="red", label="Original vector norm")
plt.plot(norm_b, color="green", label="Transformed vector norm")
#plt.legend(loc='lower right', bbox_to_anchor=(1.1, 0))
# plt.show()
