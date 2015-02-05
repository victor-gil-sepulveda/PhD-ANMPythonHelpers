import sys
import numpy as np
from functions import *

protein = sys.argv[1]

file_anm_cc = "files/" + protein + "_cc.slwevs"
file_imods_cc = "files/" + protein + "_imode_cart.evec"
protein_pdb = "files/" + protein + "_imode_model.pdb"

eigenvectors_anm_cc = parse_eigenvectors_cc(file_anm_cc)
eigenvectors_imods_cc, eigenvalues_imods_cc = parse_imods(file_imods_cc, 'cc')

# print len(eigenvectors_anm_cc[0]) # 214
# print len(eigenvectors_imods_cc[0]) # 1656

ca_positions = get_ca_positions(protein_pdb)

for mode1 in range(20):
  for mode2 in range(20):
    if mode1 == mode2:
      dot_products = []
      # contiene las tripletas correspondientes a las posiciones de los Ca
      eigenvectors_ca = get_ca_vectors(eigenvectors_imods_cc, mode1, ca_positions)

      for i in range(len(ca_positions)):
        a = np.array(eigenvectors_anm_cc[mode2][i])/norma(eigenvectors_anm_cc[mode2][i])
        b = np.array(eigenvectors_ca[i])/norma(eigenvectors_ca[i])
        dot = np.dot(a, b)
        dot_products.append(dot)

      # print dot_products
      mean = np.mean(dot_products)

      if mode1 == mode2:
        print mode1, mode2, mean
