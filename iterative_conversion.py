"""
Uses PELE++ to do the conversion. It has to run an specific test with specific input/output filenames.

Created on 03/04/2015

@author: vgil
"""
import os
import numpy

if __name__ == '__main__':
    
    scores = []
    for i in range(10):
        print "Iteration %d"%i

        # Get an initial guess again 
        # We use the template file as base and the new converted file to cover the
        # atoms undefined in the template (SD, BB ...)
        if i == 0:
            # This is the initial guess
            os.system("python /home/user/git/PhD-ANMPythonHelpers/convert_cc_pca_to_aa_cc_pca.py \
-i cc_pca_bb.nmd  -o cc_pca_aa_%d.nmd -t cc:pca -m PROPAGATE_CA -f ordered_atoms.txt"%i)
        else:
            os.system("python /home/user/git/PhD-ANMPythonHelpers/convert_cc_pca_to_aa_cc_pca.py \
-i cc_pca_bb.nmd  -a converted.nmd -o cc_pca_aa_%d.nmd \
-t cc:pca -m PROPAGATE_CA -f ordered_atoms.txt"%i)
        
        # Call PELE to do the conversion
        # PELE will produce "converted.nmd" file which is the result of converting "to convert.nmd"
        # to IC and then to CC (heavy) again
        os.system("cp cc_pca_aa_%d.vectors to_convert.vectors"%i)
        os.system("cp cc_pca_aa_%d.values to_convert.values"%i)
        os.system("~/workspace/ANMIC_AZ/Release/ANMIC_AZ")
    
        # Get the score for this conversion
        os.system("python /home/user/git/PhD-ANMPythonHelpers/plot_nmd_cc_angles_and_magnitudes.py \
--i1 ../CA/cc_pca_ca.nmd --i2 converted.nmd -o scores -f 0 -t 3 -s")
        scores.append(numpy.loadtxt("scores"))
        
    print scores
    numpy.savetxt("final_scores.txt",scores)