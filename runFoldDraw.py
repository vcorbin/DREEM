#!/usr/bin/python3
import os
import sys
from glob import glob

if __name__ == "__main__":
    mu = sys.argv[1]
    reference = sys.argv[2]

    cluster_name = mu.split('/')[-1].split('.')[0]  # e.g. "clusters_HIVwtPhil_480_480-760"
    constraints = glob(mu.split('.')[0] + '_1Dconstraints_cluster_[0-9].txt')  # [0-9] to avoid selecting constraints for VARNA
    constraints.sort()  # constraint_1, constraint_2 ...
    id_constraint = 0
    for constraint in constraints:
        id_constraint += 1
        # fold_out_name = constraint + " " + cluster_name + "_" + str(id_constraint)
        fold_out_name = cluster_name + "_" + str(id_constraint)
        viennna_out_name = cluster_name + "_" + str(id_constraint) + "_vienna.txt"
        draw_out_name = cluster_name + "Draw_" + str(id_constraint) + ".ps"
        if id_constraint == 1:  # No draw yet
            wait_option = ""
        else:
            wait_option = "-w 'ended(""structure"")'"
        os.system("bsub -J folding " + wait_option + " Fold " + reference + " " + fold_out_name + " -dms " + constraint)
        os.system("bsub -J vienna -w 'ended(""folding"")' ct2dot " + fold_out_name + " ALL " + viennna_out_name)
        os.system("bsub -J drawing -w 'ended(""vienna"")' draw " + fold_out_name + " " + draw_out_name + " -s " + constraint)
        os.system("bsub -J structure -w 'ended(""folding"")' test_structure.py " + mu + " " + viennna_out_name + " " + str(id_constraint))