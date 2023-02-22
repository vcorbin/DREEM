#!/usr/bin/python
import EM_lib as lib
import EM_files as files
import pandas
import os
import time
import sys
import numpy as np
from glob import glob

def check_for_files(filepath):
    if glob(filepath) == []:
        return False
    for filepath_object in glob(filepath):
        if not os.path.isfile(filepath_object):
            return False
    return True

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='From a cluster file, create 1D constraint files (1 file per cluster)')
    parser.add_argument('input_filename', type=str, help='File containing the clusters e.g.: clusters_Ref_0-421.tx')
    args = parser.parse_args()
    print (args)
    prefix = os.path.splitext(args.input_filename)[0]
    offset = files.load_header(args.input_filename)[2].get_coordinates()[1] + 1

    clusters = files.LoadClusters(args.input_filename)
    seq = files.load_header(args.input_filename)[2].get_sequence()

    #clust_df = pandas.read_csv(args.input_filename, header=0, index_col=0, sep="\t")
    constraints = lib.normalise_clusters(clusters)
    constraints_VARNA = lib.normalise_clusters(clusters)
    for i in range(len(constraints)):
        for j in range(len(constraints[i])):
            if seq[j] == 'G' or seq[j] == 'U' or seq[j] == 'T' or clusters[i][j] == -999:
                constraints[i][j] = -999
                constraints_VARNA[i][j] = 0

    for i in range(len(constraints)):
        out_filename = prefix + "_1Dconstraints_cluster_" + str(i + 1) + ".txt"
        out_filename_VARNA = prefix + "_1Dconstraints_cluster_for_VARNA" + str(i + 1) + ".txt"
        constraints_df = pandas.DataFrame(constraints[i])
        constraints_df.index = np.arange(offset, offset+len(constraints[i]))
        constraints_df.to_csv(out_filename, header=None, index=True, sep=" ")

        # additional outputs for use of VARNA, index dropped
        constraints_VARNA_df = pandas.DataFrame(constraints_VARNA[i])
        constraints_VARNA_df.index = np.arange(offset, offset + len(constraints[i]))
        add_zeros = pandas.DataFrame(np.zeros(offset - 1))
        constraints_VARNA_df = pandas.concat([add_zeros, constraints_VARNA_df])
        constraints_VARNA_df = constraints_VARNA_df.applymap(lambda s: float(format(s, '.6f')))
        constraints_VARNA_df.to_csv(out_filename_VARNA, header=None, index=False, sep=" ")

    constraints_filepath = prefix + "_1Dconstraints_cluster_*"
    print("Constraint files all exist or not: " + str(check_for_files(constraints_filepath)))
    while not check_for_files(constraints_filepath):
        time.sleep(1)
        print("Waiting for cluster_to_constraint to complete")
    if check_for_files(constraints_filepath):
        sys.exit()
