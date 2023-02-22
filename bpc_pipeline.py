#!/usr/bin/python
from split_responsibilitiesProb import split_reponsibilities
from os import system
import os
import numpy as np
import pandas as pd
from BasePairsCorrelation import getArray
from BasePairsCorrelation import run_bpc

def run_pipeline(inFile, args):
    run_bpc(inFile, args)
    arr = getArray(args.input_filename)
    cluster = np.nansum(arr, axis=0)/len(arr)
    inFile_prefix = os.path.splitext(inFile)[0]
    z_filename = inFile_prefix + "_Zscores.txt"
    fold_out = inFile_prefix+ "_fold_out.txt"
    expand = np.array(args.expand.split(","), dtype=int)
    dms_data = np.array([[i for i in range(expand[0] + 1, expand[0] + len(cluster) + 1)], cluster]).T.astype(str)
    dms_file = inFile_prefix+ "_dms.txt"
    f = open(dms_file, 'w')
    for row in dms_data:
        if float(row[1]) != 0:
            f.write(str(int(float(row[0]))) + " " + row[1] + "\n")
    f.close()
    os.system("Fold %s %s -x %s -dms %s" % (args.fasta, fold_out, z_filename, dms_file))
    outfile = inFile_prefix+ "_vienna.txt"
    os.system("ct2dot %s %s %s" % (fold_out, args.structure_number, outfile))

def boolean(str):
    if str == "True":
        return True
    return False

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Run a EM to cluster reads from a bam/sam file or a binary_vector '
                                                 'file')
    parser.add_argument('input_filename', type=str,
                        help='Preferably a binary_vector file (.txt) created by sam_to_binary(). If a sam/bam is given, '
                             'a binary_vector file of the same name will be created.')
    parser.add_argument('fasta', type=str, help='Name of fasta file')
    parser.add_argument('-s', '--sep', type=str, dest='sep', default='\t', help='Separation used in the bit_vector.')
    parser.add_argument('-o', '--open', type=boolean, dest='open', default=False, help='True if should not open heatmap')
    parser.add_argument('-e', '--expand', type=str, dest='expand', default='0,0',
                        help='2 integer separated by a comma. Number of 0s to put in front and behind the correlation '
                             'matrices.')
    parser.add_argument('-r', '--responsibilities', type=str, dest='resp', default=None, help='responsibilities filename '
                             'if preforming BPC for individual clusters')
    parser.add_argument('-v', '--vienna', type=int, dest='structure_number', default=1, help='The number, one-indexed, '
                             'of the structure to convert in the CT file (use -1 or ALL to convert all structures)')

    args = parser.parse_args()

    if args.resp:
        filenames = split_reponsibilities(args.resp)
    else:
        filenames = [args.input_filename]
    for f in filenames:
        run_pipeline(f, args)