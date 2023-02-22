#!/usr/bin/python3

import pandas as pd
import numpy as np
import argparse

def get_cluster(c, cluster_no=1):
    with open(c, 'r') as f:
        first_line = f.readline()
        ref_name = first_line.split('\t')[1]
        ref_name = ref_name.split('_')[0]  # e.g. "HIVwtPhil"
        second_line = f.readline()
    f.close()
    start_pos = int(first_line.split("\t")[2])
    length = int(second_line.split("\t")[2])
    end_pos = start_pos + length - 1
    df = pd.read_csv(c, sep='\t', skiprows=3)
    col = 'Cluster_' + str(cluster_no)
    cluster = pd.Series(df[col])
    cluster.index = pd.RangeIndex(start_pos, end_pos +1 )
    return cluster, ref_name


def fill_combine(combined, ref_seq, fill = -999):
    start = combined.index[0]
    end = combined.index[-1]
    length = end-start + 1
    filled_combined = combined.reindex(range(start, end+1 ), fill_value = fill)
    ref_seq = ref_seq[start:end+1]
    return filled_combined, ref_seq, start, length


def get_ref_seq(ref_file):
    with open(ref_file, 'r') as f:
        f.readline()
        ref_seq = f.readline().strip()
    f.close()
    return ref_seq


def main(f1, f2, cluster_no_1=1, cluster_no_2=1):
    left, ref_name1 = get_cluster(f1, cluster_no_1)
    right, ref_name2 = get_cluster(f2, cluster_no_2)
    assert(ref_name1 == ref_name2), "Error: Clusters have different reference"
    ref_name = ref_name1
    ref_file = ref_name + ".txt"
    combined = left.combine_first(right)
    reference = "/lab/solexa_rouskin/references/" + ref_file
    ref_seq = get_ref_seq(reference)
    filled_combined, ref_seq, start, length = fill_combine(combined, ref_seq)
    filled_combined.name = "Cluster_1"
    out_filename = "clusters_" + ref_name + "_" + str(start) + "_%d-%d_cl%d_combine_cl%d" % (start, start + length, cluster_no_1, cluster_no_2)+".txt"
    combined_df = pd.DataFrame(filled_combined ,columns = ["Cluster_1"])
    combined_df.to_csv(out_filename,sep="\t")
    with open(out_filename,'w') as f:
        f.writelines("@start\t" + ref_name + "_" + str(start) + "\t" + str(start) + '\n')
        f.writelines("@length\t" + ref_name + "_" + str(start) + "\t" + str(length) + '\n')
        f.writelines("@seq\t" + ref_name + "_" + str(start) + "\t" + ref_seq + '\n')
        combined_df.to_csv(f, sep="\t")
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Combine 2 cluster_probs, taking the value of the first cluster for ovelapping region')

    parser.add_argument('f1',type=str,
        help='String. Cluster file containing the first cluster.')

    parser.add_argument('-f2', '--input_filename_2', type=str, dest='f2', default=None,
                        help='String. Cluster file containing the second cluster.')

    parser.add_argument('-c1', '--cluster_no_1', type=int, dest='cluster_no_1', default=1,
                        help='String. Number of the cluster of interest in f1 header. If None, will use 1.')

    parser.add_argument('-c2', '--cluster_no_2', type=int, dest='cluster_no_2', default=1,
                        help='String. Number of the cluster of interest in f1 header. If None, will use 1.')

    args = parser.parse_args()

    main(args.f1, f2=args.f2, cluster_no_1=int(args.cluster_no_1), cluster_no_2=int(args.cluster_no_2))