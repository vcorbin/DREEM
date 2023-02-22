#!/usr/bin/python3
import os
import time
import sys
import EM_classes as classes
from itertools import combinations


class ChDir(object):
    """
    Step into a directory temporarily.
    """

    def __init__(self, path):
        self.old_dir = os.getcwd()
        self.new_dir = path

    def __enter__(self):
        os.chdir(self.new_dir)

    def __exit__(self, *args):
        os.chdir(self.old_dir)


def get_mu_ref(mu, trim_size):
    trim_up, trim_down = trim_size
    with open(mu, 'r') as f:
        first_line = f.readline()
        ref_name = first_line.split('\t')[1]
        ref_name = ref_name.split('_')[0] + str(".txt")  # e.g. "HIVwtPhil.txt"
    f.close()
    reference = "/lab/solexa_rouskin/references/" + ref_name  # absolute path to the reference txt
    reference_size = os.stat(reference).st_size  # get the file size in bytes
    if reference_size < 1024 and all([size == 400 for size in trim_size]):  # TODO: threshold
        return mu, reference
    else:
        trimmed_reference_file_name = "trimmed_" + ref_name  # e.g. "trimmed_HIVwtPhil.txt"
        trimmed_cluster_file_name = "trimmed_up_" + str(trim_up) + "_down_" + str(trim_down) + "_" + mu.split("/")[-1]
        return trimmed_cluster_file_name, trimmed_reference_file_name


if __name__ == '__main__':
    mu = sys.argv[1]
    trim_size = int(sys.argv[2]), int(sys.argv[3])
    trimmed_mu, trimmed_ref = get_mu_ref(mu, trim_size)
    count = 0
    path = "/".join(mu.split('/')[:-1])
    if path == "":
        path = "./"
    chdir = ChDir(path)
    chdir.__enter__()
    rel_path_mu = mu.split("/")[-1]
    rel_path_trimmed_mu = trimmed_mu.split("/")[-1]
    # ---scatterClusters.py---
    cluster_group = classes.clusterGroup(rel_path_mu)
    nclusters = len(cluster_group.cluster_names)
    if nclusters == 2:
        os.system("bsub -J scatter scatterClusters.py " + rel_path_mu)
    elif nclusters >= 3:
        for combo in combinations(cluster_group.cluster_names, 2):
            c1, c2 = combo
            os.system("bsub -J scatter scatterClusters.py " + rel_path_mu + " -c1 " + c1 + " -c2 " + c2)
    # ---scatterClusters.py---
    os.system("bsub -J trim trimFold.py " + rel_path_mu + " " + str(trim_size[0]) + " " + str(trim_size[1]))
    while not os.path.isfile(rel_path_trimmed_mu) or not os.path.isfile(trimmed_ref):
        time.sleep(1)
        count += 1
    if os.path.isfile(rel_path_trimmed_mu) and os.path.isfile(trimmed_ref):
        print("waiting time: " + str(count))
        os.system("bsub -J c2c -w 'ended(""trim"")' cluster_to_constraint.py " + rel_path_trimmed_mu)
        os.system("bsub -w 'ended(""c2c"")' runFoldDraw.py " + rel_path_trimmed_mu + " " + trimmed_ref)
    chdir.__exit__()