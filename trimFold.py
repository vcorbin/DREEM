#!/usr/bin/python3
import os
import time
import sys
from glob import glob


#  TODO: take care of trim_size overflow； add trim_size argument to parallel, bic_helper
def trim_ref(mu, trim_size):
    trim_up, trim_down = trim_size
    with open(mu, 'r') as f:
        first_line = f.readline()
        ref_name = first_line.split('\t')[1]
        ref_name = ref_name.split('_')[0] + str(".txt")  # e.g. "HIVwtPhil.txt"
        start_pos = int(first_line.split("\t")[2])
        second_line = f.readline()
        length = int(second_line.split("\t")[2])
        end_pos = start_pos + length - 1
    f.close()
    reference = "/lab/solexa_rouskin/references/" + ref_name  # absolute path to the reference txt
    reference_size = os.stat(reference).st_size  # get the file size in bytes
    trimmed_reference_file_name = "trimmed_" + ref_name  # e.g. "trimmed_HIVwtPhil.txt"
    if reference_size < 1024 and all([size == 400 for size in trim_size]):  # No trimming if the reference < 1k AND no trimming option is given
        return mu, reference
    with open(reference, "r") as f:
        with open(trimmed_reference_file_name, "w") as f1:
            ref_string = ""
            for count, line in enumerate(f, -1):
                if count != -1:
                    ref_string += line
                else:
                    f1.write(line)
            trimed_start_pos = max(1, start_pos - trim_up)
            trimed_end_pos = min(len(ref_string), end_pos + trim_down)
            trimed_ref_string = ref_string[(trimed_start_pos - 1): trimed_end_pos]
            f1.write(trimed_ref_string)
    f.close()
    f1.close()

    # now rewrite the reference file and cluster file
    trimmed_cluster_file_name = "trimmed_up_" + str(trim_up) + "_down_" + str(trim_down) + "_" + mu.split("/")[-1]
    start = -1
    if trimed_start_pos == 1:  # trim_size too big, exceeding the reference start
        trim_up = start_pos - 1
    print("Starting position of trimming:" + str(trim_up))
    with open(mu, "r") as f:
        with open(trimmed_cluster_file_name, "w") as f1:
            for count, line in enumerate(f, 0):
                if count == 0:  # start
                    new_line = line.split('\t')
                    new_line[-1] = str(trim_up + 1) + '\n'
                    new_line = '\t'.join(new_line)
                    f1.write(new_line)
                elif line.split('\t')[0] == 'Position':
                    start = count + 1
                    f1.write(line)
                elif count == start:
                    new_line = line.split('\t')
                    new_line[0] = str(int(new_line[0]) - start_pos + trim_up + 1)  # Updated Position for cluster
                    new_line = '\t'.join(new_line)
                    start += 1
                    f1.write(new_line)
                else:
                    f1.write(line)
    f.close()
    f1.close()
    return trimmed_cluster_file_name, trimmed_reference_file_name


if __name__ == "__main__":
    mu = sys.argv[1]
    trim_size = int(sys.argv[2]), int(sys.argv[3])
    trim_ref(mu, trim_size)

    #     os.system("bsub -J c2c cluster_to_constraint.py " + mu)
    #     os.system("bsub -w 'ended(""c2c"")' runFoldDraw.py " + mu + " " + reference)
    # else:
    #     raise ValueError
