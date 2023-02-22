#!/usr/bin/python3

import sys
import EM_files as files
import os
from Clustering import post_process

# What it does: Breaks up each EM run into it own bsub -n 9 job.
# How to run it:
# Step 1:
# python parallel.py <bitvector filename> <all the other options you would use on Clustering.py except -s >
# Example: python parallel.py 170210Rou_D17-2031HIVwt_480_780.bam_bitvector.txt -r 1000 -n 2 -i 10

# Step 2: to recombine the runs afterward (that is, to put the d1 with the others and to collect the LLs and BIC scores), run this:
# python parallel.py combine <dir_name>
# (for example: python parallel.py combine 170210Rou_D17-2031HIVwt_480_780.bam_bitvector_Nclusters_2_Nruns_2_maxits_10)
# Make sure the dir_name is the directory name of the MULTIPLE clusters, not the nclusters_1 directory.
# Note: you must wait until all the runs are done (count your lsf emails lol or check the directory) until you can run the combine command.

if __name__ == "__main__":

    if sys.argv[1] == "combine":
        workDir = sys.argv[2]
        sub_dir_prefix = sys.argv[3]
        EM_its = int(sys.argv[4])
        deletions = sys.argv[5]
        is_d1 = sys.argv[6]
        trim_size = int(sys.argv[7]), int(sys.argv[8])
        post_process(workDir, sub_dir_prefix, EM_its, deletions, is_d1, trim_size)
        sys.exit()

    print(" Running Exact_clustering in Parallel")

    args, workDir, prefix, seeds = files.set_up(convert_args=False)

    if not args.threading:
        threading = '--no_threading'
    else:
        threading = ''

    if args.remove_jackpotting:
        remove_jackpotting = '-j'
    else:
        remove_jackpotting = ''

    sub_dir_prefix = "run_"
    # run d1
    job_dir_prefix = workDir + "/" + sub_dir_prefix
    job_dir = job_dir_prefix + "d1"
    print(job_dir)

    os.system('bsub -n 20 -J ' + job_dir + ' -q ' + args.queue +
              ' Clustering.py %s -d %s -n %s -i %s -r %s -m %s -M %s -I %s -q %s -P %f -s %s -W %s -del %s -dmuts %s %s %s'
              % (args.input_filename, 1, 1, args.max_its, args.nreads, args.muts_threshold,
                 args.mask_threshold, args.mask_index, args.na_threshold, args.max_mut_prob, args.seed, job_dir, args.deletions, args.dmuts,
                 threading, remove_jackpotting))

    command = "bsub -q " + args.queue + " -w 'ended(\"" + job_dir + "\")' -J combine " + "parallel.py combine %s %s %d %s %s %d %d" % (workDir,
                                                                                                  sub_dir_prefix,
                                                                                                  args.EM_its, args.deletions, "yes", args.trimUp, args.trimDown)
    os.system(command)
    # if correct_for_duplicate is on, need to wait for d1 to be finish to use its result for the correction.
    # TO DO
    # run the other runs
    seed = args.seed
    for i in range(args.EM_its):
        job_dir = job_dir_prefix + str(i + 1)
        os.system('bsub -q ' + args.queue + ' -n 12 -J ' + job_dir +
                  ' Clustering.py %s -d %s -n %s -i %s -r %s -m %s -M %s -I %s -q %s -P %f -s %s -W %s -del %s -dmuts %s %s %s'
                  % (args.input_filename, args.ncluster, 1, args.max_its, args.nreads, args.muts_threshold,
                     args.mask_threshold, args.mask_index, args.na_threshold, args.max_mut_prob, seeds[i],
                     job_dir, args.deletions, args.dmuts, threading, remove_jackpotting))

    command = "bsub -q " + args.queue + " -w 'ended(\"" + job_dir_prefix + "*\")' -J combine " + "parallel.py combine %s %s %d %s %s %d %d" % (workDir,
                                                                                                  sub_dir_prefix,
                                                                                                  args.EM_its, args.deletions, "no", args.trimUp, args.trimDown)
    print(command)
    os.system(command)

    print("Done", workDir)
    if args.recur == "yes":
        if args.ncluster <= 2:
            is_d1 = 1
        elif args.ncluster > 2:
            is_d1 = 0

        # bic_helper
        command = "bsub -w 'ended(""combine"")'" + \
                  " bic_helper.py %s %d %d %f %d %f %s %s %f %s %s %d %d %s %d %d %d %s" % (args.input_filename,
                                                                                args.ncluster, args.max_its, args.nreads, args.muts_threshold,
                                                                                args.mask_threshold, args.mask_index,
                                                                                args.na_threshold, args.max_mut_prob, args.deletions,
                                                                                workDir, args.EM_its, is_d1, args.queue, args.trimUp, args.trimDown, args.dmuts, remove_jackpotting)
        print("Current working directory： " + workDir)
        os.system(command)