#!/usr/bin/python3
import EM_files as files
import sys
import os


def update_bic(input_filename, ncluster, max_its, nreads, muts_threshold, mask_threshold, mask_index,
               na_threshold, max_mut_prob, deletions, workDir, EM_its, is_d1, queue, trim_size, remove_jackpotting, dmuts=3):
    """
    Additional parameters:
    :param workDir:
    :param EM_its:
    :param is_d1: if this is the first run where d1 is considered
    :param last_bic:
    """
    sub_dir_prefix = "run_"
    job_dir_prefix = workDir + "/" + sub_dir_prefix
    prefix = os.path.splitext(os.path.basename(input_filename))[0]
    # read the bics from loglikelihood.txt
    best_run = files.clean_up(workDir, sub_dir_prefix, EM_its)
    loglik_file = workDir + '/loglikelihoods.txt'
    with open(loglik_file, 'r') as f:
        for count, line in enumerate(f, -1):
            if is_d1 == 1:
                if count == 0:
                    last_bic = line.split('\t')[-2]  # for the first task with, bic_last is bic of d1
                    # last_bic = 1e8
                    # this is only for debugging as with small num_iter the d1 performs the best
            try:
                if count == int(best_run):
                    cur_bic = line.split('\t')[-2]
            except ValueError:
                print("d1 is the best run!")
                sys.exit()
    f.close()
    if is_d1 == 0:  # get the last best bic
        prev_ncluster = ncluster - 1
        prev_workDir = files.get_dir('.', prev_ncluster, EM_its, max_its,
                                prefix=prefix)
        prev_best_run = files.clean_up(prev_workDir, sub_dir_prefix, EM_its)
        prev_loglik_file = prev_workDir + '/loglikelihoods.txt'
        print("previous working_dir: " + prev_workDir)
        with open(prev_loglik_file, 'r') as f:
            for count, line in enumerate(f, -1):
                if count == int(prev_best_run):
                    last_bic = line.split('\t')[-2]
        f.close()

    print("Is this the first run?: " + str(is_d1))
    print("last bic: " + str(last_bic))
    print("current bic: " + str(cur_bic))

    if float(last_bic) <= float(cur_bic):
        print("bic stops improving. Stop with optimal number of cluster: " + str(ncluster-1))
        sys.exit()
    else:
        # run next iteration
        print("bic is improving. Continue parallel.")
        ncluster += 1
        workDir = files.get_dir('.', ncluster, EM_its, max_its,
                                prefix=prefix)  # update workDir with ncluster += 1
        job_dir_prefix = workDir + "/" + sub_dir_prefix
        print("Updated working directory: " + workDir)
        print("Updated number of clusters: " + str(ncluster))
        # run parallel with ncluster += 1
        os.system('bsub -n 12 -J ' + workDir +
                  " parallel.py %s -d %d -n %d -i %d -r %f -m %s -M %s -I %s -q %s -P %f -del %s -W %s -recur yes -queue %s -trimUp %d -trimDown %d -dmuts %s %s"
                  % (input_filename, ncluster, EM_its, max_its, nreads, muts_threshold,
                     mask_threshold, mask_index, na_threshold, max_mut_prob,
                     deletions, workDir, queue, trim_size[0], trim_size[1], dmuts, remove_jackpotting))


if __name__ == "__main__":
    # TODO: pass the args between programs in command line
    input_filename = sys.argv[1]
    ncluster = int(sys.argv[2])
    max_its = int(sys.argv[3])
    nreads = float(sys.argv[4])
    muts_threshold = int(sys.argv[5])
    mask_threshold = float(sys.argv[6])
    mask_index = sys.argv[7]
    na_threshold = float(sys.argv[8])
    max_mut_prob = float(sys.argv[9])
    deletions = sys.argv[10]
    workDir = sys.argv[11]
    EM_its = int(sys.argv[12])
    is_d1 = int(sys.argv[13])
    queue = sys.argv[14]
    trim_size = int(sys.argv[15]), int(sys.argv[16])
    dmuts = int(sys.argv[17])
    try:
        remove_jackpotting = sys.argv[18]
    except IndexError:
        remove_jackpotting = ""
    update_bic(input_filename, ncluster, max_its, nreads, muts_threshold, mask_threshold, mask_index,
               na_threshold, max_mut_prob, deletions, workDir, EM_its, is_d1, queue, trim_size, remove_jackpotting, dmuts=dmuts)


