#!/usr/bin/python3

import EM_binary as EM
import EM_files as files
import time
import multiprocessing
import numpy
from jackPotting import jackPotting
from glob import glob
from PrettyPlots import run_PrettyPlots
from EM_lib import remove_jackpotting


# Important: For now max_threshold is an argument that does not seem to be used. Needs to be implemented in load_binary.

def runMultipleItersEM(x, refs, d=2, threading=False, EM_its=10, init_clusters=None, max_its=200,
                       seeds=None, workDir="./", sub_dir_prefix="run_", dmuts=3):
    if seeds is None:
        print("Warning: no seed provided, the seeds for the individual runs will be generated randomly.")
        seeds = numpy.random.randint(4294967295,size=EM_its)

    nprocess = min(EM_its, multiprocessing.cpu_count())
    # Execute the 1-cluster run
    print("Run 1 cluster")
    d1_lls, d1_clusters, d1_cps = EM.runEM(x, 1, refs, threading=threading, max_its=min(20, max_its),
                                           workDir=workDir + '/run_d1', init_clusters=init_clusters, dmuts=dmuts)

    print("Start multiprocessing with %d parallel processes" % nprocess)
    mprocesses = []
    pool = multiprocessing.Pool(processes=nprocess)
    for its in range(EM_its):
        print(its + 1, seeds[its])
        workDir_its = workDir + '/' + sub_dir_prefix + str(its + 1)
        print(threading)
        mprocesses.append(pool.apply_async(EM.runEM, args=(x, d, refs),
                                           kwds=dict(threading=threading, max_its=max_its,
                                                     workDir=workDir_its,
                                                     seed=seeds[its], init_clusters=init_clusters)))
    pool.close()
    pool.join()

    all_logliks = []
    for its, process in enumerate(mprocesses):
        print("getting the log-likelihood for run %d" % (its + 1))
        [lls, mus, cps] = process.get()
        all_logliks.append(lls)
    print("End of multiprocessing")

    # Save the likelihoods
    #all_logliks = [d1_lls] + all_logliks
    #files.write_lls(all_logliks, workDir)


def post_process(workDir, sub_dir_prefix, EM_its, deletions, is_d1="no", trim_size=0):
    if is_d1 == "yes":
        best_run = "d1"
    else:
        best_run = files.clean_up(workDir, sub_dir_prefix, EM_its)  # Combine and save the lls and bics
    best_run_dir = workDir + '/' + sub_dir_prefix + best_run
    print(best_run_dir)
    cluster_files = glob(best_run_dir + '/' + "clusters*.txt")  # avoid taking constraint files as cluster
    cluster_files.sort(key=len, reverse=False)
    mu = cluster_files[0]  # the path of the cluster file
    responsibilities = best_run_dir + '/' + "responsibilities.txt"
    files.mutation_rate_histogram(workDir + '/' + sub_dir_prefix + "d1" + "/responsibilities.txt", deletions)
    jackPotting(responsibilities, mu, deletions)
    run_PrettyPlots(mu, ymax=0.3)

    #os.system("bsub -J tf run_trimFold.py " + mu + " " + str(trim_size[0]) + " " + str(trim_size[1]))


if __name__ == "__main__":

    t0 = time.time()
    print(t0)
    args, workDir, prefix, seeds = files.set_up()
    x, refs = files.load_binary(args.input_filename, muts_threshold=args.muts_threshold,
                                max_mut_prob=args.max_mut_prob, nreads=args.nreads, mask_index=args.mask_index,
                                na_threshold=args.na_threshold, deletions=args.deletions, dmuts=args.dmuts)
    assert x is not None, "Binary vectors matrix could not be loaded."
    assert refs is not None, "List of references could not be loaded."

    x = x.collapse_similar_reads()
    if args.remove_jackpotting:
        remove_jackpotting(x)

    print("Number of reads:", refs['refs'].shape[0])
    print("Number of unique reads:", x.bitMatrix.shape[0])

    if args.threading:
        print("Threading on")

    sub_dir_prefix = "run_"

    if args.EM_its > 1:  # This is never used
        runMultipleItersEM(x, refs, d=args.ncluster, threading=args.threading, EM_its=args.EM_its,
                              max_its=args.max_its, workDir=workDir, sub_dir_prefix = sub_dir_prefix, seeds=seeds, dmuts=args.dmuts)
        post_process(workDir, sub_dir_prefix, args.EM_its, args.deletions) #args.filename

    else:   # args.EM_its = 1 always (for Clustering)
        lls, clusters, cps = EM.runEM(x, args.ncluster, refs, workDir=workDir, threading=args.threading,
                                      max_its=args.max_its, seed=args.seed, dmuts=args.dmuts)

    t1 = time.time()
    print("Run time:", t1-t0)

