#!/usr/bin/python

import EM_plots as plots
import os
import pandas
import EM_classes as classes
import EM_files as files
import numpy
import scipy
import denominator as denom
from EM_lib import replace_nans


def jackPotting(resp, mu, deletions, cluster_probs=None, remove_reads=True, workDir=None):
    """

    :param resp: either an object of class bitvect_resp or a responsibilities file name 
    :param mu: either a 2d numpy.array of mutation probabilities for each cluster or a mu file name
    :param cluster_probs: array containing the relative abundance of each cluster or a cluster_prob file name
    :param remove_reads: if true all the reads with question marks will be removed instead of the question marks being
    probabilistically assigned to 0 or 1
    :param workDir: name of the directory in which the files and plot generated will be saved
    :return: 
    """

    if workDir is None:
        if type(mu) == str:
            workDir = os.path.dirname(mu)
            if workDir == "":
                workDir = "."
        else:
            workDir = "."
    filename = workDir + '/jackpotting'

    if type(resp) == str:
        resp = files.loadResp(resp)
    x = resp.to_indbitvect(deletions)
    if type(mu) is str:
        print("Loading mus")
        mu = files.LoadClusters(mu)
    if cluster_probs is not None:
        if type(cluster_probs) is str:
            print("Loading cluster_probs")
            cluster_probs = files.load_cluster_prob(cluster_probs)
    else:
        cluster_probs = numpy.sum(resp.responsibilities.T * resp.N, axis=1) / resp.N.sum()

    jackpot = jackPotting_calculation(x, resp, mu, cluster_probs, remove_reads=remove_reads)

    plots.PlotJackPotting(jackpot, filename)
    jackpot.to_csv(filename + '.txt')


def jackPotting_calculation(x, resp, mu, cluster_probs, remove_reads=True):
    """

    :param x: object of class ind_bit_vect
    :param mu: 2d numpy array of mutation probabilities for each cluster (shape is K*I)
    :param cluster_probs: array containing the relative abundance of each cluster
    :param remove_reads: if true all the reads with question marks will be removed instead of the question marks being
    probabilistically assigned to 0 or 1
    :param plot_name: string, file name without extension that will be used as prefix for the various plots
    :return: an object of class jackPotting
    """

    assert mu.shape[0] == len(cluster_probs), "mus and cluster_probs don't have compatible sizes: %d vs %d" \
                                              % (mu.shape[0], len(cluster_probs))
    assert x.bitMatrix.shape[1] == mu.shape[1], "The sizes of bitMatrix and mu are not compatible: %d vs %d" \
                                                % (x.bitMatrix.shape[1], mu.shape[1])

    I = x.bitMatrix.shape[1]
    K = len(cluster_probs)

    # Replacing the nan with appropriate 0 or 1.
    new_x = x.copy()
    print("Replacing nans")
    replace_nans(new_x, resp, mu, cluster_probs, remove_reads=remove_reads)
    print("Done replacing the nans")
    nreads = x.N.sum()

    bitdict = {}
    for row, bitvect in enumerate(new_x.bitMatrix.values):
        key = ''.join(map(str, map(int, new_x.bitMatrix.values[row])))
        try:
            bitdict[key][0] += new_x.N[row]
        except KeyError:
            bitdict[key] = [new_x.N[row], bitvect, key]
    N_obs, bitvects, bitStrings = zip(*bitdict.values())
    N_obs = numpy.array(N_obs)
    bitvects = numpy.array(bitvects)
    bitStrings = numpy.array(bitStrings)

    N = len(N_obs)
    Px = numpy.zeros((N, I, K))  # Log of probability of a base in vector Xn given mus in cluster k
    L_corr_clust = numpy.zeros(K)
    for k in range(K):
        Px[:, :, k] = scipy.stats.bernoulli.logpmf(bitvects, mu[k])
        L_corr_clust[k] = numpy.log(denom.smart_denominator(mu[k])[0])
    Px = numpy.nansum(Px, axis=1)
    log_prob = L_corr_clust - numpy.log(cluster_probs)
    Px = numpy.subtract(Px, log_prob)
    Px = numpy.sum(numpy.exp(Px), axis=1)

    print("Px is calculated")

    N_theory = Px * nreads
    nmuts = numpy.nansum(bitvects, axis=1)
    jackpot = classes.jackPotting(N_theory, N_obs, nmuts, bitStrings)

    return jackpot


if __name__ == "__main__":
    import argparse

    """
    List of arguments and corresponding variables

    Mandatory:
                     resp --> resp
                       mu --> mu

    Optional:
       -R, --remove_reads --> remove_reads
      -C, --cluster_probs --> cluster_probs
            -W, --workDir --> workDir
    """

    parser = argparse.ArgumentParser(
        description='Make a plot describing the jackpotting.')

    parser.add_argument('resp', type=str, help='responsibilities filename.')
    parser.add_argument('mu', type=str, help='Cluster parameters filename.')
    parser.add_argument('-R', '--remove_reads', dest='remove_reads', action='store_true',
                        help="If specified, the reads with '?' are removed. Otherwise the '?' are" +
                             " probabilistically assigned to 0 or 1.")
    parser.add_argument('-C', '--cluster_probs', type=str, dest='cluster_probs', default=None,
                        help="cluster probs filename.")
    parser.add_argument('-W', '--workDir', type=str, dest='workDir', default=None,
                        help='Directory to record the plots.')

    args = parser.parse_args()
    resp = args.resp
    mu = args.mu
    cluster_probs = args.cluster_probs
    remove_reads = args.remove_reads
    workDir = args.workDir

    jackPotting(args.resp, args.mu, cluster_probs=args.cluster_probs, remove_reads=args.remove_reads,
                workDir=args.workDir)
