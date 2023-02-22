#!/usr/bin/python

import os
import pandas
import EM_files as files
import matplotlib
import pylab
import numpy
import scipy
from collections import defaultdict
from matplotlib import pyplot as plt
import denominator as denom
import scipy
matplotlib.use('Agg')  # Remove the need to run a X server

def jackPotting(bitVector, mu, responsibilities, remove_reads=False):
    workDir = os.path.dirname(mu)
    filename = workDir + '/jackpotting'
    print filename
    x = PlotJackPotting(bitVector, mu, responsibilities, plot_name=filename, remove_reads=remove_reads)
    x = pandas.DataFrame(x)
    x.to_csv(filename+'.txt', sep="\t", header=['N_theory', 'N_obs', 'N_sd', 'nmut', 'bitString'], index_label=False)

def PlotJackPotting(x, mu, responsibilities, plot_name=None, remove_reads=False):
    # remove_reads will remove all the reads with question marks in them instead of probabilistically assign the
    # question marks to 0 or 1.

    if type(x) is str:
        print "Loading bitvectors"
        x, labels = files.load_binary(x)
    if type(mu) is str:
        print "Loading mus"
        mu = files.LoadClusters(mu)
    if type(responsibilities) is str:
        print "Loading responsibilities"
        responsibilities = files.loadResp(responsibilities).responsibilities

    cluster_probs = numpy.sum(responsibilities.T * x.N, axis=1) / x.N.sum()
    I = x.bitMatrix.shape[1]
    K = len(cluster_probs)

    # Replacing the ? with appropriate 0 or 1.
    new_x = x.copy()
    print "Replacing the nans"
    replace_nans(new_x, mu, cluster_probs,remove_reads=remove_reads)
    print "Done replacing the nans"
    nreads = x.N.sum()

    bitdict = {}
    for row, bitvect in enumerate(new_x.bitMatrix.values):
        key = ''.join(map(str,map(int,new_x.bitMatrix.values[row])))
        try:
            bitdict[key][0] += new_x.N[row]
        except KeyError:
            bitdict[key] = [new_x.N[row], bitvect, key]
    Ns, bitvects, bitStrings = zip(*bitdict.values())
    Ns = numpy.array(Ns)
    bitvects = numpy.array(bitvects)
    bitStrings = numpy.array(bitStrings)

    N = len(Ns)
    Px = numpy.zeros((N, I, K))  # Log of probability of a base in vector Xn given mus in cluster k
    L_corr_clust = numpy.zeros(K)
    for k in xrange(K):
        Px[:, :, k] = scipy.stats.bernoulli.logpmf(bitvects, mu[k])
        L_corr_clust[k] = numpy.log(denom.smart_denominator(mu[k])[0])
    Px = numpy.nansum(Px, axis=1)
    log_prob = L_corr_clust - numpy.log(cluster_probs)
    Px = numpy.subtract(Px,log_prob)
    Px = numpy.sum(numpy.exp(Px),axis=1)

    print "Px is calculated"

    N_theory = numpy.ceil(Px * nreads)
    sd = numpy.sqrt(N_theory)
    print "sd is calculated"
    N_dist = (Ns - N_theory) / sd
    print "N_dist is calculated"

    nmuts = numpy.nansum(bitvects,axis=1)
    ordind = N_theory.argsort()
    N_theory = N_theory[ordind]
    Ns = Ns[ordind]
    nmuts = nmuts[ordind]
    N_dist = N_dist[ordind]
    bitStrings = bitStrings[ordind]
    print "N_dist is ordered"

    # Plots
    pylab.figure()
    pylab.scatter(range(len(N_dist)), N_dist)
    pylab.title('Deviation from the expected number of occurences')
    pylab.ylabel('Deviation (standard deviations)')
    pylab.xlabel('Reads (ordered by expected number of occurences)')
    if plot_name is not None:
        pylab.savefig(plot_name + '_sds.png')
        pylab.close()
    pylab.figure()
    pylab.scatter(N_theory, N_dist)
    pylab.title('Deviation from the expected number of occurences vs expected number of occurences')
    pylab.ylabel('Deviation (standard deviations)')
    pylab.xlabel('N_expected')
    pylab.xscale('log')
    if plot_name is not None:
        pylab.savefig(plot_name + '_sds_vs_N.png')
        pylab.close()
    return zip(N_theory,Ns,N_dist,nmuts,bitStrings)

def replace_nans(x, mu, cluster_probs,remove_reads=False):
    I = mu.shape[1]
    K = mu.shape[0]
    new_bitM = []
    new_N = []

    assert mu.shape[0] == len(cluster_probs), "mus and cluster_probs don't have compatible sizes: %d vs %d" % (mu.shape[0], len(cluster_probs))
    assert x.bitMatrix.shape[1] == mu.shape[1], "The sizes of bitMatrix and mu are not compatible: %d vs %d" % (x.bitMatrix.shape[1], mu.shape[1])

    mut_probs = numpy.zeros(I)
    for iclust in xrange(K):
        mut_probs = mut_probs + mu[iclust] * cluster_probs[iclust]

    for i, bitvect in enumerate(x.bitMatrix.values):
        cond = numpy.isnan(bitvect)
        if len(bitvect[cond]) == 0:
            new_bitM.append(bitvect)
            new_N.append(x.N[i])
        else:
            if not remove_reads:
                for j in xrange(x.N[i]):
                    nan_mut_probs = mut_probs[cond]
                    bitvect[cond] = numpy.ones(len(nan_mut_probs)) * (numpy.random.rand(len(nan_mut_probs)) <= nan_mut_probs)
                    new_bitM.append(bitvect)
                    new_N.append(1)
            else:
                pass

    x.bitMatrix = pandas.DataFrame(new_bitM)
    x.N = numpy.array(new_N)

if __name__ == "__main__":

    import argparse

    '''
    List of arguments and corresponding variables

    Mandatory:
                bitVector --> bitVector
                       mu --> mu
         responsibilities --> responsibilities

    Optional:
                  -R, --remove_reads --> remove_reads
    '''

    parser = argparse.ArgumentParser(
        description='Make a plot describing the jackpotting.')

    parser.add_argument('bitVector', type=str, help='Bitvector filename.')
    parser.add_argument('mu', type=str, help='Cluster parameters filename.')
    parser.add_argument('responsibilites', type=str, help='Responsibilites filename.')
    parser.add_argument('-R', '--remove_reads', dest='remove_reads', action='store_true',
                        help="If specified, the reads with '?' are removed. Otherwise the '?' are" +
                        " probabilistically assigned to 0 or 1.")

    args = parser.parse_args()
    bitVector = args.bitVector
    mu = args.mu
    responsibilities = args.responsibilites
    remove_reads = args.remove_reads

    jackPotting(bitVector, mu, responsibilities, remove_reads=remove_reads)
