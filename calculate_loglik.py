# Contain useful function that could be used in different stage of the EM
import numpy
import scipy.misc, scipy.stats
from EM_classes import ind_bit_vect
from EM_classes import bitvect_resp
from EM_lib import calculate_cluster_probs
from EM_lib import obs_cluster_probs_to_real
from EM_files import loadResp
from EM_files import LoadClusters
from EM_binary import Expectation_step
from denominator import smart_denominator


def calculate_ll(bitvects, mu, obs_cluster_probs=None, return_cluster_probs = False, dmuts=3):
    """
    Calculate the log-likelihood given a set of reads, mutation probabilities and cluster_probs. 
    Not used in the EM where the ll is calculated while updating the responsibilities. 
    Used to compare different solutions.
    :param bitvects: ind_bit_vect, bitvect_resp or responsibilities file name
    :param mu: list of list containing the mutation probs for the diverse clusters or file containing that data
    :param obscluster_probs: list of the cluster proportions (of the observed data, not the real data). Needs to be 
                          specified only if bitvects is a ind_bit_vect, otherwise it can be calculated from the 
                          responsibilities information
    :param return_cluster_probs: default False. If True the real_cluster_probs will be returned in addition to the 
                                 loglik. Used more to debug than anything
    :return: float, log-likelihood
    """
    # Check the format of the input
    data_format = type(bitvects)
    assert data_format is ind_bit_vect or data_format is bitvect_resp or data_format is str, \
        "Error: bitvects must be a ind_bit_vect, a bitvect_resp or a responsibilities file name."
    assert type(mu) is numpy.ndarray or type(mu) is str, \
        "Error: mu must be a numpy.ndarray or a cluster file (mutation rates)."

    # set up the data in the right format

    if type(mu) is str:
        mu = LoadClusters(mu)
    else:
        pass

    if obs_cluster_probs is not None:
        cluster_probs = obs_cluster_probs_to_real(obs_cluster_probs, mu, dmuts=dmuts, normalise=False)

    if data_format is ind_bit_vect:
        assert obs_cluster_probs is not None, \
            "Error: if bitvects is given as a ind_bit_vect object then cluster_probs must be specified."
        x = bitvects
    else:
        if data_format is str:
            resp = loadResp(bitvects)  # load bitvects into a bitvect_resp object
        else:  # bitvects is therefore a bitvect_resp object
            resp = bitvects
        x = resp.to_indbitvect()
        if obs_cluster_probs is None:
            cluster_probs = calculate_cluster_probs(resp.responsibilities, x, mu, dmuts=dmuts, normalise=False)
        else:
            print("cluster_prob has been specified, so it will be used and not calculated from the responsibilities.")

    _, loglik = Expectation_step(x, mu, numpy.log(cluster_probs))

    if return_cluster_probs:
        return loglik, cluster_probs / sum(cluster_probs)
    else:
        return loglik



















