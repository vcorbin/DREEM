import numpy as np
import scipy
from denominator import smart_denominator


def ll(BVM, Ns, mu, pi_obs, dmuts=3):
    """
    Calculate the loglikelihood
    :param BVM: bitvector matrix as a numpy 2d array
    :param Ns: number of time each read is observed (numpy 1d array)
    :param mu: mu parameters (numpy 2d array)
    :param pi_obs: list of the observable pi. same length as the mu
    :return: loglikelihood
    """
    # TODO: check mu[0] is the same length as BVM.shape[1], length(mu) = length(pi_obs) = K
    K = len(mu)
    N = len(BVM)
    log_Px = np.zeros((N, K))
    log_denom = np.zeros(K)

    for k in range(K):
        log_denom[k] = np.log(smart_denominator(mu[k], dmuts=dmuts)[0])
        log_Pxik = scipy.stats.bernoulli.logpmf(BVM, mu[k])  # contain the logs of probs for each base of each read for each cluster of shape (N, I)
        log_Px[:, k] = log_Pxik.sum(axis=1) - log_denom[k]
    reads_log_Px = scipy.misc.logsumexp(np.add(np.log(pi_obs), log_Px), axis=1)
    loglik = np.dot(reads_log_Px, Ns).sum()
    return loglik
