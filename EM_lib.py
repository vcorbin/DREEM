# Contain useful function that could be used in different stage of the EM
import numpy
import pandas
try:
    # non-website import
    from denominator import smart_denominator
except:
    from .denominator import smart_denominator


def bitString_to_bitVect(bitString, deletions='1', unknowns = 0):
    """
    Transform a bitString into a bitVect numpy array
    :param bitString: a string containing '0', '1', 'A', 'C', 'T', 'G', 'D', '?' and '.'
    :param deletions: Needs to be '0', '1' or '2'. If '0' it will be replaced by 0.0, if '1' it will be replaced by 1.0, 
    if '2' it will be replaced by NA. 
    :return: a numpy array where the '0', and ('?', '.') of the bitString have been replaced by 0 and numpy.nan 
             respectively, while the rest are converted to 1.
    """
    assert deletions in ['0', '1', '2'], "Error: deletion must be '0', '1' or '2'"
    bitVect = numpy.array(list(bitString))
    bitVect[bitVect == '?'] = '2'
    bitVect[bitVect == '.'] = '2'
    bitVect[bitVect == 'N'] = '2'
    bitVect[bitVect == 'A'] = '1'
    bitVect[bitVect == 'C'] = '1'
    bitVect[bitVect == 'G'] = '1'
    bitVect[bitVect == 'T'] = '1'
    bitVect[bitVect == 'D'] = deletions
    num_mutations = "".join(bitVect).count('1')
    bitVect = bitVect.astype(float)  # Want it to be floats in order to avoid rounding of division in later calculations
    bitVect[bitVect == 2.0] = unknowns
    return bitVect, num_mutations


def bitVect_to_bitString(bitVect):
    """
    Transform a bitVect numpy.array into a bitString
    :param bitVect: numpy.array containing 0, 1 and numpy.nan
    :return: a string version of the bitVect containing '?' instead of numpy.nan.
    """
    bitVectTemp = bitVect.copy()
    bitVectTemp[bitVectTemp == numpy.nan] = 2.0
    bitVectTemp = bitVectTemp.astype(int).astype(str)
    bitVectTemp[bitVectTemp == "2"] = "?"
    bitString = "".join(bitVectTemp)
    return bitString


def expand_df(df, front, behind):
    """
    Expands a dataFrame by adding 0s
    :param df: pandas dataFrame
    :param front: int. Number of 0s in front of the df
    :param behind: int. Number of 0s behind of the df
    :return: new dataFrame with the correct dimension
    """
    new_df = pandas.DataFrame(numpy.zeros((front + df.shape[0] + behind, front + df.shape[1] + behind)))
    new_df.loc[front:front+df.shape[0]-1, front:front+df.shape[1]-1] = df.values
    return new_df


def normalise_clusters(clust_df):
    """
    Normalises the cluster signals using the median of the 10 highest signals. signals higher than the median become 1.
    Transforms a cluster signals into a 1D constraint files (1 for each clusters).
    :param clust_df: pandas.DataFrame. Contain the signal for the clusters.
    :return: pandas.DataFrame of same dimensions as clust_df
    """
    Nmax = 10  # Number of peaks used to calculate the normalization factor
    norm_df = clust_df.copy()
    for i in range(len(norm_df)):
        #print column
        #norm = numpy.median(numpy.sort(norm_df[column].values)[-Nmax:])
        norm = numpy.sort(norm_df[i])[-5]
        norm_df[i] = norm_df[i] / norm
    norm_df[norm_df > 1.0] = 1.0
    return norm_df


#TODO: Have to take the basis into account.
def remove_jackpotting(x, mu_max=None):
    """
    Calculates the maximum probability of each read to be in the data, and uses it to remove jackpotted reads. 
    :param x: object of class ind_bit_vect
    :return: modify x. The parameter x.N is modified
    """
    if mu_max is None:
        mu_max = numpy.nansum(x.N * x.bitMatrix.values.T, axis=1) / \
                numpy.nansum(x.N * numpy.isfinite(x.bitMatrix.values).T, axis=1)
    N_total = numpy.sum(x.N)
    new_N = (x.bitMatrix * mu_max + (1.0 - x.bitMatrix) * (1.0 - mu_max)).product(axis=1) * N_total
    N_max = numpy.maximum(1.0, new_N)
    x.N = numpy.array(numpy.minimum(x.N, N_max))


def calculate_observed_cluster_probs(responsibilities, x):
    """
    calculate the  observed cluster_probs using the responsibilities
    :param responsibilities: numpy array of shape (N,K). can be extracted form a bitvect_resp object with 
                             bitvect_resp.responsibilities
    :param x: indbitvect object
    :return: a list of length K containing the abundance of each cluster in the observed data
    """
    cluster_probs = numpy.sum(responsibilities.T * x.N, axis=1) / x.N.sum()
    return cluster_probs


def obs_cluster_probs_to_real(obs_cluster_probs, mu, dmuts=3, normalise=True):
    """
    calculate the real cluster_prob from the observed cluster_prob
    :param obs_cluster_probs: list of length K. observed cluster_prob
    :param mu: numpy.array of shape (K, I). mutation probability for each cluster 
    :return: non-normalised (unless specified) real cluster_prob, used to calculate log-likelihood
    """
    K = len(obs_cluster_probs)
    assert mu.shape[0] == K, \
        "Error: mu doesn't have the right dimension. It's length should be the same as the length of obs_cluster_probs."
    real_cluster_probs = numpy.array([obs_cluster_probs[k] / smart_denominator(mu[k], dmuts=dmuts)[0] for k in range(K)])
    if normalise:
        real_cluster_probs = real_cluster_probs / sum(real_cluster_probs)
    return real_cluster_probs


def real_cluster_probs_to_obs(real_cluster_probs, mu, dmuts=3, normalise=True):
    """
    calculate the observed cluster_prob from the real cluster_prob
    :param real_cluster_probs: 
    :param mu: 
    :param normalise: 
    :return: 
    """
    K = len(real_cluster_probs)
    assert mu.shape[0] == K, \
        "Error: mu doesn't have the right dimension. It's length should be the same as the length of obs_cluster_probs."
    obs_cluster_probs = numpy.array([real_cluster_probs[k] * smart_denominator(mu[k], dmuts=dmuts)[0] for k in range(K)])
    if normalise:
        obs_cluster_probs = obs_cluster_probs / sum(obs_cluster_probs)
    return obs_cluster_probs


def calculate_cluster_probs(responsibilities, x, mu, dmuts=3, normalise=True):
    """
    calculate the real cluster_probs (including unobserved data) using the responsibilities
    :param responsibilities: numpy array of shape (N,K). can be extracted form a bitvect_resp object with 
                             bitvect_resp.responsibilities
    :param x: indbitvect object
    :param mu: numpy array of shape (K, I) containing the mutation probabilities
    :return: a list of length K containing the real abundance of each cluster
    """
    obs_cluster_probs = calculate_observed_cluster_probs(responsibilities, x)
    cluster_probs = obs_cluster_probs_to_real(obs_cluster_probs, mu, dmuts=dmuts, normalise=normalise)
    return cluster_probs


# TODO: fix. Doesn't deal with question marks in the right way. Actually crashes if remove_reads=False
def replace_nans(x, resp, mu, cluster_probs, remove_reads=True):
    # remove_reads will remove all the reads with question marks in them instead of probabilistically assign the
    # question marks to 0 or 1.
    I = mu.shape[1]
    K = mu.shape[0]
    new_bitM = []
    new_N = []

    assert mu.shape[0] == len(cluster_probs), "mus and cluster_probs don't have compatible sizes: %d vs %d" \
                                              % (mu.shape[0], len(cluster_probs))
    assert x.bitMatrix.shape[1] == mu.shape[1], "The sizes of bitMatrix and mu are not compatible: %d vs %d" \
                                                % (x.bitMatrix.shape[1], mu.shape[1])
    mut_probs = numpy.zeros((K, I))
    for iclust in range(K):
        mut_probs[iclust] = mu[iclust] * cluster_probs[iclust]

    for i, bitvect in enumerate(x.bitMatrix.values):
        cond = numpy.isnan(bitvect)
        if len(bitvect[cond]) == 0:
            new_bitM.append(bitvect)
            new_N.append(x.N[i])
        else:
            if not remove_reads:
                nan_mut_probs = 0
                for k in range(K):
                    nan_mut_probs += resp.responsibilities[i][k] * mut_probs[k]
                for j in range(x.N[i]):
                    nan_mut_probs = nan_mut_probs[cond]
                    bitvect[cond] = numpy.ones(len(nan_mut_probs)) * \
                                    (numpy.random.rand(len(nan_mut_probs)) <= nan_mut_probs)
                    new_bitM.append(bitvect)
                    new_N.append(1)
            else:
                pass

    x.bitMatrix = pandas.DataFrame(new_bitM)
    x.N = numpy.array(new_N)




















