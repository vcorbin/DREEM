#!/usr/bin/python
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
import scipy.stats
import os
try:
    from . import EM_lib as lib
    from . import EM_files as files
except:
    import EM_lib as lib
    import EM_files as files

np.set_printoptions(threshold=np.inf)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def run_bpc(f, args):
    expand = np.array(args.expand.split(","), dtype=int)

    outfile_prefix= os.path.splitext(f)[0]

    correlation_matrix= get_cov(f)
    text = np.empty((len(correlation_matrix), len(correlation_matrix)), dtype=str)
    ref = ['A']*len(correlation_matrix)
    for i in range(len(correlation_matrix)):
        for j in range(len(correlation_matrix)):
            text[i][j] = "%s , %s" % (ref[i], ref[j])

    print("Starting Z scores")

    z_scores = get_Zscores(correlation_matrix)
    mod_Z = z_scores.copy()
    mod_Z[mod_Z<0] = 0.0
    mod_Z = -mod_Z

    print("Finished mod Z")

    # Expand the matrices
    correlation_df = lib.expand_df(pd.DataFrame(correlation_matrix), expand[0], expand[1])
    mod_Z_df = lib.expand_df(pd.DataFrame(mod_Z), expand[0], expand[1])
    z_scores_df = lib.expand_df(pd.DataFrame(z_scores), expand[0], expand[1])

    correlation_df.to_csv(outfile_prefix + "_covariance.txt", sep=" ", header=False, index=False)
    mod_Z_df.to_csv(outfile_prefix + "_Zscores.txt", sep=" ", header=False, index=False)

    print("Starting plots")

    # Change the scale of the correlation matrix. There should be a better way to do that without changing the
    # values of the matrix, but using plot options instead.
    #ymax = max(correlation_matrix.max(),abs(correlation_matrix.min()))
    #correlation_matrix[correlation_matrix>0.0] = correlation_matrix[correlation_matrix>0.0] * ymax /correlation_matrix.max()
    #correlation_matrix[correlation_matrix<0.0] = correlation_matrix[correlation_matrix<0.0] * ymax / abs(correlation_matrix.min())

    data = [go.Heatmap(z=correlation_df.values, text=text, hoverinfo="x+y+text")]
    fig = go.Figure(data=data)
    plotly.offline.plot(fig, filename=outfile_prefix + '_covariance.html', auto_open=args.open)

    data = [go.Heatmap(z=z_scores_df.values, text=text, hoverinfo="x+y+text")]
    fig = go.Figure(data=data)
    plotly.offline.plot(fig, filename=outfile_prefix + '_Z_scores.html', auto_open=args.open)


def get_cov(filename):
    x = files.load_binary(filename)[0]
    correlation_matrix = cov(x)
    np.fill_diagonal(correlation_matrix, 0)
    return correlation_matrix


def cov(x):  # Warning: modifies the object x
    N_bases = x.bitMatrix.shape[1]
    x.bitMatrix = x.bitMatrix.replace(np.nan,0)
    k_tests = x.get_Ntotal()
    averages_of_bases = x.get_pop_average()
    reads = x.bitMatrix.values
    cov_transposed = (reads - averages_of_bases).T * np.sqrt(x.N)

    covar = cov_transposed.T

    cov_matrix = np.dot(cov_transposed, covar)

    cov_matrix = np.divide(cov_matrix, k_tests)

    std = x.get_std()

    std = std.reshape(1, N_bases)
    std_t = std.reshape((N_bases, 1))
    std_matrix = np.dot(std_t, std)

    with np.errstate(divide='ignore', invalid='ignore'):
        cov_matrix = np.true_divide(cov_matrix, std_matrix)

    return cov_matrix


def get_Zscores(M):
    """
    Calculate the Z-scores for a dataFrame.
    :param M: numpy.ndarray (NxN Matrix)
    :return: numpy.ndarray (NxN Matrix)
    """
    zscores = scipy.stats.mstats.zscore(M, axis=0, ddof=0)
    return zscores

def boolean(str):
    if str == "True":
        return True
    return False


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Run a EM to cluster reads from a bam/sam file or a binary_vector '
                                                 'file')
    parser.add_argument('input_filename', type=str,
                        help='Preferably a binary_vector file (.txt) created by sam_to_binary(). If a sam/bam is given, '
                             'a binary_vector file of the same name will be created.')
    parser.add_argument('-s', '--sep', type=str, dest='sep', default='\t', help='Separation used in the bit_vector.')
    parser.add_argument('-o', '--open', type=boolean, dest='open', default=False, help='False if should not open heatmap')
    parser.add_argument('-e', '--expand', type=str, dest='expand', default='0,0',
                        help='2 integer separated by a comma. Number of 0s to put in front and behind the correlation '
                             'matrices.')

    args = parser.parse_args()

    run_bpc(args.input_filename, args)