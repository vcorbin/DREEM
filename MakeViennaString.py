#!/usr/bin/python
try:
    from . import EM_files as files
except:
    import EM_files as files
import numpy


def makeViennaStringFromFile(cluster_filename, cluster_index, open_threshold, closed_bases=None):
    mu = files.LoadClusters(cluster_filename)[cluster_index]
    return makeViennaString(mu, open_threshold, closed_bases=closed_bases)

def makeViennaString(mu,open_threshold, closed_bases=None):
    vienna_string = []
    for base in mu:
        if base > open_threshold:
            vienna_string.append('x')
        else:
            vienna_string.append('.')
    vienna_string = numpy.array(vienna_string)
    if closed_bases is not None:
        vienna_string[closed_bases] = '|'
    return ''.join(vienna_string)

if __name__ == "__main__":
    import argparse

    '''
    List of arguments and corresponding variables:
    
    Mandatory:
                cluster_filename --> cluster_filename
                cluster_index --> cluster_index
                threshold --> threshold
    
    Optional:
                -c, --closed --> closed_bases
    '''

    parser = argparse.ArgumentParser(
        description='Produce a Vienna format string from a list of mutation probabilities.')

    parser.add_argument('cluster_filename', type=str,
                        help='Name of the cluster file containing the mutation Probs for each cluster.')

    parser.add_argument('cluster_index', type=int,
                        help='Which cluster that will be used to generate the vienna string.')

    parser.add_argument('open_threshold', type=float,
                        help='Float between 0 and 1. Threshold for the mut probs above which a base is consideredopen.')

    parser.add_argument('-c', '--closed', type=str, dest='closed_bases', default=None,
                        help='List of position that are forced to be closed, each separated by a comma.')

    args = parser.parse_args()

    # Variable conversion
    if args.closed_bases is not None:
        closed_bases = map(int,args.closed_bases.split(','))
    else:
        closed_bases = args.closed_bases

    print (makeViennaStringFromFile(args.cluster_filename, args.cluster_index, args.open_threshold, closed_bases))
