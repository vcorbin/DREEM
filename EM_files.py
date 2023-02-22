import numpy
import pandas
import os
import sys
import time
import subprocess
import arguments
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
try:
    from . import EM_lib as lib
    from . import EM_classes as EMClass
except:
    import EM_lib as lib
    import EM_classes as EMClass

def web_read_bitvector(file, comment='@'):
    '''
    Read a bitvector file into memory for use in the website TODO: Just noticed the comment class. (see function load_header() below)
    Eventually, should be integrated into comment class, but I need to keep track of ERRORS to report
    in an email to the user.

    @param: filename - bitvector file name.
    @param: comment - character which indicates which lines to skip

    @return: regions - list of dictionaries, for each reference region in the bitvector,
                       containing reference, start, end, length, sequence, and number of reads (reads)
    @return: bitvector - dataframe containng the string bitvectors
    @return: ERRORS - A list of errors generated when trying to parse the bitvector header
    '''
    bitvector = pandas.read_csv(file, sep='\t', comment=comment)

    ERRORS = []
    regions = []
    refs = []
    with open(file) as f:
        for line in f:
            if line[0] != comment:
                # Only read comment lines to generate dictionary
                break
            # Clean up line for processing
            line = line.strip(comment)
            line = line.strip('\n')
            line = line.split('\t')
            # Populate dictionary matching key (line[0]) to value (line[2])
            ref = line[1]
            # If we are looking at a new reference region, add in to the list and create a new region
            if ref not in refs:
                if len(refs) != 0:
                    # If this is not the first region we are looking at, add the previous region to regions
                    regions.append(region)
                refs.append(ref)
                region = {'reference':ref}
            try:
            # Convert integer values to integers
                value = int(line[2])
            except:
                value = line[2]
            region[line[0]] = value
    # Add the last region to the list of regions

    # TODO: Modify this to accept multiple reference regions within the same bitvector file - QUESTION, is this every necessary??
    if 'start' not in region:
        ERRORS.append('Bitvector heading does not contain the start of the region. Please add to the header the following tab delimited line \n @start\t<ref>\it<start>, eg @start\tHIV\t455')
    if 'length' not in region:
        ERRORS.append('Bitvector heading does not contain the length of the region. Please add to the header the following tab delimited line containing the length of each bitvector \n @length\t<ref>\it<length>, eg @start\tHIV\t100')
    if 'seq' not in region:
        ERRORS.append('Bitvector heading does not contain the sequence of the region. Please add to the header the following tab delimited line containing the sequence of the bitvector. The sequence length must match the length of the bitvector  \n @length\t<ref>\it<seq>, eg @start\tHIV\tNCNNTCTGGTTAGACCAGANCTGAGCCTGGGAGCNCTCTGGCTAACTAGGGAACCCACNGCTNAAGCCNCAATAAAGCTTGCCTTGAGTGCTCAAANTAG')
    if 'seq' in region and 'length' in region:
        if region['length']!=len(region['seq']):
            ERRORS.append('The length of the sequence provided does not match the length of the bitvectors. Please ensure these values are accurate.')
    if not ERRORS:
        region['end'] = region['start']+region['length']
        region['nreads'] = len(bitvector['Reference']==region['reference'])
    regions.append(region)
    return regions, bitvector, ERRORS


def set_up(convert_args=True):
    args = arguments.parse_args(convert_args=convert_args)
    prefix = os.path.splitext(os.path.basename(args.input_filename))[0]

    if args.workDir is not None:
        workDir = args.workDir
    else:
        workDir = get_dir('.', args.ncluster, args.EM_its, args.max_its, prefix=prefix)

    if args.seed is None:
        args.seed = numpy.random.randint(4294967295)

    write_log(workDir, args)

    numpy.random.seed(seed=args.seed)
    seeds = numpy.random.randint(4294967295, size=args.EM_its)

    return args, workDir, prefix, seeds


def clean_up(dir_name, sub_dir_prefix, EM_its):
    list_sub_dir = [x[1] for x in os.walk(dir_name)][0]

    loglik_file = open(dir_name + '/loglikelihoods.txt', 'w')
    loglik_file.write('Runs\tLog_likelihood\tBIC_score\t \n')

    lls = []
    bics = []
    run_labels = []

    # Looking for run_d1:
    ll, bic = _get_ll_bic_run('d1', dir_name, sub_dir_prefix, list_sub_dir)
    lls.append(ll)
    bics.append(bic)
    run_labels.append('d1')

    # Combining the rest of the runs:
    for i in range(EM_its):
        run_index = i + 1
        ll, bic = _get_ll_bic_run(run_index, dir_name, sub_dir_prefix, list_sub_dir)
        lls.append(ll)
        bics.append(bic)
        run_labels.append(str(run_index))

    max_i = numpy.argmax(lls)
    best_run = [""] * len(lls)
    best_run[max_i] = "best run"
    for i in range(len(lls)):
        loglik_file.write(run_labels[i] + '\t' + str(lls[i]) + '\t' + str(bics[i]) + '\t' + best_run[i] + '\n')

    loglik_file.close()
    return run_labels[max_i]


def _get_ll_bic_run(run_index, workDir, sub_dir_prefix, list_sub_dir):
    sub_dir = sub_dir_prefix + str(run_index)
    if sub_dir in list_sub_dir:
        ll, bic = load_run_likelihood(workDir + '/' + sub_dir + "/loglikelihood.txt")
    else:
        print("Warning: run %s couldn't be found." % run_index)
        ll, bic = '-', '-'
    return ll, bic


def load_run_likelihood(ll_file):
    if type(ll_file) is str:
        ll_file = open(ll_file, 'r')
    ll = float(ll_file.readline())
    bic = float(ll_file.readline())
    ll_file.close()
    return ll, bic


def SaveRun(x, lls, run_resp, run_clusters, ref_labels, outDir):
    '''
    Saves one run of cluster data, updates the loglikelihood text file, outputs the cluster mutation probabilities,
    outputs the cluster responsibilities.
    @param: x  -
    @param: run_logL - float - the loglikelihood of this run_logL
    @param: run_resp -

    TO DO: FILL OUT DOCUMENTATION
    '''
    outDir = outDir + '/'
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    N_reads = ref_labels.shape[0]
    N_unique_reads = x.bitMatrix.shape[0]
    nreads_file = open(outDir + '/nreads.txt', 'w')
    nreads_file.write('Number of reads: ' + str(N_reads) + '\n')
    nreads_file.write('Number of unique reads: ' + str(N_unique_reads))
    nreads_file.close()

    run_logL = lls[-1]
    loglik_file = open(outDir + 'loglikelihood.txt', 'w')
    loglik_file.write(str(run_logL) + '\n')
    loglik_file.write(str(bic(run_logL, len(run_clusters), len(run_clusters[0]), x.N.sum())))
    loglik_file.close()

    all_lls_file = open(outDir + 'all_lls.txt', 'w')
    all_lls_file.write("Log_likelihood\n")
    for ll in lls:
        all_lls_file.write(str(ll) + "\n")
    all_lls_file.close()

    # Can be generalised later to allow plotting of multiple reference labels and co-ordinate systems
    gene, start = x.refs_info.get_coordinates()
    end = start + run_clusters.shape[1]

    SaveClusters(run_clusters, outDir + 'clusters_' + str(gene) + '_' + str(start) + '-' + str(end) + '.txt',
                 x.refs_info)

    SaveResponsibilities(x, run_resp, run_clusters, outDir)


def bic(ll, K, I, N):
    return numpy.log(N) * K * I - 2 * ll


def SaveResponsibilities(x, run_resp, run_clusters, outDir):
    # bitMatrix = x.bitMatrix.copy(deep=False)
    # bitMatrix.replace([0,1,numpy.nan],['0','1','?'],inplace=True)
    # bitMatrix = map("".join,numpy.asarray(bitMatrix))
    df = pandas.DataFrame(numpy.asarray(run_resp).astype(float),
                     columns=numpy.core.defchararray.add('Cluster_',numpy.char.array(range(1, len(run_clusters) + 1)).decode('UTF-8'))).assign(N=x.N).assign(
        bit_vector=x.bitStrings)
    df.to_csv(outDir + 'responsibilities.txt', sep="\t")

def hardAssignment(responsibilities):
    '''
    Assign bitStrings to clusters based on the probability distribution across clusters.
    @param: responsibilities - DataFrame consisting of the headings     Cluster_1   Cluster_2   ... Cluster_n   N   bit_vector
    @return: labels - DataFrame consisting of the heading Cluster   N   bitString, where cluster is 1,2,...n and bitString
    '''
    cluster_cols =  [col for col in responsibilities if col.startswith('Cluster')]
    probs = responsibilities[cluster_cols].as_matrix()
    # Randomly sample clusters from the probability distribution. Add 1 to make the clusters 1-indexed
    labels = numpy.array([numpy.random.choice(list(range(len(cluster_cols))), p=p) for p in probs])+1
    hardAssign = {'Cluster':labels, 'bitString':responsibilities.bit_vector, 'N': responsibilities.N}

    return pandas.DataFrame(hardAssign)

def SaveClusters(run_clusters, outFilename, refs_info):
    # Create index of data frame, containing Position, Cluster i for i in Clusters

    gene, start = refs_info.get_coordinates()
    end = start + run_clusters.shape[1]
    index = numpy.core.defchararray.add('Cluster_',numpy.char.array(range(1, len(run_clusters) + 1)).decode('UTF-8'))
    index = numpy.concatenate((['Position'], index), 0)

    # Create numpy array of positions and cluster probabilities
    dataset = numpy.asarray(run_clusters)
    position_index = numpy.expand_dims(numpy.array(range(start, end)), 1).T
    dataset = numpy.concatenate((position_index, dataset), 0)

    # Write comments:
    outF = open(outFilename, "w")
    outF.write(refs_info.text)
    outF.close()

    # Convert data set to dataframe and save to clusters text file.
    df =pandas.DataFrame(dataset, index=index).T
    df['Position'] = df['Position'].astype(int)

    df.to_csv(outFilename, sep="\t", index=False, mode="a")


def LoadClusters(filename):
    """
    Load a cluster file into a list of clusters
    :param filename: Name of the cluster file
    :return: a list of lists (clusters)
    """
    clusters = EMClass.clusterGroup(filename)
    return clusters.aslist()


def loadResp(filename):
    return EMClass.bitvect_resp(filename)


def save_cluster_prob():
    # to write
    return


# Not sure what this function is for, but shouldn't be used
def load_cluster_prob(filename):
    F = open(filename, 'r')
    cluster_prob = []
    for line in F:
        line.strip()
        cluster_prob.append(float(line))
    return cluster_prob


def save_jackpot(jackpot, filename):
    """
    Save a jackpot object to a tab sep file
    :param jackpot: Object of class jackPotting
    :param filename: Name of the file to which jackpot will be saved
    :return: void
    """
    jackpot.to_csv(filename)


def load_jackpot(filename):
    """
    Load a jackpot object from a tab sep file
    :param filename: Name of the file from which the jackpot will be loaded
    :return: An object of class jackPotting
    """
    df = pandas.read_csv(filename, sep="\t", header=0)
    return EMClass.jackPotting(df["N_theory"].values, df["N_obs"].values, df["nmuts"].values, df["bitStrings"].values)

def get_Nmut(file_name):
    """
    Get the Number of mutations column from the bitvector file
    :param file_name: path to bitvector file
    :return: pd.Series of N_mutations
    """
    df = pd.read_csv(file_name, sep='\t', skiprows=2, usecols=['N_Mutations'], index_col=False)
    return df['N_Mutations']

def MAD(y):
    """
    :param y: pandas Series
    :return: median absolute deviation
    """
    return abs(y - y.median()).median()

def get_threshold(s, m):
    """
    :param s: pandas Series
    :param m: number of MAD
    :return:
    """
    return m * MAD(s) /0.6745 + s.median()


# TODO : in labels use binary string instead of the string from bitvectors
def load_binary(binary_filename, muts_threshold=0, max_mut_prob=-1.0, mask_index=[], nreads=numpy.inf,
                reverse_column_order=False, na_threshold=0.2, nan_valid=False, deletions='0', web_version=False,
                dmuts=3):
    """

    :param binary_filename:
    :param muts_threshold:
    :param max_mut_prob: maximum number of mutation allowed in a read. 
    :param mask_index:
    :param nreads:
    :param reverse_column_order:
    :param na_threshold:
    :param dmuts:
    :param nan_valid: If True, nans around mutations are allowed and transformed to 0s. If False, the read is discarded.
    :return:
    """
    print(muts_threshold, na_threshold, deletions, dmuts)
    if type(binary_filename) is not str:
        print("Error: binary file must be a string. A %s was given" % type(binary_filename))
        sys.exit()

    sep = "\t"
    binF, header, refs_info = load_header(binary_filename)
    # TODO: change the max_mut_prob to N MAD
    if max_mut_prob == -1.0:
        Nmut = get_Nmut(binary_filename)
        Nmut_threshold = get_threshold(Nmut, 3)  # MAD/0.6745 = std of normal
    else:
        Nmut_threshold = max_mut_prob
    print("Reads with more than " + str(Nmut_threshold) + " mutations will be dropped.")

    #if len(header.strip().split(sep)) < 4:
    #    sep = " "
    #    if len(header.strip().split(sep)) < 4:
    #        print("ERROR: The binary-vectors file is neither tab nor space delimited.")
    #        sys.exit()
    if reverse_column_order:
        nmuts_i = 3
        refs_i = 2
    else:
        nmuts_i = 2
        refs_i = 5
    bit_vectors = {}
    N_occur = {}
    unique_bitStrings = []  # Used to maintain the same order each time a bitvector file is loaded
    labels = {'refs': [], 'bitStrings': []}
    bitVect_keys = {}  # This will contain modified bitStrings to be used as keys for bit_vectors and N_occurs
    # in case the original bitstring contain question marks around 1s that will be transformed into
    # 0s. At the moment it contains only problematic bitString in order to minimize memory usage at
    # the expanse of code clarity. In order to make code clearer could include all bitStrings.
    N_discard = 0
    N_invalid = 0
    for i, line in enumerate(binF):
        line = line.strip().split(sep)
        bitString = line[1]
        try:
            line[nmuts_i] = float(line[nmuts_i])
        except ValueError:
            sys.stderr.write("Error: Unable to parse line: " + str(line))
            continue
        for ind in mask_index:
            if bitString[ind] == '1':
                bitString[ind] = '0'  # We mask before checking that the read is valid, so the masked mutations are
                                       # not involved in the exclusion principle
                line[nmuts_i] -= 1
        try:
            N_occur[bitString] += 1
        except KeyError:
            try:  # Check the reason the key doesn't exist is not because of transformed question marks into 0
                N_occur[bitVect_keys[bitString]] += 1  # should not happen anymore as these reads aren't to be loaded.
            except KeyError:
                if deletions == '2' and 'D' in bitString:
                    N_discard += 1
                    continue
                bitVect, num_mutations = lib.bitString_to_bitVect(bitString, deletions)
                if (num_mutations < muts_threshold or (bitString.count('.') + bitString.count('?') +
                                                       bitString.count('N')) > len(bitString) * na_threshold):
                    N_discard += 1
                    continue
                if num_mutations > Nmut_threshold:
                    N_discard += 1
                    continue
                bitString = line[1]
                # Check read is valid. Reads with nan around 1s are thrown away.
                valid_bitVect = MutProxTest(bitVect, dmuts, nan_valid=nan_valid)
                if valid_bitVect is False:  # If bitVect is not a valid bitVect
                    N_discard += 1
                    N_invalid += 1
                    continue
                if not numpy.allclose(valid_bitVect, bitVect, equal_nan=True):
                    if nan_valid == True:
                        bitVect_keys[line[1]] = bitString
                    else:
                        print("Error: a read with '?' around a mutation sneaked in.")
                        return
                N_occur[bitString] = 1
                bit_vectors[bitString] = valid_bitVect.copy()
                unique_bitStrings.append(bitString)

        try:
            labels["refs"].append(line[refs_i])
        except IndexError:
            labels["refs"].append("ref")
        labels["bitStrings"].append(bitString)
        if len(labels["refs"]) >= nreads:
            break

    nsize = len(set(map(len, bit_vectors.values())))
    if nsize > 1:
        print("Error: The binary vectors don't all have the same size. %d different sizes detected." % nsize)
        return None, None

    x = EMClass.ind_bit_vect(bit_vectors, N_occur, bitStrings=numpy.array(unique_bitStrings), refs_info=refs_info,
                             filename=binary_filename)

    # Mask unwanted bases:
    #if len(mask_index) > 0:
    #    assert len(refs_info.refs) == 1, "Error: There is more than one reference in the header. There need to be only" \
    #                                     " one if mask_index is used, otherwise there might be more than one " \
    #                                     "coordinate offset possible."
    #    ref = list(refs_info.refs.keys())[0]
    #    mask_index = numpy.array(mask_index) - refs_info.refs[ref]['start']
    #    x.bitMatrix[mask_index] = 0

    labels["refs"] = pandas.DataFrame(numpy.asarray(labels["refs"]))
    labels["bitStrings"] = pandas.DataFrame(numpy.asarray(labels["bitStrings"]))

    binF.close()

    print("Bitvectors loaded: %d reads, %d kept, %d discarded %d invalid." % (i + 1, x.N.sum(), N_discard, N_invalid))

    if web_version:
        return x, labels, Nmut_threshold, x.N.sum(), N_discard, N_invalid
    else:
        return x, labels

def load_header(filename):
    """
    Loads the header and the comments in a file (e.g.: bitvector file). It assumes there is at least a header.
    The comment lines must start with '@'.
    :param filename: String, name of the bitVector file
    :return: file object located at the first bitvector, the header (str), and comments (object of class comments)
    """

    binF = open(filename, 'r')
    comments = EMClass.comments(binF)

    return binF, comments.header, comments


def splice_binary(binary_filename, out_filename, index_list, insep=' '):
    binF = open(binary_filename, 'r')
    outF = open(out_filename, 'w')
    header = binF.next()
    outF.write(header)
    for line in binF:
        line = line.split(insep)
        assert len(line) == 5, \
            'ERROR: The line was not parsed properly. Only %d fields were generated instead of 5.' % len(line)
        bitvect = numpy.array(list(line[1]))[index_list]
        line[1] = ''.join(bitvect)
        bitvect[bitvect == '?'] = '2'
        bitvect = bitvect.astype(float)
        bitvect[bitvect == 2.0] = numpy.nan
        line[2] = str(int(numpy.nansum(bitvect)))
        outF.write('\t'.join(line))
    binF.close()
    outF.close()


def MutProxTest(bitVect, dmuts, nan_valid=False):
    """
    Check that a bitVect doesn't have a mutation dmuts bases from another one.
    If bitVect is valid, a new bitVect is returned where nan located dmuts bases around 1s have been replaced by 0s.
    :param bitVect: numpy.array containing 0, 1 and numpy.nan
    :param dmuts: integer. Number of bases around mutations that are required to be clear of any mutations
    :param nan_valid: boolean. If False, a nan around a 1 is considered as non-valid and return False. If True, the nan is converted into a 0.
    :return: numpy.array containing 0, 1 and numpy.nan
    """
    valid_bitVect = bitVect.copy()
    for i in numpy.where(bitVect == 1.0)[0]:
        sub_bitVect = bitVect[max(0, i - dmuts):min(len(bitVect), i + dmuts + 1)]
        if nan_valid:
            if numpy.nansum(sub_bitVect) > 1.0:
                return False
            else:
                valid_bitVect[max(0, i - dmuts):min(len(bitVect), i + dmuts + 1)] = 0
                valid_bitVect[i] = 1
        else:
            if sub_bitVect.sum() > 1.0 or not numpy.isfinite(sub_bitVect.sum()):
                return False
            else:
                pass
    return valid_bitVect


def LoadSimInfo(filename):
    def clust_parse(x):
        return numpy.asarray(map(float, x.strip().split(" ")))

    inF = open(filename, 'r')
    line = inF.next()
    cluster_probs = clust_parse(inF.next())
    Nclust = len(cluster_probs)
    line = inF.next()
    line = inF.next()
    real_cluster_probs = clust_parse(inF.next())
    line = inF.next()
    line = inF.next()
    clusters = []
    for i in range(Nclust):
        clusters.append(clust_parse(inF.next()))
    clusters = numpy.asarray(clusters)
    line = inF.next()
    line = inF.next()
    real_clusters = []
    for i in range(Nclust):
        real_clusters.append(clust_parse(inF.next()))
    real_clusters = numpy.asarray(real_clusters)

    inF.close()
    return clusters, real_clusters, cluster_probs, real_cluster_probs


def get_dir(workDir, d, EM_its, max_its, prefix=''):
    return workDir + '/' + prefix + '_Nclusters_%d_Nruns_%d_maxits_%d' % (d, EM_its, max_its)


# TODO: It seems that the option ref is broken and that it will automatically choose the sequence from the last ref.
def get_sequence(fasta, coordinates, ref=None, sequence=None):
    if sequence is None:
        fasta = fasta if fasta is not None else '/lab/solexa_rouskin/references/HIVwtPhil.txt'
        F = open(fasta, 'r')
        seq = ''
        refseqs = {}
        line = F.next().strip()
        assert line[0] == ">", "ERROR: The reference fasta file must start with a reference line ('>')"
        ref = line.strip(">")
        refseqs[ref] = ''
        for line in F:
            line = line.strip()
            if line[0] == ">":
                ref = line.strip(">")
                refseqs[ref] = ''
            else:
                refseqs[ref] += line
        F.close()
        if ref is None:
            assert len(refseqs.keys()) == 1, \
                "ERROR: More than one reference are present in the fasta file. Reference needs to be specified."
            ref = refseqs.keys()[0]
        else:
            if ref not in refseqs.keys():
                print("Error: The specified reference can not be found in the fasta file.")
            else:
                pass
        if coordinates is None:
            sequence = refseqs[ref]
        else:
            sequence = refseqs[ref][coordinates[0]:coordinates[1]]
    else:
        if any(x is not None for x in [fasta, coordinates]):
            print("WARNING: -S (--sequence) is provided, -F (--fasta) and -C (--coordinates) are ignored.")
    return sequence


def write_log(runs_dir, args):
    if not os.path.exists(runs_dir):
        os.makedirs(runs_dir)
    filename = runs_dir + '/log.txt'
    F = open(filename, "w")
    F.write("Start of the runs: " + get_date() + "\n\n")
    F.write("Ran from: " + get_current_dir() + "\n\n")
    F.write("Git commit: " + get_git() + "\n\n")
    F.write("Command line: " + get_command() + "\n\n")
    F.write("full list of options:\n")
    args_str = ''
    for arg, value in vars(args).items():
        args_str += arg + " = " + str(value) + "\n"
    F.write(args_str)
    F.close()


def write_lls(all_logliks, workDir):
    # all_logliks is a list of list. Each element is the list of lls for a run.
    all_logliks = homogeneize_dataFrame(all_logliks).T
    last_logliks = all_logliks.values[all_logliks.shape[0] - 1]
    loglik_file = open(workDir + '/loglikelihoods.txt', 'w')
    loglik_file.write('Runs\tLog_likelihood\n')
    loglik_file.write('d1\t' + str(last_logliks[0]) + '\n')
    for irun, loglik in enumerate(last_logliks[1:]):
        loglik_file.write(str(irun + 1) + '\t' + str(loglik) + '\n')
    loglik_file.close()
    header = ['Run_d1']
    for i in range(all_logliks.shape[1] - 1):
        header.append('Run_' + str(i + 1))
    all_logliks.to_csv(workDir + '/all_logliks.txt', sep='\t', header=header)


def get_git():
    filepath = sys.argv[0]
    try:
        filepath = subprocess.check_output("which", filepath)
    except TypeError:
        pass

    if 'RNA_structure' not in filepath:
        filepath = subprocess.check_output(["pwd", "-P"]).strip()
    else:
        pass

    pathlist = numpy.array(filepath.split("/"))
    ind = numpy.where(pathlist == 'RNA_structure')[0][-1]
    pathlist = pathlist[:ind + 1]

    gitpath = '/'.join(pathlist)

    try:
        return str(subprocess.check_output(["git", "--git-dir=" + gitpath + "/.git/", "describe", "--always"]).strip())
    except subprocess.CalledProcessError:
        return "Could not find the git path '" + gitpath + "'."


def get_current_dir():
    return str(subprocess.check_output(["pwd", "-P"]).strip())


def get_date():
    return time.strftime("%d %b %Y, %H:%M:%S")


def get_command():
    return (' '.join(sys.argv))


def homogeneize_dataFrame(all_logliks):
    lengths = numpy.asarray([len(x) for x in all_logliks])
    lengths = max(lengths) - lengths
    for i, l in enumerate(lengths):
        all_logliks[i] = all_logliks[i] + all_logliks[i][-1:] * l
    return pandas.DataFrame(all_logliks)


def mutation_rate_histogram(resp_file, deletions):
    vals = {}
    df = pd.read_csv(resp_file, sep='\t')
    for index, row in df.iterrows():
        N = row['N']
        n_mut = lib.bitString_to_bitVect(row['bit_vector'], deletions)[1]
        if n_mut in vals:
            vals[n_mut] = vals[n_mut] + N
        else:
            vals[n_mut] = N
    out_file = os.path.join(os.path.split(resp_file)[0], "mutation_histogram.png")
    plt.bar(list(vals.keys()), vals.values())
    plt.savefig(out_file)


def find_nread(bitfile):
    """
    Find the number of reads in a bitvector file
    :param bitfile: bitvector file
    :return: int, number of reads
    """
    bitF, _, __ = load_header(bitfile)
    nline = 0
    for _ in bitF:
        nline += 1
    bitF.close()
    return nline