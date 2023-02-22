#!/usr/bin/python

import sys

sys.path.append('..')
import os
import argparse
import numpy

try:
    from helper import *
except:
    from .helper import *

try:
    from coordinates import get_coordinates
except:
    from .coordinates import get_coordinates

try:
    from coordinates import load_fasta
except:
    from .coordinates import load_fasta

try:
    from .EM_classes import split_filename
except:
    from EM_classes import split_filename

"""
TODO:
-Plot read coverage as png file, also?
-Save everything in a new directory?
-Get a better consensus sequence?
-Number of mutations histogram/ coverage
"""


def convert_to_bitvector(inFile, outFile, ref_trim, paired=False, pre_index='', index_length=0, target_length=None,
                         nreads=np.inf,
                         use_multiprocessing=True, strict_overlap=False, save_reads=None, verbose=False, ref_seq=None,
                         write_log=None):
    '''
    Generates a bit vector with for the target gene in the region of interest and saves bitvector to outFile.

    @param: inFile - String which is the input sam/bam file
    @param: outFile - String which is the output file name
    @param: ref_trim - Boolean mapping the target gene to target region, gene: [start,stop)
    @param: pre_index - The sequence preceding the 'unique' index
    @param: index_length - The length of the sequence of the 'unique' index
    @param: target_length - The length of the original target sequence, which might be padded for sequencing. None if unknown
    @param: nreads - maximum number of lines to read
    @param: multiprocessing - Boolean which is true if we should use multiple processors
    @param: strict_overlap - Boolean which is true iff we want both reads to completely overlap in the region of interest
    @param: save_reads - a tuple containing the coordinates such that if a read overlaps these coordinates, the read is saved. None if save no reads.
    @param: ref_seq - UNUSED
    @param: write_log - A function for writing comments to a log, otherwise None

    

    @return: dataframe containing bitvectors, titled outFile.
    '''

    # Ensure logging function works and is backward-compatible
    if write_log is None:
        def display(msg):
            print(msg)
    else:
        try:
            write_log('Begin convert to bitvector')

            def display(msg):
                write_log(msg)
                print(msg)
                return
        except:
            def display(msg):
                print(msg)

    if paired:  # parse as a Miseq (paired read) file
        try:
            from miseq_to_bitvector import collect_pairs as collect_reads
            from miseq_to_bitvector import process_paired_bitvector as process_bitvector
        except:
            from .miseq_to_bitvector import collect_pairs as collect_reads
            from .miseq_to_bitvector import process_paired_bitvector as process_bitvector
        pair_name = 'pairs of reads'
    else:  # parse as a Miseq (unpaired read) file
        try:
            from hiseq_to_bitvector import collect_reads
            from hiseq_to_bitvector import process_unpaired_bitvector as process_bitvector
        except:
            from .hiseq_to_bitvector import collect_reads
            from .hiseq_to_bitvector import process_unpaired_bitvector as process_bitvector
        pair_name = 'reads'

    sam = pysam.AlignmentFile(inFile)
    ref_lengths = dict(zip(sam.references, sam.lengths))
    sam.close()

    plots_prefix = os.path.splitext(outFile)[0]

    # WARNING: in case of multiple reference, this look like it will create multiple bitvector files with the same names,
    # erasing the previous ones. Same with the figures.

    for reference in ref_trim.keys():
        ref_start = ref_trim[reference][0]
        ref_stop = ref_trim[reference][1]
        display('Collecting matching {3} aligning to "{0}" in the region of interest {1} - {2}...'.format(reference,
                                                                                                          ref_start,
                                                                                                          ref_stop,
                                                                                                          pair_name))
        reads, coverage, coverage_plot_html = collect_reads(inFile, nreads, reference, ref_start, ref_stop,
                                                            strict_overlap,
                                                            save_reads=save_reads, verbose=verbose,
                                                            web_version=plots_prefix + '_read_coverage.html',
                                                            target_length=target_length)
        print("Read collected")
        if len(reads) == 0:
            display('No {0} found matching the region of interest. Terminating process.'.format(pair_name))
            return 0, None, None, None, coverage_plot_html
        display('Processing {0} {1}...'.format(len(reads), pair_name))

        if not use_multiprocessing:  # do not use multiprocessing
            bitvectors, all_mutated_bases = generate_bitvectors(process_bitvector, reads,
                                                                {reference: [ref_start, ref_stop]}, ref_lengths,
                                                                pre_index, index_length)
        else:  # use multiprocessing
            import multiprocessing

            nprocess = multiprocessing.cpu_count()
            mprocesses = []
            pool = multiprocessing.Pool(processes=nprocess)

            subpair_size = len(reads) // (nprocess) + 1
            display('Creating bitvectors in parallel with {0} processes'.format(nprocess))
            # Break lists into sublists for parallel processing
            for its in range(nprocess):
                subpairs = reads[subpair_size * its:min(subpair_size * (its + 1), len(reads))]
                mprocesses.append(pool.apply_async(generate_bitvectors, args=(
                    process_bitvector, subpairs, {reference: [ref_start, ref_stop]}, ref_lengths, pre_index,
                    index_length)))
            pool.close()
            # Wait for results
            display('Waiting for results...')
            pool.join()
            display('Done with multiprocessing. Fetching results...')
            bitvector_frames = []
            all_mutated_bases = []
            for its, process in enumerate(mprocesses):
                bitvector_df, mutated_bases = process.get()
                bitvector_frames.append(bitvector_df)
                all_mutated_bases += mutated_bases
            bitvectors = pd.concat(bitvector_frames, ignore_index=True)

        # Generate consensus sequence based on all reads studied

        sequence, DMS_changes_html = plot_DMS_changes(all_mutated_bases, ref_start, ref_stop,
                                                      web_version=plots_prefix + '_DMS_mutations.html')

        if ref_seq is None:
            ref_seq = ''.join(sequence)

        # Plot population average
        pop_avg_html = plot_population_average(bitvectors,
                                               columns=['Combined_vector', 'Mutations_vector', 'Deletions_vector'],
                                               ref_start=ref_start, ref_stop=ref_stop,
                                               web_version=plots_prefix + '_populationAverage.html',
                                               sequence=list(ref_seq))

        # Plot mutations per base
        mut_per_base_html, AC_percent = plot_muts_per_base(bitvectors['Mutations_vector'], list(ref_seq),
                                                           web_version=plots_prefix + '_muts_per_base.html')

        # Plot mutation per read distribution
        basevectors = bitvectors[
            ['Query_name', 'Bases_vector', 'N_mutations', 'N_deletions', 'Coverage', 'Reference', 'Index']]
        mut_dist_html = plot_mutation_distribution(basevectors,
                                                   web_version=plots_prefix + '_mutation_distribution.html')
        display('Saving output file to {0}.'.format(outFile))

        ref_length = ref_stop - ref_start
        ref_index = '{0}_{1}'.format(reference, ref_start)

        header = ''
        header += '@start\t{0}\t{1}\n'.format(ref_index, ref_start)  # Add the region start
        header += '@length\t{0}\t{1}\n'.format(ref_index, ref_length)  # Add the region length
        header += '@seq\t{0}\t{1}\n'.format(ref_index, ref_seq)  # Add the region sequence

        # TODO
        # directory = 'analysis_'+outFile
        # if not os.path.exists(directory):
        #     os.makedirs(directory)

        bases_file = open(outFile, 'w')
        bases_file.write(header)
        basevectors.to_csv(bases_file, sep="\t", index=False,
                           header=['Query_name', 'Bases_vector', 'N_mutations', 'N_deletions', 'Coverage', 'Reference',
                                   'Index'])
        bases_file.close()

        # Save indelvectors:
        indelvectors = bitvectors[
            ['Query_name', 'Indels_vector', 'N_mutations', 'N_deletions', 'Coverage', 'Reference', 'Index']]

        indels_file = open(outFile + '_indels.txt', 'w')
        indels_file.write(header)
        indelvectors.to_csv(indels_file, sep="\t", index=False,
                            header=['Query_name', 'Indels_vector', 'N_mutations', 'N_deletions', 'Coverage',
                                    'Reference',
                                    'Index'])

        # Plot the indel signal
        indels_html = plot_indels(indelvectors, start=ref_start, stop=ref_stop, sequence=list(ref_seq),
                                  web_version=plots_prefix + '_indels_distribution.html')

        # Convert bitvectors to strings to save as output

        def stringify_column(array, whichvector='bit'):
            new_array = []
            for i, vector in enumerate(array):
                vector = np.array(list(vector))
                if whichvector == 'bit' or whichvector == 'mut':
                    vector[vector == 'A'] = '1'
                    vector[vector == 'T'] = '1'
                    vector[vector == 'G'] = '1'
                    vector[vector == 'C'] = '1'
                if whichvector == 'bit' or whichvector == 'del':
                    vector[vector == 'D'] = '1'
                new_array.append(''.join(vector))
            return new_array

        display('Done.')

        # TODO: Replace None with muts_per_base_html
        return len(
            reads), DMS_changes_html, pop_avg_html, mut_dist_html, mut_per_base_html, coverage_plot_html, AC_percent


def generate_bitvectors(process_bitvector, reads_pairs, ref_trim, ref_lengths, pre_index, index_length):
    '''
    Generate bitvectors from two lists of Miseq reads, where pairs[i] is a tuple of paired alignedRead objects.
    Mutations are denoted by a '1', non-mutations by a '0' and poor quality positions by '?'
    A mutation is valid only if it occurs in both read, where there is overlap.

    @param: read_pairs - list of Miseq read tuples. Each read is a AlignedRead object
    @param: ref_trim - Boolean mapping the target gene to target region, gene: [start,stop)
    @param: ref_lengths -
    @param: pre_index - The sequence preceding the 'unique' index
    @param: index_length - The length of the sequence of the 'unique' index

    @return: a dataframe containing the query name, bit vector, number of mutatons, reference name, start position
    '''

    bitvectors = {'Query_name': [],
                  'Combined_vector': [], 'Mutations_vector': [], 'Bases_vector': [], 'Indels_vector': [],
                  'Deletions_vector': [],
                  'N_mutations': [], 'N_deletions': [],
                  'Coverage': [],
                  'Reference': [],
                  'Index': []}

    all_mutated_bases = []
    ref_name = list(ref_trim.keys())[0]
    ref_start = ref_trim[ref_name][0]
    for pair in reads_pairs:
        query_name, bitvector, mutvector, delvector, basevector, indelvector, n_muts, n_dels, coverage, index, mutated_bases = process_bitvector(
            pair, pre_index, index_length, ref_trim, ref_lengths)
        all_mutated_bases += mutated_bases

        bitvectors['Query_name'].append(query_name)
        bitvectors['Combined_vector'].append(bitvector)
        bitvectors['Mutations_vector'].append(mutvector)
        bitvectors['Deletions_vector'].append(delvector)
        bitvectors['Bases_vector'].append(basevector)
        bitvectors['Indels_vector'].append(indelvector)
        bitvectors['N_mutations'].append(n_muts)
        bitvectors['N_deletions'].append(n_dels)
        bitvectors['Coverage'].append(coverage)
        bitvectors['Reference'].append('{0}_{1}'.format(ref_name, ref_start))
        bitvectors['Index'].append(index)
    bitvectors = pd.DataFrame(bitvectors,
                              columns=['Query_name', 'Combined_vector', 'Mutations_vector', 'Deletions_vector',
                                       'Bases_vector', 'Indels_vector', 'N_mutations', 'N_deletions', 'Coverage',
                                       'Index', 'Reference'])

    return bitvectors, all_mutated_bases


def __get_ref_info_from_args(ref_coord, primers, reference, fasta_ref):
    if ref_coord is None:
        if primers is None or fasta_ref is None:
            display(
                "Error: When coordinates are not provided, a set of primers (-p) and a reference fasta file (-r) must "
                "be specified")
            sys.exit()
        else:
            # Fasta_ref can only have one sequence in it, it will be checked when calculating the coordinates
            primer_1, primer_2 = primers.strip().split(",")
            ref_coord, ref_seq, ref_name = get_coordinates(primer_1, primer_2, fasta_ref, buffering=buffer)
            ref_seq = ref_seq[ref_coord[0]:ref_coord[1]]
            if reference is not None:
                if reference != ref_name:
                    display(
                        "Warning: The reference name (--ref_name) provided doesn't match the reference found in the "
                        "reference fasta file (-r). The name provided '%s' will be ignored and the reference from the"
                        " file '%s' will be used instead" % (reference, ref_name))
                else:
                    pass
            else:
                reference = ref_name

    else:
        ref_coord = [int(x) for x in ref_coord.strip().split(",")]
        # Check a name for the reference can be found
        if fasta_ref is None:
            if reference is None:
                display(
                    "Error: Neither a reference fasta file (-r) nor a reference name (--ref_name) have been specified.")
                sys.exit()
            else:
                ref_seq = None
        else:
            ref_fasta_content = load_fasta(fasta_ref)
            if reference is None:
                assert len(ref_fasta_content.keys()) == 1, "Error: the fasta file %s contains %d references, " \
                                                           "when it should contain only one" % (
                                                               fasta_ref, len(ref_fasta_content.keys()))
                reference = list(ref_fasta_content.keys())[0]
                seq = ref_fasta_content[reference]
                assert len(seq) > 0, "Error: there is no sequence corresponding to the name %s" % reference
            else:
                if reference not in ref_fasta_content:
                    if len(ref_fasta_content.keys()) != 1:
                        display("Error: The reference provided by --ref_name (%s) was not found in the reference fasta "
                                "file %s" % (reference, fasta_ref))
                        sys.exit()
                    else:
                        display(
                            "Warning: The reference name (--ref_name) provided doesn't match the reference found in the "
                            "reference fasta file (-r). The name provided '%s' will be ignored and the reference from the"
                            " file '%s' will be used instead" % (reference, list(ref_fasta_content.keys())[0]))
                        reference = list(ref_fasta_content.keys())[0]
            ref_seq = ref_fasta_content[reference]
            assert len(ref_seq) >= ref_coord[1], "Error: the coordinates specified by -c (%d, %d) don't match the " \
                                                 "reference sequence length (%d)" \
                                                 % (ref_coord[0], ref_coord[1], len(ref_seq))
            ref_seq = ref_seq[ref_coord[0]:ref_coord[1]]

    return ref_coord, ref_seq, reference


def make_bitvector_filename(bam_file, ref_coord):
    """
    Make the name of a bitvector file from the name of the bam file and the reference coordinates
    :param bam_file: str, name of the file used to create the bitvector
    :param ref_coord: list of 2 integers, coordinates (with respect to the reference sequence used in the alignment) of 
    the bitvectors that will be written in the bitvector file
    :return: the name of the bitvector file (with same path as the bam_file and .txt extension)
    """
    ref_str = '_'.join([str(ref_coord[0]), str(ref_coord[1])])
    prefix = split_filename(bam_file).noext
    prefix = '{0}_{1}'.format(prefix, ref_str)
    outFile = '{0}_{1}'.format(prefix, 'bitvectors.txt')
    return outFile


########################################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert aligned DMS-modified reads to bitvectors')

    parser.add_argument('-f', '--input_file', type=str, dest='inFile', required=True,
                        help='String, required. File containing aligned reads to be converted to bitvectors')

    strict_parser = parser.add_mutually_exclusive_group(required=False)
    strict_parser.add_argument('-so', '--overlap', dest='strict_overlap', action='store_true')
    parser.set_defaults(strict_overlap=False)

    parser.add_argument('-r', '--fasta_ref', type=str, dest='fasta_ref', default=None,
                        help='Name of the reference fasta file. If not specified, a reference name (--ref_name) and '
                             'a set of coordinates (-c) need to be input. As of now only one reference must be present '
                             'in the fasta.')

    parser.add_argument('--ref_name', type=str, dest='reference', default=None,
                        help='String, required if -r is not specified. Enter the name of the gene/reference under '
                             'consideration.')

    parser.add_argument('-c', '--coord', type=str, dest='ref_coord', default=None,
                        help='Comma delimited integer pair, required if -r is not specified. Start (inclusive) and end '
                             '(exclusive) positions of region of interest, delimited by comma, eg -c "500,600". If '
                             'provided it will ignore any primer sequence given')

    parser.add_argument('-p', '--primers', type=str, dest='primers', default=None,
                        help='sequence of two primers separated by ",". required if --coord is not specified')

    parser.add_argument('--paired', type=int, dest='paired', default=1,
                        help='Integer (0,1), default: 1 (paired). Select "0" to parse as a Hiseq (unpaired reads) file, or "1" to '
                             'parse as Miseq (paired read) file for paired reads. Reads must be sorted by pairs for a '
                             'Miseq file.')

    parser.add_argument('-i', '--pre-index', type=str, dest='pre_index', default='',
                        help='String, default: "". Enter the sequence which precedes the index, or omit if there is no '
                             'index.')

    parser.add_argument('-l', '--index_length', type=int, dest='index_length', default=0,
                        help='Integer, default: 0. Enter the length of the index, or omit if there is no index.')

    parser.add_argument('-t', '--target_length', type=int, dest='target_length', default=None,
                        help='Integer, default: Inferred from read. Enter the length of the insert to filter for reads which are on target and are not excessively soft-clipped.')

    parser.add_argument('-o', '--output_file', type=str, dest='outFile', default=None,
                        help='String, default: None. Name of the bit-vector file if one is created. If not specified, '
                             'the name of the bam file is used.')

    parser.add_argument('-m', '--multiprocessing', type=int, dest='use_multiprocessing', default=1,
                        help='Integer (0,1), default: 1. Enter 1 to use multiprocessing to speed up code, or 0 for no '
                             'parallelization.')

    parser.add_argument('-n', '--nreads', type=str, dest='nreads', default='1000000',
                        help='Integer or "all", default: 1 million. Maximum number of reads to analyse.')

    parser.add_argument('-v', '--verbose', type=int, dest='verbose', default=0,
                        help='Integer (0,1), default: 0. Enter 1 to output extra details to shell.')

    args = parser.parse_args()
    print(args)

    # --------------------------------------------- Variables conversion ----------------------------------
    inFile = args.inFile
    paired = args.paired
    fasta_ref = args.fasta_ref
    primers = args.primers
    reference = args.reference
    ref_coord = args.ref_coord
    outFile = args.outFile
    index_length = args.index_length
    pre_index = args.pre_index
    use_multiprocessing = args.use_multiprocessing
    if args.nreads == "all":
        nreads = numpy.inf
    else:
        nreads = int(args.nreads)
    verbose = args.verbose
    strict_overlap = args.strict_overlap
    target_length = args.target_length

    buffer = 3

    # Check all options required are there
    # Check the coord can be found
    ref_coord, ref_seq, reference = __get_ref_info_from_args(ref_coord, primers, reference, fasta_ref)
    print(reference, ref_seq)
    ref_trim = dict()
    ref_trim[reference] = ref_coord

    if outFile is None:
        outFile = make_bitvector_filename(inFile, ref_coord)

    convert_to_bitvector(inFile, outFile, ref_trim, paired=paired, pre_index=pre_index, index_length=index_length,
                         target_length=target_length,
                         nreads=nreads, use_multiprocessing=use_multiprocessing, strict_overlap=strict_overlap,
                         ref_seq=ref_seq)


################################################################################################ TESTS ################################################################################################################################
def small_test():
    inFile = 'small_sample.sam'
    outFile = 'small_sample.bv'
    print('Testing small file')
    convert_to_bitvector(inFile, outFile, ref_trim={'ref_gene': [0, 100]}, paired=True, pre_index='', index_length=0,
                         nreads=1000, use_multiprocessing=False, strict_overlap=False, save_reads=None, verbose=False)


#small_test()


def paired_test():
    inFile = 'test_files/170117Rou_D17-802Jan2017l2.sam'
    outFile = 'test_files/miseq_test'
    display('Testing paired file')
    convert_to_bitvector(inFile, outFile, ref_trim={'HIV': [455, 855]}, paired=True, pre_index='', index_length=0,
                         nreads=1000, use_multiprocessing=False, strict_overlap=False, save_reads=None, verbose=False)


def paired_index_test():
    inFile = 'test_files/index_BAM.txt'
    outFile = 'test_files/miseq_index_test'
    display('Testing paired index file')
    convert_to_bitvector(inFile, outFile, ref_trim={'VegF': [50, 250]}, paired=True, pre_index='GAGTCCCAGGTCCGACGATCC',
                         index_length=10, nreads=10, use_multiprocessing=True, strict_overlap=False, save_reads=None,
                         verbose=False)


def unpaired_test():
    inFile = 'test_files/170117Rou_D17-802Jan2017l2.sam'
    outFile = 'test_files/hiseq_test'
    display('Testing unpaired file')
    convert_to_bitvector(inFile, outFile, ref_trim={'HIV': [455, 555]}, paired=False, pre_index='', index_length=0,
                         nreads=1000, use_multiprocessing=False, strict_overlap=False, save_reads=None, verbose=False)


def unpaired_index_test():
    inFile = 'test_files/index_BAM.txt'
    outFile = 'test_files/hiseq_index_test'
    display('Testing unpaired index file')
    convert_to_bitvector(inFile, outFile, ref_trim={'VegF': [50, 250]}, paired=False, pre_index='GAGTCCCAGGTCCGACGATCC',
                         index_length=10, nreads=10, use_multiprocessing=True, strict_overlap=False, save_reads=None,
                         verbose=False)
