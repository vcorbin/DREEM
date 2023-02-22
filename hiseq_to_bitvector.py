import pysam
import numpy as np
import pandas as pd
import multiprocessing
import time

try:
    from .helper import *
    from .alignedRead import Read
    from .readQuality import is_softclipped
except:
    from helper import *
    from alignedRead import Read
    from readQuality import is_softclipped


def collect_reads(inFile, nreads, reference, ref_start, ref_stop, strict_overlap, target_length=None, save_reads=None, verbose=False,
                  web_version=True):
    '''
    Collects all the relevant reads corresponding to a region of interest, and outputs a list containing the reads.
    @param: inFile - the name of the input sam/bam file
    @param: nreads - the maximum number of reads to allow
    @param: reference -  the name of the reference gene
    @param: ref_start - the start (inclusive) of the region of interest
    @param: ref_stop - the end (exclusive) of the region of interest
    @param: strict_overlap - Boolean which is true iff we want the read to completely cover the region of interest

    @output:
        None
    '''
    sam = pysam.AlignmentFile(inFile, 'rb')

    all_reads = []

    processed_reads = 0

    coverage = {}
    softclip = 0

    for read in sam:
        if read.reference_name != reference:
            continue
              
        if read.is_unmapped: continue

        # ignore read if it doesn't overlap with the trimmed region
        if read.reference_start > ref_stop or read.reference_end < ref_start:
            continue

        if strict_overlap:
            # Only allow reads that completely cover the region of interest
            if not (read.reference_start <= ref_start and read.reference_end >= ref_stop):
                continue
        
        # Remove bad reads
        if is_softclipped(read, target_length):
            softclip += 1
            print(read.cigartuples)
            continue

        try:
            coverage[(read.reference_start, read.reference_end)] += 1
        except KeyError:
            coverage[(read.reference_start, read.reference_end)] = 1

        all_reads.append(Read(read))

        processed_reads += 1
        if processed_reads > nreads: break

    sam.close()

    print('{0} reads were discarded due to softclipping'.format(softclip))
    coverage_plot = plot_hiseq_read_coverage(coverage, web_version)

    return (all_reads, coverage, coverage_plot)


def process_unpaired_bitvector(read, pre_index, index_length, ref_trim, ref_lengths):
    '''
    #TODO: Write documentation
    '''

    mutvector, delvector, basevector, indel_vector = create_bitvectors(read, ref_trim, ref_lengths)
    bitvector = mutvector + delvector
    if pre_index != '' and index_length != 0:
        index = get_index(pre_index, index_length, read.get_original_query_sequence())
    else:
        index = ''

    n_muts = np.nansum(mutvector)
    n_dels = np.nansum(delvector)

    coverage = len(np.where(basevector != '.')[0])

    return read.query_name, bitvector, mutvector, delvector, ''.join(
        basevector),''.join(indel_vector), n_muts, n_dels, coverage, index, read.mutated_bases
