import pysam
import numpy as np
import pandas as pd
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

try:
	from .helper import *
	from .alignedRead import Read
	from .readQuality import is_softclipped
except:
	from helper import *
	from alignedRead import Read
	from readQuality import is_softclipped



def collect_pairs(inFile, nreads, reference, ref_start, ref_end, strict_overlap, target_length=None, save_reads=None, verbose=False,
				  web_version=True):
	'''
	Given a Miseq aligned sam/bam file, this method collects reads and their corresponding pairs
	which align to a reference within a start and end frame.

	@param: inFile - the name of the input sam/bam file
	@param: nreads - the maximum number of reads to allow
	@param: reference -  the name of the reference gene
	@param: ref_start - the start (inclusive) of the region of interest
	@param: ref_end - the end (exclusive) of the region of interest
	@param: strict_overlap - Boolean which is true iff we want the read to completely cover the region of interest
	@param: target_length - 

	'''
	sam = pysam.AlignmentFile(inFile, 'r')

	# Initialize list containing corresponding pairs of reads.
	pairs = []

	# Keep track of unpaired reads
	unpaired_read_count = 0

	i = 0  # initialize count to get read1
	read1, read2 = None, None
	num_pairs = 0

	read_starts = {}
	read_ends = {}

	coverage = {}

	softclip = 0
	naccept = 0
	for nread, read in enumerate(sam):
		if i % 2 == 0:
			read1 = read
			i += 1  # Got read1, next time, grab read2
			continue

		## Catch read2 this time ##

		# Ensure successive reads are proper pairs
		if read.query_name.split('/')[0] != read1.query_name.split('/')[0]:  # if read does not match previous read, then discard, start afresh and grab a read2 next time
			read1 = read
			i = 1  # got read1, next time, grab read2
			if read1.reference_name == reference:
				unpaired_read_count += 1
			continue

		read2 = read
		i = 0  # got read2, next time, grab read1

		# grab two paired reads and compare them

		try:
			assert read1.reference_name == read2.reference_name  # Paired reads must be from the same reference
			# Remove reads with the wrong reference
			if read1.reference_name != reference:
				continue  # Ignore reads which do not align to the desired reference.

		except AssertionError:
			if verbose:
				print('Paired reads align to different references: %s and %s \t %s \t %s \t %s' % (
					read1.reference_name, read2.reference_name, read1.query_name.split('/')[0], read1.reference_start,
					read2.reference_start))
			continue

		try:
			assert read1.query_name.split('/')[0] == read2.query_name.split('/')[0]
		except AssertionError:
			if verbose:
				print('Unpaired reads slipped through filter')


		# Perform quality control:

		# Check if either read is soft-clipped.
		if is_softclipped(read1, target_length) or is_softclipped(read2, target_length):
			softclip += 1
			continue

		try:
			coverage[(read1.reference_start, read2.reference_start, read1.reference_end, read2.reference_end)] += 1
		except KeyError:
			coverage[(read1.reference_start, read2.reference_start, read1.reference_end, read2.reference_end)] = 1

		# Ignore read pairs if neither read aligns to the desired region on the reference.
		if ((read1.reference_start > ref_end or read1.reference_end < ref_start) and
				(read2.reference_start > ref_end or read2.reference_end < ref_start)):
			continue

		if strict_overlap:
			# Ensure that both reads entirely cover the region of interest
			if not (read1.reference_start <= ref_start and read1.reference_end >= ref_end and read2.reference_start <= ref_start and read2.reference_end >= ref_end):
				continue
		naccept += 1
		print(nread, naccept)
		pairs.append([Read(read1), Read(read2)])

		num_pairs += 1

		if num_pairs > nreads:
			break

	sam.close()
	print('{0} reads were discarded due to softclipping'.format(softclip))
	coverage_plot = plot_miseq_read_coverage(coverage, web_version)
	print('{0} unpaired reads found'.format(unpaired_read_count))
	return (pairs, coverage, coverage_plot)


def combine_pair_xor(read1, read2, ref_trim, ref_lengths, verbose=False):
	'''
	Return a bitvector which is the bitvector in agreement between two paired reads.
	@param: read1 -
	@param: read2 -

	'''
	# Create binary vectors
	mutvector1, delvector1, basevector1, indelvector1 = create_bitvectors(read1, ref_trim, ref_lengths)
	bitvector1 = mutvector1 + delvector1

	mutvector2, delvector2, basevector2, indelvector2 = create_bitvectors(read2, ref_trim, ref_lengths)
	bitvector2 = mutvector2 + delvector2

	# mutations and deletions must be present in both reads.
	nan_vals = np.intersect1d(np.argwhere(np.isnan(mutvector1)), np.argwhere(np.isnan(mutvector2)))

	mutvector = np.logical_and(mutvector1, mutvector2).astype(np.float)
	mutvector[nan_vals] = np.nan

	delvector = np.logical_and(delvector1, delvector2).astype(np.float)
	delvector[nan_vals] = np.nan

	bitvector = np.logical_and(bitvector1, bitvector2).astype(np.float)
	bitvector[nan_vals] = np.nan

	basevector = np.array(['.'] * len(basevector1))
	indelvector = np.array(['.'] * len(indelvector1))

	string_vectors = [(basevector1, basevector2), (indelvector1, indelvector2)]
	for i, vector in enumerate([basevector, indelvector]):
		# Determine regions where neither read provides coverage.
		uncovered = np.intersect1d(np.where(string_vectors[i][0] == '.'), np.where(string_vectors[i][1] == '.'))

		string_vectors[i][0][np.where(string_vectors[i][0] == '.')] = '?'
		string_vectors[i][1][np.where(string_vectors[i][1] == '.')] = '?'

		vector[np.where(string_vectors[i][0] == string_vectors[i][1])] = string_vectors[i][0][np.where(string_vectors[i][0] == string_vectors[i][1])]

		# In the case of poor quality (?) use the information from the other read. If both are low quality, this puts a ?
		vector[np.where(string_vectors[i][0] == '?')] = string_vectors[i][1][np.where(string_vectors[i][0] == '?')]
		vector[np.where(string_vectors[i][1] == '?')] = string_vectors[i][0][np.where(string_vectors[i][1] == '?')]

		# Fill in the positions where only one read covers the region
		vector[uncovered] = '.'


	return (bitvector, mutvector, delvector, basevector, indelvector, read1.mutated_bases + read2.mutated_bases)


def process_paired_bitvector(pair, pre_index, index_length, ref_trim, ref_lengths):
	'''
	#TODO: Write documentation
	'''
	(read1, read2) = pair
	# Ensure both reads come from the same reference
	assert read1.reference_name == read2.reference_name, 'Both paired reads do not come from the same reference'
	(bitvector, mutvector, delvector, basevector, indelvector, mutated_bases) = combine_pair_xor(read1, read2, ref_trim, ref_lengths)

	read1_index = get_index(pre_index, index_length, read1.get_original_query_sequence())
	read2_index = get_index(pre_index, index_length, read2.get_original_query_sequence())

	# TODO: Get a consensus index from both reads.
	if read1_index == '':
		index = read2_index
	elif read2_index == '':
		index = read1_index
	else:
		index = read1_index
	# TODO: Get a consensus index from both reads.

	coverage = len(np.where(basevector != '.')[0])
	n_muts = np.nansum(mutvector)
	n_dels = np.nansum(delvector)

	return (read1.query_name, bitvector, mutvector, delvector, ''.join(basevector), ''.join(indelvector), n_muts, n_dels, coverage, index,
			mutated_bases)
