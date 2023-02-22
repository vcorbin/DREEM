#--------------------- Helper functions for parsing MD string and converting positions to quality ------------------------------------
import pysam
import numpy as np
import pandas as pd
import multiprocessing
import time 

def is_softclipped(read, insert_length):
	'''
	Quality control measure which determines whether a read is too softclipped to be considered a proper alignment.
	What's the rule for discarding these?

	If the insert length is greater than 300, ensure the read is no more than 20% softclipped
	Otherwise, do nothing

	@param: read - pysam read 
	@param: insert_length - The length of the original target sequence, which might be padded for sequencing. None if unknown
	
	@return: True if the read is softclipped and should be discarded, False otherwise
	'''

	REQUIRED_MATCH_PERCENT = 0.8 # We require 80% of the insert_length to NOT be soft clipped
	MIN_INSERT_LENGTH = 100

	# If no insert length is provided, default behaviour is to require MATCH_PERCENT of the entire read.
	# Length of entire read is inferred by pysam from the CIGAR string (includes hard clipped bases)
	if insert_length is None:
		insert_length = read.infer_read_length()


	# If the insert length is less than some minimum target length, do not do any filtering.
	if insert_length<MIN_INSERT_LENGTH:
		return False

	matches = 0
	clips = 0
	for (op, count) in read.cigartuples:
		if op == 4:
			clips += count
		if op == 0:
			matches += count

	# The number of matches must be a significant proportion of the insert_length
	# If the number of matches is fewer than this amount, the read is soft clipped
	if matches < insert_length*REQUIRED_MATCH_PERCENT:
		return True
		
	return False