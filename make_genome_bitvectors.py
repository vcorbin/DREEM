from make_bitvector import convert_to_bitvector
import numpy as np
import argparse

def make_genome_bitvector(fastaFile, samFile, paired=True,target_length=None, regionsFile=None, nreads=np.inf,
							 use_multiprocessing=True, strict_overlap=False, verbose=False):
	'''
	@param: fastaFile - FASTA file containing the names and transcripts of the genes you want to make bitvectors for
	@param: samFile - SAM file containing the aligned reads of the sample
	@param: regionsFile - Tab separated file, each line containing custom regions of interest of the format: Gene	Start	End
	@param: verbose - Display additional comments as code progresses
	'''

	import time 
	import os 

	#Change current path to path of sample
	if '/' in samFile:
		path = '/'.join(samFile.split('/')[:-2])
	print(path)
	# os.chdir(path)

	start_time = time.time()

	sampleName = samFile[:-4]
	i = 0
	prefix = sampleName
	while os.path.exists(prefix):
		i += 1
		prefix = '{0}_{1}'.format(sampleName,i)
	os.mkdir(prefix)

	# Helper function for writing to a log
	def write_log(msg, start_time=start_time, prefix=prefix):
		log = open('{0}/log.txt'.format(prefix),'a+')
		log.write('[{0:.5}]:{1}\n'.format(time.time()-start_time, msg))
		log.close()

	# Generate a dictionary of genes: sequence from the sequence FASTA file
	sequences = {}
	regions = {} #Create a dictionary of each gene to the region of interest
	with open(fastaFile, 'r+') as f:
		for line in f:
			line = line.strip('\n')
			if len(line) == 0:
				continue
			if line[0] ==  '>':
				name = line[1:]
			else:
				if name is None:
					assert False, 'ERROR reading fastaFile. Remove the line break on/before: \n {0}.'.format(line)
				sequences[name] = line.upper() 
				regions[name] = [0,len(line)] #Create a bitvector from the entire region
				name = None 
	# If custom regions of interest have been defined, ignore the fastaFile regions and use those instead
	if regionsFile is not None:
		try:
			regions = {}
			with open(regionsFile,'r') as f:
				for line in f:
					# line = gene \t start \t end 
					line = line.split('\t')
					regions[line[0]] = [int(line[1]), int(line[2])]
		except:
			print('Could not read regions from input file. Using entire genomic sequence as region of interest.')

	# For each gene in the fasta file, create a new directory and run make_bitvector
	for gene in sequences:
		os.mkdir('{0}/{1}'.format(prefix,gene))
		write_log('Start make_bitvector on {0}'.format(gene))

		start,end = regions[gene][0],regions[gene][1]
		ref_trim = {gene:[start,end]}
		ref_seq = sequences[gene][start:end]
		outFile = '{0}/{1}/{1}_bitvector.bv'.format(prefix,gene) 
		convert_to_bitvector(samFile, 
							outFile=outFile,
							ref_trim=ref_trim, 
							ref_seq=ref_seq, 
							paired=paired, 
							target_length=target_length, 
							nreads=nreads,
							use_multiprocessing=use_multiprocessing, 
							strict_overlap=strict_overlap, 
							write_log=write_log,
							verbose=verbose)
		write_log('Finished make_bitvector on {0}.'.format(gene))

	# Save csv file containing all the genes, their population averages and their mutation percentage on each original base
	write_log('Creating summary file of all genes.')
	for gene in sequences:
		# Read 
		pass 
	# 

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Create bitvectors from aligned DMS-modified reads for whole genome samples.')

	parser.add_argument('-f', '--fastaFile', type=str, dest='fastaFile', required=True,
						help='String, required. Fasta file containing reference sequences for genes in sample.')

	parser.add_argument('-s', '--samFile', type=str, dest='samFile', required=True,
						help='String, required. File containing aligned reads to be converted to bitvectors')


	parser.add_argument('-m', '--multiprocessing', type=int, dest='use_multiprocessing', default=1,
						help='Integer (0,1), default: 1. Enter 1 to use multiprocessing to speed up code, or 0 for no '
							 'parallelization.')
	
	parser.add_argument('-p', '--paired', type=int, dest='paired', default=1,
					help='Integer (0,1), required. Select "0" to parse as a Hiseq (unpaired reads) file, or "1" to '
						 'parse as Miseq (paired read) file for paired reads. Reads must be sorted by pairs for a '
						 'Miseq file.')

	parser.add_argument('-t', '--target_length', type=int, dest='target_length', default=None,
					help='Integer, required. Enter the length of the insert to filter for reads which are on target and are not excessively soft-clipped.')

	parser.add_argument('-n', '--nreads', type=int, dest='nreads', default=np.inf,
					help='Integer, required. Maximum number of reads to use to make bitvector.')

	args = parser.parse_args()


	make_genome_bitvector(args.fastaFile, args.samFile, paired=args.paired, target_length=args.target_length, regionsFile=None, nreads=args.nreads,
							 use_multiprocessing=args.use_multiprocessing, strict_overlap=False, verbose=False)

