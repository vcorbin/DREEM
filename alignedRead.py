import numpy as np

class Read:
	'''
	A copy of the arguments from a AlignedSegment object relevant to this analysis. Add attributes of mutations, mutation_vector and deletions_vector
	mutated_bases is a list of tuples (position on reference, original base, mutated base)

	Attributes:
	self.reference_name -  the name of the reference in the BAM file
	self.reference_start - the start position of the read (relative to the reference)
	self.reference_end - the end position of the read (relative to the reference)
	self.MD -  the MD tag of the read
	self.cigartuples - a list (op, length) that constitute the read. Refer to pysam documentation

	'''
	def __init__(self, pysam_read):
		self.reference_name = pysam_read.reference_name
		self.reference_start = pysam_read.reference_start 
		self.reference_end = pysam_read.reference_end
		self.MD = pysam_read.get_tag("MD")
		self.cigartuples = pysam_read.cigartuples

		self.query_name = pysam_read.query_name
		self.query_sequence = np.array(list(pysam_read.query_sequence))
		self.query_alignment_length = pysam_read.query_alignment_length
		self.query_alignment_qualities = pysam_read.query_alignment_qualities

		self.original_query_sequence = pysam_read.query_sequence
		self.aligned_pairs = pysam_read.get_aligned_pairs(with_seq=True)
		
		# Generate a list of cigar values from cigar tuple excluding softclips (4) eg [(4,5),(2,1),[1,1]]->[2,1]
		# 0 - alignment (match or mismatch), 1 - insertion, 2 - deletion, 3 - skipped region in reference, 4- softclip 
		self.cigarlist = np.hstack(np.array([[tup[0]]*tup[1] for tup in self.cigartuples if tup[0]!=4]))
		
		# Compute positions of insertions. Insert positions are calculated as the position immediately preceding the insertion
		INSERT_INDEX = 1
		self.insert_positions = np.where(self.cigarlist==INSERT_INDEX)[0]
		# Subtract the index+1 of each insert, since each insert adds an artificial increase in the index e
		# eg [0,1,0,1,0,0,1] -> inserts=[1,3,6], should be [0,1,3] (1-1,3-2,6-3)
		offset =  np.array([i+1 for i in range(len(self.insert_positions))],dtype='int64')+1
		self.insert_positions -= offset
		self.insert_positions += self.reference_start
		# Remove the inserts from the cigar list, since they distort subsequent positions
		self.cigarlist = self.cigarlist[self.cigarlist!=INSERT_INDEX]

		# muts - list of tuples of (mutation_base, mutations_position)
		# dels - list of tuples of (deletion_base, deletion_position)
		self.muts, self.dels = self.get_positions()
		self.number_mutations = len(self.muts)
		self.number_deletions = len(self.dels)

		# Generate a list of mutations, deletions relative to reference coordinates
		self.mutation_positions = np.array([e[1] for e in self.muts],dtype='int64') + self.reference_start
		self.deletion_positions = np.array([e[1] for e in self.dels],dtype='int64') + self.reference_start

		# indel_positions - list of positions that are both insertions and deletions (that is, a deletion followed by an insertion/s)
		# inmut_positions - list of positions that are both insertions and mutations (that is, a mutation followed by an insertion/s)
		self.indel_positions = np.intersect1d(self.deletion_positions, self.insert_positions)
		self.inmut_positions = np.intersect1d(self.mutation_positions, self.insert_positions)

		# Generate a list of (position, reference base, mutated base)
		self.aligned_sequence = self.generate_aligned_reference_sequence()

		# mutated_bases - list of tuples of (position relative to reference, original base, inserted base)
		self.mutated_bases = list(zip(self.mutation_positions, [e[0] for e in self.muts], self.aligned_sequence[[e[1]for e in self.muts]]))

		# Generate a list of (position, deleted base)
		self.deleted_bases = [(e[1],e[0]) for e in self.dels]

		# We want to discard all bases with low quality, but require a much higher standard for accepting a mutation than a WT base.
		# This is because our prior that a base is mutated is much lower, so we need more evidence to believe a mutation

		# Generate a list of mutated positions with low quality (<99%)
		below_99 = np.where(self.combine_quals() < 20)[0] + self.reference_start
		self.poor_quality_muts = np.intersect1d(self.mutation_positions,below_99)

		# Generate a list of WT bases with low quality
		self.poor_quality_WT = np.where(self.combine_quals() < 7)[0] + self.reference_start


		return None 

	def generate_aligned_reference_sequence(self):
		# Reference sequence should not include inserted bases 
		insert_pos = -1 # so that positions are zero indexed
		insert_index = []

		initial_softclip = self.cigartuples[0]
		for cigar_op in self.cigartuples:
			(op, length) = cigar_op
			if op != 1 and op !=2: #ignore insertions and deletions
				insert_pos += length
			if op == 1:
				for j in range(length):
					insert_index.append(insert_pos+j+1)
			
		remove_inserted = np.array([True]*len(self.query_sequence))
		remove_inserted[insert_index] = False
		if initial_softclip[0] == 4: # The first set of bases were softclipped
			remove_inserted[np.array(range(initial_softclip[1]))] = False #Throw away the softclipped bases from the reference string 
		reference_sequence = np.array(list(self.query_sequence))[remove_inserted]

		# add deletions to sequence
		for d in self.dels:
			(base, pos) = d
			reference_sequence = np.insert(reference_sequence, pos, 'D')

		return reference_sequence

	def parse_md(self):
		'''
		Parse MD string, separating text from numbers 
		12T3^CC0C -> [12,3,0] & ['T','^CC','C']

		@return: text - list of text portions of MD string
		@return: nums - list of integer portions of MD string
		'''
		a = list(self.MD)
		for i in range(len(a)):
			try: a[i]=int(a[i])
			except ValueError: pass

		text= []
		nums = []
		if isinstance(a[0], int):
			text_mode = False
		else:
			text_mode = True
			nums = [0]

		curr = a[0]
		for i in range(1, len(a)):
			if text_mode:
				if type(a[i-1])==type(a[i]):
					curr += a[i]
				else:
					text.append(curr)
					curr=a[i]
					text_mode = not text_mode
			else:
				if type(a[i-1])==type(a[i]):
					curr = curr*10+a[i]
				else:
					nums.append(curr)
					curr=a[i]
					text_mode = not text_mode
		if isinstance(curr, int):
			nums.append(curr)
		else:
			text.append(curr)
		return text, nums

	def get_positions(self):
		'''
		Finds the mutations, deletions and their respective positions, 
		zero-indexed from the start of the read

		@return: muts - list of tuples of (mutation_base, mutations_position)
		@return: dels - list of tuples of (deletion_base, deletion_position)
		'''
		text, nums = self.parse_md()
		muts = []
		dels = []
		pos = 0 
		for i in range(len(text)):
			pos += nums[i]
			if text[i][0]!='^':
				for j in range(len(text[i])):
					muts.append((text[i][j], pos+j))
					pos += 1
			else:
				dl = text[i]
				for d in range(1, len(dl)):
					dels.append((dl[d],pos))
					pos += 1
		return muts, dels

	def combine_quals(self):
		'''
		Computes the quality scores of each sequenced position
		@ params:   
			read - aligned read segment 
			dels - list of tuples (deletion_base, deletion_position)
			use_deletions - Bool:  True if deletions are considered as mutations
		@ output:
			v_pos - list of quality score relative to the start of the read
		'''
		d_pos = set([i for (j,i) in self.dels]) # list of positions containing deletions, relative to read start
		del_score = 30
		v_pos = [] # valid positions by quality score
		index = 0
		for i in range(self.reference_end - self.reference_start):
			if i in d_pos: 
				v_pos.append(del_score)
				continue
			v_pos.append(self.query_alignment_qualities[index])
			index += 1
		return np.array(v_pos)

	def get_original_query_sequence(self):
		return self.original_query_sequence

	def get_aligned_sequence(self):
		return self.aligned_sequence