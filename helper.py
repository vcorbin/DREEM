#--------------------- Helper functions for parsing MD string and converting positions to quality ------------------------------------
import pysam
import numpy as np
import pandas as pd
import multiprocessing
import time 

try:
	from .alignedRead import Read
except:
	from alignedRead import Read
	
import plotly 
import plotly.plotly as py
import plotly.graph_objs as go

def create_bitvectors(read, ref_trim, ref_lengths):    
	'''
	Create a mutations bitvector and deletions bitvector given a Read object and reference region of interest

	@param: read - 
	@param: ref_trim - 
	@param: ref_lengths - 

	@return: mut_vector
	@return: del_vector
	'''

	# Initialize empty vector
	mut_vector = np.array([np.nan]*ref_lengths[read.reference_name])
	del_vector = np.array([np.nan]*ref_lengths[read.reference_name])
	base_vector = np.array(['.']*ref_lengths[read.reference_name])

	# mark off all positions covered by the read
	pos = range(read.reference_start, read.reference_end)
	mut_vector[pos] = 0
	del_vector[pos] = 0
	base_vector[pos] = '0'
	# Identify positions where mutations/deletions occur

	if read.number_mutations > 0:
		mut_vector[read.mutation_positions] = 1
	if read.number_deletions > 0:
		del_vector[read.deletion_positions] = 1

	for (pos, original, mutated) in read.mutated_bases:
		base_vector[pos] = mutated

	
	for (base, pos) in read.dels:
		base_vector[pos+read.reference_start] = 'D'

	indel_vector = np.copy(base_vector)
	# Add information about insertions to bitvector
	indel_vector[read.insert_positions] = 'I'
	indel_vector[read.indel_positions] = 'X'
	indel_vector[read.inmut_positions] = 'M'

	if len(read.poor_quality_muts) > 0:
		# Remove poor quality mutations information
		mut_vector[read.poor_quality_muts] = np.nan
		del_vector[read.poor_quality_muts] = np.nan
		base_vector[read.poor_quality_muts] = '?'
		indel_vector[read.poor_quality_muts] = '?'

	if len(read.poor_quality_WT) > 0:
		# Remove poor quality WT base information
		mut_vector[read.poor_quality_WT] = np.nan
		del_vector[read.poor_quality_WT] = np.nan
		base_vector[read.poor_quality_WT] = '?'
		indel_vector[read.poor_quality_WT] = '?'

	# Trim vector to region of interest
	mut_vector = mut_vector[ref_trim[read.reference_name][0]:ref_trim[read.reference_name][1]]
	del_vector = del_vector[ref_trim[read.reference_name][0]:ref_trim[read.reference_name][1]]
	base_vector = base_vector[ref_trim[read.reference_name][0]:ref_trim[read.reference_name][1]]
	indel_vector = indel_vector[ref_trim[read.reference_name][0]:ref_trim[read.reference_name][1]]

	return mut_vector, del_vector, base_vector, indel_vector
	
def get_index(pre_index, index_length, read_sequence):
	'''
	Returns the unique index in a read, or None if no index is found
	@param: pre_index - A string which is the sequence immediately preceding the unique index
	@param: index_length - An integer which is the length of the unique index
	@param: read_sequence - A list which is the sequence of the read to search
	'''
	string_sequence = ''.join(read_sequence)
	if pre_index not in string_sequence:
		return ''
	index_start = string_sequence.index(pre_index)+len(pre_index)
	return read_sequence[index_start:min(index_start+index_length,len(read_sequence)-1)] + '?'*(index_start+index_length-len(read_sequence)-1)


def plot_simple_coverage(coverage, web_version=False):
	'''
	@param: coverage - numpy array showing the total coverage at each position.
	'''

	data = [go.Histogram(x=coverage)]

	if web_version != True:
		plotly.offline.plot(data , auto_open=False, filename=web_version)
	return plotly.offline.plot(data , auto_open=False, output_type='div', show_link=False)

def plot_miseq_read_coverage(coverage, web_version=False):

	rects = []
	for key in coverage:
		(r1_start, r2_start, r1_stop, r2_stop) = key
		count = coverage[key]
		start = max(r1_stop, r2_stop) 
		stop = min(r1_start, r2_start)
		width = stop - start
		uncovered_start = min(r1_stop, r2_stop) 
		uncovered_stop = max(r1_start, r2_start)
		uncovered_width = uncovered_stop - uncovered_start
		rects.append((count, start, width, uncovered_start, uncovered_width))

	rects.sort(key=lambda x: x[2])

	current_y = 0

	shapes = [] # keeps track of dictionaries of Rectangles

	min_x = np.inf
	max_x = -np.inf

	x = []
	y = []
	text = []

	for (count, start, width, uncovered_start, uncovered_width) in rects:
		min_x = min(min_x, start)
		max_x = max(max_x, start+width+uncovered_width)
		# create a rectangle object
		rect = {
				'type': 'rect',
				'x0': start,
				'x1': start+width,
				'y0': current_y,
				'y1': current_y + count,
				'line': {
					'color': 'rgba(50, 171, 96, 0.6)',
				},
				'fillcolor': 'rgba(50, 171, 96, 0.6)',
			}

		# add this to the shapes
		shapes.append(rect)

		# add the uncovered region if it exists
		if uncovered_width > 0:
			# create a rectangle object
			rect = {
					'type': 'rect',
					'x0': uncovered_start,
					'x1': uncovered_start+uncovered_width,
					'y0': current_y,
					'y1': current_y + count,
					'line': {
						'color': 'rgba(220,60,60,1) ',
					},
					'fillcolor': 'rgba(220,60,60,1) ',
				}

			# add this to the shapes
			shapes.append(rect)

		# add the annotations to the x_data and y_data for the trace
		x.append(0.5*(start+stop))
		y.append(current_y + 0.5*count)
		text.append('%d - %d : %d reads' %(start, stop, count))

		current_y += count

	trace0 = go.Scatter(
		x=x,
		y=y,
		text=text,
		mode='markers',
		opacity = 0,
	)
	data = [trace0]
	layout = {
		'xaxis': {
			'range': [0,1000],
			'showgrid': False,
		},
		'yaxis': {
			'range': [0, current_y+10]
		},
		'shapes': shapes,
	}
 
	fig = {
		'data': data,
		'layout': layout,
	}

	if web_version != True:
		print('Plotting read coverage for paired reads')
		plotly.offline.plot(fig, auto_open=False, filename=web_version)
		
#		Trying to get plots to save as static plots as well.
		# import plotly.offline as offline
		# from selenium import webdriver
		# driver = webdriver.PhantomJS(executable_path="phantomjs.exe")
		# driver.set_window_size(1000, 500)
		# driver.get(web_version)
		# driver.save_screenshot(web_version[:-5]+'.svg')

	return plotly.offline.plot(fig, auto_open=False,output_type='div', show_link=False)

		
def plot_hiseq_read_coverage(coverage, web_version=False):
	'''
	@param: coverage is a dictionary from (read_start, read_stop) : count
	'''

	current_y = 0

	shapes = [] # keeps track of dictionaries of Rectangles

	min_x = np.inf
	max_x = -np.inf

	x = []
	y = []
	text = []

	rects = []
	#collect rects in order of width
	for key in coverage:
		(start, stop) = key
		count = coverage[key]
		rects.append((start, stop, count, stop-start))
	rects.sort(key=lambda x: (x[0], x[3]))

	for elt in rects:
		(start, stop, count, width) = elt
		min_x = min(min_x, start)
		max_x = max(max_x, stop)
		# create a rectangle object
		rect = {
				'type': 'rect',
				'x0': start,
				'x1': stop,
				'y0': current_y,
				'y1': current_y + count,
				'line': {
					'color': 'rgba(50, 171, 96, 1)',
				},
				'fillcolor': 'rgba(50, 171, 96, 0.6)',
			}

		# add this to the shapes
		shapes.append(rect)

		# add the annotations to the x_data and y_data for the trace
		x.append(0.5*(start+stop))
		y.append(current_y + 0.5*count)
		text.append('%d - %d : %d reads' %(start, stop, count))

		current_y += count

	import plotly 
	import plotly.plotly as py
	import plotly.graph_objs as go

	trace0 = go.Scatter(
		x=x,
		y=y,
		text=text,
		mode='markers',
		opacity = 0,
	)
	data = [trace0]
	layout = {
		'xaxis': {
			'range': [min_x-10, max_x+10],
			'showgrid': False,
		},
		'yaxis': {
			'range': [0, current_y+10]
		},
		'shapes': shapes,
	}

	fig = {
		'data': data,
		'layout': layout,
	}

	if web_version==True:
		return plotly.offline.plot(fig, auto_open=False,output_type='div', show_link=False)
	else:
		print('Plotting read coverage for unpaired reads')
		plotly.offline.plot(fig, auto_open=False, filename=web_version)
		return None 

def plot_DMS_changes(all_mutated_bases, ref_start, ref_stop, web_version=True):
	'''
	@param: DMS_dictionary is a dictionary containing {position: [DMS_mutations]}
	'''

	from collections import defaultdict
	DMS_dictionary = defaultdict(list)
	sequence = {}
	for i in range(ref_start,ref_stop):
		sequence[i] = ['N']
	for (key, original, value) in all_mutated_bases:
		DMS_dictionary[key].append(value)
		if key in sequence:
			sequence[key].append(original)

	for key in sequence:
		sequence[key] = max(set(sequence[key]), key=sequence[key].count)

	# Consolidate sequence into list
	sequence_tup = sequence.items()
	sequence_tup = sorted(sequence_tup, key=lambda x: x[0])
	sequence = [s[1] for s in sequence_tup] # get letters only

	DMS_changes = {}
	DMS_changes['Position'] = list(range(ref_start+1, ref_stop+1))
	DMS_changes['A'] = [DMS_dictionary.get(pos,[]).count('A') for pos in range(ref_start,ref_stop)]
	DMS_changes['T'] = [DMS_dictionary.get(pos,[]).count('T') + DMS_dictionary.get(pos,[]).count('U') for pos in range(ref_start,ref_stop)]
	DMS_changes['C'] = [DMS_dictionary.get(pos,[]).count('C') for pos in range(ref_start,ref_stop)]
	DMS_changes['G'] = [DMS_dictionary.get(pos,[]).count('G') for pos in range(ref_start,ref_stop)]

	color = {'G':'goldenrod', 'U':'mediumturquoise', 'T':'mediumturquoise', 'A':'crimson', 'C':'royalblue', 'N':'grey'}
	
	trace1 = go.Bar(
		x=DMS_changes['Position'],
		y=DMS_changes['A'],
		name='A',
		marker=dict(
		color=color['A'],
		)

	)
	trace2 = go.Bar(
		x=DMS_changes['Position'],
		y=DMS_changes['T'],
		name='T',
		marker=dict(
		color=color['T'],
		)
	)

	trace3 = go.Bar(
		x=DMS_changes['Position'],
		y=DMS_changes['C'],
		name='C',
		marker=dict(
		color=color['C'],
		)
	)
	trace4 = go.Bar(
		x=DMS_changes['Position'],
		y=DMS_changes['G'],
		name='G',
		marker=dict(
		color=color['G'],
		)
	)

	data = [trace1, trace2, trace3, trace4]

	layout = go.Layout(
		barmode='stack',
		title='RT Inserted Bases',
		xaxis=dict(
			title='Position',
			),
		yaxis=dict(
			title='Number of Mutations',
			),
		)

	fig = go.Figure(data=data, layout=layout)
	if web_version != True:
		plotly.offline.plot(fig, auto_open=False, filename=web_version)	
		print('Plotted and saved DMS changes')
	return sequence, plotly.offline.plot(fig, auto_open=False,output_type='div', show_link=False)

def plot_mutation_distribution(basevectors, web_version, include_dels=False, include_coverage=False):
	'''
	Plot the distribution of mutations/coverage over all reads. 
	@param: basevectors - pandas df of basevectors
	@param: include_dels - True if deletions should be counted as mutations	
	@param: include_coverage - True if you wish to normalize by read coverage.
	'''
	import plotly
	import plotly.graph_objs as go

	if include_dels:
		mutations = basevectors['N_mutations'] + basevectors['N_deletions']
	else:
		mutations = basevectors['N_mutations']
	distribution = max(basevectors['Coverage'])*mutations /basevectors['Coverage']

	if include_coverage:	
		distribution = max(basevectors['Coverage'])*mutations/basevectors['Coverage']
	else:
		distribution = mutations

	data = [go.Histogram(x=distribution)]
	layout = go.Layout(
		title='Distribution of Mutations Per Read',
		xaxis=dict(
			title='Number of Mutations',
			),
		yaxis=dict(
			title='Number of Reads',
			),
		)
	fig = go.Figure(data=data, layout=layout)

	if web_version != True:
		print('Plotted and saved mutation distribution')
		plotly.offline.plot(fig , auto_open=False, filename=web_version)	
	return plotly.offline.plot(fig , auto_open=False,output_type='div', show_link=False)

def plot_population_average(bitvectors, columns, ref_start, ref_stop, web_version, sequence=None):
	# Compare replicates and generate a plot coloured by bases 
	colors = {'G':'goldenrod', 'U':'mediumturquoise', 'T':'mediumturquoise', 'A':'crimson', 'C':'royalblue', 'N':'grey'}

	import plotly
	import plotly.plotly as py
	import plotly.graph_objs as go
	from plotly import tools
	import numpy as np

	if sequence is None:
		sequence = ['N' for i in range(ref_start, ref_stop)]
	sequence = list(sequence)
	profile_fig = tools.make_subplots(rows=len(columns), cols=1, shared_yaxes=True, print_grid=False, subplot_titles=tuple(columns))
	col_num = 0
	for col in columns:
		bitreads = np.array(bitvectors[col].tolist())
		y = np.nansum(bitreads, axis=0)/bitreads.shape[0]
		col_num += 1
		trace = go.Bar(
				x=list(range(ref_start+1, ref_stop+1)),
				y= np.nansum(bitreads, axis=0)/bitreads.shape[0],
				text = [nucleotide for nucleotide in sequence],
				name = col,
				marker=dict(
						color=list(map(colors.get,sequence)),
				),   
		)

		profile_fig.append_trace(trace,col_num,1)
	layout = go.Layout(
	title='Reactivity Profile for Population Average',
	xaxis=dict(
		title='Position',
		),
	yaxis=dict(
		title='Mutation Probability',
		),
	)
	# profile_fig['layout'] = layout
	if web_version != True:
		print('Plotted and saved population average')
		plotly.offline.plot(profile_fig , auto_open=False, filename=web_version)
	return plotly.offline.plot(profile_fig , auto_open=False,output_type='div', show_link=False)

def plot_muts_per_base(mutvectors, sequence, web_version):
	'''
	Plot the number of mutations occurring on each nucleotide.
	@param: mutvectors is a column of a pandas array containing mutation vectors
	@param: sequence is a list containing he reference sequence
	@param: web_version is True if do not save any plots, otherwise the path to which plots should be saved
	'''
	if sequence is None:
		sequence = ['N' for i in range(ref_start, ref_stop)]
	sequence = list(sequence)


	bitreads = np.array(mutvectors.tolist())
	muts = np.nansum(bitreads, axis=0) # Counts the number of mutations at each base

	A = [i for i, letter in enumerate(sequence) if letter == 'A']
	C = [i for i, letter in enumerate(sequence) if letter == 'C']
	T = [i for i, letter in enumerate(sequence) if letter == 'T']
	G = [i for i, letter in enumerate(sequence) if letter == 'G']
	N = [i for i, letter in enumerate(sequence) if letter == 'N']

	A_count = np.nansum(muts[A])
	C_count = np.nansum(muts[C])
	T_count = np.nansum(muts[T])
	G_count = np.nansum(muts[G])
	N_count = np.nansum(muts[N])
	
	total = np.sum([A_count,C_count,T_count,G_count,N_count])
	AC_percent = (A_count+C_count)*100.0/total

	# Compare replicates and generate a plot coloured by bases 
	colors = {'G':'goldenrod', 'U':'mediumturquoise', 'T':'mediumturquoise', 'A':'crimson', 'C':'royalblue', 'N':'grey'}

	import plotly
	import plotly.plotly as py
	import plotly.graph_objs as go
	from plotly import tools
	
	trace = go.Pie(labels=['A','T','C','G','N'], 
					values=[A_count,T_count,C_count,G_count,N_count],
					 marker=dict(colors=[colors['A'],colors['T'],colors['C'],colors['G'],colors['N']])
				)

	layout = go.Layout(
		title='Percentage of DMS Modifications Per Base',
		)

	fig = go.Figure(data=[trace], layout=layout)

	if web_version != True:
		print('Plotted and saved mutations per base')
		plotly.offline.plot(fig , auto_open=False, filename=web_version)
	return plotly.offline.plot(fig , auto_open=False,output_type='div', show_link=False), AC_percent

def plot_indels(basevector, start, stop, sequence, web_version):
	'''
	Generate plot showing the signal occurring on mutated positions, deleted positions, insert positions, deletions followed by insertions and mutations following by insertions
	@param: basevector - DataFrame containing basevectors
	@param: start - start of region of interest
	@param: stop - stop of region of interest
	@param: sequence -  reference sequence corresponding to region of interest
	'''
	# Map mutations to 1, deletions to 2, insertions to 3, del-ins by 4, mut-ins by 5
	TRANSLATE_MAP =  {'.':np.nan, '?':np.nan, 'A':1.0, 'T':1.0, 'C':1.0, 'G':1.0, '0': 0.0, 'D':2.0, 'I':3.0, 'X':4.0, 'M':5.0}
	translated = np.array([list(map(TRANSLATE_MAP.get,list(bv))) for bv in basevector['Indels_vector']], dtype='float64')
	muts, dels, ins, delins, mutins = np.where(translated==1), np.where(translated==2), np.where(translated==3), np.where(translated==4), np.where(translated==5)

	colors = {'G':'goldenrod', 'U':'mediumturquoise', 'T':'mediumturquoise', 'A':'crimson', 'C':'royalblue', 'N':'grey'}

	all_bv = {}
	blank = np.zeros(translated.shape)
	blank[muts] = 1
	blank = blank.tolist()
	all_bv['Mutations_Only'] = blank #[''.join(x) for x in blank]

	blank = np.zeros(translated.shape)
	blank[dels] = 1
	blank = blank.tolist()
	all_bv['Deletions_Only'] = blank

	blank = np.zeros(translated.shape)
	blank[ins] = 1
	blank = blank.tolist()
	all_bv['Insertions_Only'] = blank

	blank = np.zeros(translated.shape)
	blank[delins] = 1
	blank = blank.tolist()
	all_bv['Deletions_Insertions'] = blank

	blank = np.zeros(translated.shape)
	blank[mutins] = 1
	blank = blank.tolist()
	all_bv['Mutations_Insertions'] = blank


	all_bv = pd.DataFrame(all_bv)


	indels_html = plot_population_average(all_bv, ['Mutations_Only','Deletions_Only','Insertions_Only','Deletions_Insertions','Mutations_Insertions'], start, stop, web_version, sequence)

# basevector = pd.read_csv('small_sample.bv', sep='\t', comment='@')
# start= 0
# stop = 100
# web_version='PLOT.html'
# sequence='NNNNGNNGNTGCACCTGATGNAATGGGAAAATATATCAAATCNTTCGTTGAGCGAGTTCTCNNAANTGANCANATGTCGACGGATCCTTTTTTAGGGAAG'
# plot_indels(basevector, start, stop, sequence, web_version)

def plots_from_basevector(basevector, start, stop, sequence, web_version):
	'''
	Generate plots from basevector, not BAM file.
	@param: basevector - DataFrame containing basevectors
	@param: start - start of region of interest
	@param: stop - stop of region of interest
	@param: sequence -  reference sequence corresponding to region of interest
	'''
	COMBINED_MAP =  {'.':np.nan, '?':np.nan, 'A':1.0, 'T':1.0, 'C':1.0, 'G':1.0, 'D':1.0, '0': 0.0}
	DELETIONS_MAP = {'.':np.nan, '?':np.nan, 'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0, 'D':1.0, '0': 0.0}
	MUTATIONS_MAP = {'.':np.nan, '?':np.nan, 'A':1.0, 'T':1.0, 'C':1.0, 'G':1.0, 'D':0.0, '0': 0.0}

	# Convert basevectors to get combined_vectors, mutation_vector and deletion_vector
	basevector['Combined_vector'] = [list(map(COMBINED_MAP.get,list(bv))) for bv in basevector['Bases_vector']]
	basevector['Mutations_vector'] = [list(map(MUTATIONS_MAP.get,bv)) for bv in basevector['Bases_vector']]
	basevector['Deletions_vector'] = [list(map(DELETIONS_MAP.get,bv)) for bv in basevector['Bases_vector']]

	# Plot population average
	pop_avg_html = plot_population_average(basevector, ['Combined_vector', 'Mutations_vector', 'Deletions_vector'], start, stop, web_version='web_version.html', sequence=sequence)
	# Plot muts per base
	muts_per_base_html = plot_muts_per_base(basevector['Mutations_vector'], sequence, web_version=web_version)
	#Plot muts per read
	muts_per_read_html = plot_mutation_distribution(basevector, web_version=web_version, include_dels=False, include_coverage=False)

	return pop_avg_html, muts_per_base_html, muts_per_read_html
	