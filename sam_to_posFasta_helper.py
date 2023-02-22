#--------------------- Helper functions for parsing MD string and converting positions to quality ------------------------------------
import pysam
import numpy as np
import pandas as pd
import multiprocessing
import time 
from alignedRead import Read

import plotly 
import plotly.plotly as py
import plotly.graph_objs as go

def create_posFasta(read, ref):    
	'''
	Create a mutations bitvector and deletions bitvector given a Read object and reference region of interest

	@param: read - 
	@param: ref - 
	'''

	# Initialize empty vector
	mut_vector = np.array([np.nan]*ref_lengths[read.reference_name])
	del_vector = np.array([np.nan]*ref_lengths[read.reference_name])
	base_vector = np.array(['.']*ref_lengths[read.reference_name])

	# mark off all positions covered by the read
	pos = range(read.reference_start, read.reference_end)
	# Identify positions where mutations/deletions occur

	if read.number_deletions > 0:
		del_vector[read.deletion_positions] = 1

	for (pos, original, mutated) in read.mutated_bases:
		base_vector[pos] = mutated

	# Trim vector to region of interest

	return mut_vector, del_vector, base_vector
	
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
	return plotly.offline.plot(data , auto_open=False, output_type='div')

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

	return plotly.offline.plot(fig, auto_open=False, output_type='div')

		
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
		return plotly.offline.plot(fig, auto_open=False, output_type='div')
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

	trace1 = go.Bar(
	    x=DMS_changes['Position'],
	    y=DMS_changes['A'],
	    name='A'
	)
	trace2 = go.Bar(
	    x=DMS_changes['Position'],
	    y=DMS_changes['T'],
	    name='T'
	)

	trace3 = go.Bar(
	    x=DMS_changes['Position'],
	    y=DMS_changes['C'],
	    name='C'
	)
	trace4 = go.Bar(
	    x=DMS_changes['Position'],
	    y=DMS_changes['G'],
	    name='G'
	)

	data = [trace1, trace2, trace3, trace4]
	layout = go.Layout(
	    barmode='stack'
	)

	fig = go.Figure(data=data, layout=layout)
	if web_version == True:
		return sequence, plotly.offline.plot(fig, auto_open=False, output_type='div')
	else:
		print('Plotted and saved DMS changes')
		plotly.offline.plot(fig, auto_open=False, filename=web_version)
		return sequence, None 

def plot_mutation_distribution(basevectors, web_version, include_dels=False):
	'''
	Plot the distribution of mutations/coverage over all reads. 
	@param: basevectors - pandas df of basevectors
	@param: include_dels - True if deletions should be counted as mutations
	'''
	import plotly
	import plotly.graph_objs as go

	if include_dels:
		mutations = basevectors['N_mutations'] + basevectors['N_deletions']
	else:
		mutations = basevectors['N_mutations']
	distribution = max(basevectors['Coverage'])*mutations/basevectors['Coverage']
	data = [go.Histogram(x=distribution)]

	if web_version == True:
		return plotly.offline.plot(data , auto_open=False, output_type='div')
	else:
		plotly.offline.plot(data , auto_open=False, filename=web_version)
		return None

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
	profile_fig = tools.make_subplots(rows=len(columns), cols=1, shared_yaxes=True, print_grid=False)
	col_num = 0
	for col in columns:
		bitreads = np.array(bitvectors[col].tolist())
		col_num += 1
		trace1 = go.Bar(
				x=list(range(ref_start+1, ref_stop+1)),
				y= np.nansum(bitreads, axis=0)/bitreads.shape[0],
				#text = [str(t) for t in zip(np.nansum(bitreads, axis=0), sequence)],
				name = col,
				marker=dict(
						color=list(map(colors.get,sequence)),
				),   
		)

		profile_fig.append_trace(trace1,col_num,1)
	
	# profile_fig['layout'].update(title=title)
	if web_version == True:
		return plotly.offline.plot(profile_fig , auto_open=False, output_type='div')
	else:
		plotly.offline.plot(profile_fig , auto_open=False, filename=web_version)
		return None
