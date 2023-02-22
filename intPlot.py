#/Users/vcorbin/libs/anaconda3/envs/py27/bin/python

import pandas as pd
from plotly.offline import plot

def intPlot(filename):
	x = pd.DataFrame.from_csv(filename,sep="\t")
	plot([{'x':x.index,'y':x[col],'name':col} for col in x.columns])

if __name__ == "__main__":

	import argparse

	parser = argparse.ArgumentParser(description='Takes a file of the type clusters.txt as an argument, and makes an offline interactive plot using plotly.offline')
	parser.add_argument('input_filename',type=str,
		help='clusters.txt file')
	args = parser.parse_args()

	intPlot(args.input_filename)