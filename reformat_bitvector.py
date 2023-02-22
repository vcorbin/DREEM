import pandas as pd
import argparse

def main():
	parser = argparse.ArgumentParser(description='Convert aligned DMS-modified reads to bitvectors')

	parser.add_argument('-f','--input_file',type=str,dest='inFile', required=True,
	help='Bitvector file to change into new format')


	args = parser.parse_args()
	print(args)

	#--------------------------------------------- Variables conversion ----------------------------------
	inFile = args.inFile


	df = pd.read_csv(inFile, sep='\t', names=['Query_name', 'Binary_vector', 'N_mutations', 'Reference_name', 'Start_position'])
	del df['Start_position']
	df.to_csv(inFile, sep='\t', index=False)