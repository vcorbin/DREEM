import pandas as pd 
import numpy as np
import sys

def stringify_vec(vec, whichvector='bit'):
        vector = np.array(list(vec))
        if whichvector == 'bit' or whichvector == 'mut':
            vector[vector=='A'] = '1'
            vector[vector=='T'] = '1'
            vector[vector=='G'] = '1'
            vector[vector=='C'] = '1'
            vector[vector=='.'] = '?'
        if whichvector == 'bit' or whichvector == 'del':
            vector[vector=='D'] = '1'
        return "".join(vector)

if __name__ == "__main__":
    filename = sys.argv[1]
    df = pd.read_csv(filename, sep="\t", header=3)#, header = ['Query_name','Bases_vector','N_mutations','N_deletions','Coverage','Reference','Index'])
    bases = df['Bases_vector']
    df['Bases_vector'] = df['Bases_vector'].apply(stringify_vec)
    df.rename(columns={'Bases_vector': 'Binary_Vector'}, inplace=True)
    df.to_csv(filename[:-4] + "_zero_ones" ".txt", sep="\t", index=False)