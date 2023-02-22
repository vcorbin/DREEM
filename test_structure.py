#!/usr/bin/python3
import sys
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon, ranksums, mannwhitneyu, ttest_1samp


class Prediction(object):
    def __init__(self, mu, dot, n_cluster = 1):
        self.mu = mu
        self.dot = dot
        self.n_cluster = n_cluster
        df_cluster = pd.read_csv(mu, sep='\t', skiprows=3)
        self.cluster = df_cluster['Cluster_'+str(n_cluster)]

    def pos(self):
        with open(self.mu, 'r') as f:
            first_line = f.readline()
            start_pos = int(first_line.split("\t")[2])
            second_line = f.readline()
            length = int(second_line.split("\t")[2])
            end_pos = start_pos + length - 1
        f.close()
        return start_pos, end_pos

    def dot_structures(self):
        start_pos, end_pos = self.pos()
        search_list = ['(',')']
        dot_structures = []
        with open(self.dot, "r") as f:
            for count,line in enumerate(f):
                if all([item in line for item in search_list]):
                    dot_structures.append(line[start_pos:end_pos+1])
                # if all([nuc in line for nuc in ['A','T']]):
                #     original = line[:-1]  # original RNA sequence
        f.close()
        return dot_structures

    def paired(self, s):
        paired = self.cluster.where([i == '(' or i == ')' for i in s]).dropna()
        return np.asarray(paired)

    def unpaired(self, s):
        unpaired = self.cluster.where([i == '.' for i in s]).dropna()
        return np.asarray(unpaired)

    def sort_pvals(self, test=mannwhitneyu):
        pvals = []
        for s in self.dot_structures():
            paired, unpaired = self.paired(s), self.unpaired(s)
            #     print(paired.median(),unpaired.median())
            pval = test(paired, unpaired, use_continuity=True, alternative='less').pvalue
            pvals.append(pval)

        print("Sort P-values for mannwhitneyu test: ")
        print(np.argsort(pvals) + 1)
        print("Best P-value: " + str(np.min(pvals)))
        return pvals


if __name__ == '__main__':
    mu = sys.argv[1]
    dot = sys.argv[2]
    n_cluster = sys.argv[3]

    pre = Prediction(mu, dot, n_cluster)
    pvals = pre.sort_pvals()
    pvals_argsort = np.argsort(pvals) + 1
    with open("vienna_prediction_cluster.txt", "a") as f:
        f.write("P values sorted for cluster " + str(n_cluster) + ':\n[ ' + str(pvals_argsort).strip('[]') + ' ]\n')
        f.write("The best P value: " + str(np.min(pvals)) + '\n')
    f.close()
