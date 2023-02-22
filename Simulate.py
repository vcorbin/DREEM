import matplotlib
matplotlib.use('agg')
import numpy
import random


def Simulated_to_binary(clusters,clusters_prob,outFilename,nreads=1000000,read_length=50):

    if len(clusters) != len(clusters_prob): return("Error: the number of clusters doesn't match the number of cluster probabilities provided.")
    full_read_length = len(clusters[0])
    for clust in clusters:
        if len(clust) != full_read_length: return("Error: the elements of 'clusters' don't have the same length.")
    if sum(clusters_prob) != 1: return("Error: the cluster probabilitites don't add up to 1.")

    # Record the clusters information and probabilities
    outLogName = outFilename.split(".")
    outLogName = ".".join(outLogName[:(len(outLogName)-1)]) + "_clusters_info.txt"
    outL = open(outLogName,"w")
    outL.write("Clusters_prob:\n")
    probs = map(str,clusters_prob)
    outStr = " ".join(probs) + "\n\n"
    outL.write(outStr)
    outL.write("Clusters mutations prob:\n")
    for clust in clusters:
        clust = map('{0:.8f}'.format,clust)
        outStr = " ".join(clust) + '\n'
        outL.write(outStr)
    outL.close()

    outF = open(outFilename,"w")
    seq = "".join(["N"] * full_read_length)
    #outstr = "@start\tref\t0\n@length\tref\t%d\n@seq\tref\t%s\n" % (full_read_length, seq)
    outstr = "@ref\tsimulation;simulation\tseq\n@coordinates:length\t0,%d:%d\n" % (full_read_length - 1, full_read_length)
    outF.write(outstr)
    outF.write("Query_name\tBit_vector\tN_Mutations\n")

    for j in range(nreads):

        read = numpy.empty(full_read_length,dtype='str')
        read.fill('?')
        
        # Choose the cluster from which the read will be drawn 
        x = random.random()
        clust_prob = 0.0
        for i in range(len(clusters_prob)):
            clust_prob += clusters_prob[i]
            if x <= clust_prob:
                rand_cluster = clusters[i]
                break
            else:
                pass

        # Choose where the aligned read starts between 1-read_length and full_read_length-1
        #start_pos = random.randint(1-read_length,full_read_length - 1)
        #end_pos = min(full_read_length,start_pos + read_length)
        #start_pos = max(0,start_pos)
        start_pos = 0
        end_pos = start_pos + full_read_length

        # Fill the corresponding region with 0s:
        read[start_pos:end_pos] = '0'

        # Choose which locations in the read will contain mutations
        read[start_pos:end_pos][numpy.random.random(end_pos-start_pos)<=rand_cluster[start_pos:end_pos]] = '1'
        nmut = numpy.count_nonzero(read=='1')

        # Write the mutated read onto the binary file
        outStr = str(j+1) + "\t" + "".join(read) + "\t" + str(nmut) + "\n"
        outF.write(outStr)

    outF.close()


def LoadSimInfo(filename):
    """
    Load the cluster_probs and mus from the info file created by Simulated_to_binary()
    :param filename: Name of the file
    :return: pis and mus
    """
    def clust_parse(x):
        return numpy.asarray(list(map(float, x.strip().split(" "))))
    F = open(filename, "r")
    next(F)
    line = next(F)
    pi = clust_parse(line)
    next(F)
    next(F)
    mu = []
    for line in F:
        mu.append(clust_parse(line))
    mu = numpy.asarray(mu)
    return mu, pi