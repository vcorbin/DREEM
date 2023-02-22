import matplotlib
import matplotlib
matplotlib.use('Agg') # Remove the need to run a X server
import pylab
import os
import numpy
try:
    from . import EM_files as files
except:
    import EM_files as files
from matplotlib import pyplot as plt

def PlotRun(x,run_logL, run_lls, cps, run_resp,run_clusters,labels,outDir,ref_coords, times, cluster_probs=None):
    
    outDir = outDir + '/'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    #Can be generalised later to allow plotting of different co-ordinates
    genes = ref_coords.keys()
    start = list(ref_coords.values())[-1][0]
    start = list(ref_coords.values())[-1][1]

    PlotLLs(run_lls, outDir + "log_liks.png")
    PlotTimes(times, outDir + "timestamps.png")
    PlotCPs(numpy.array(cps), outDir + "cp_progression.png")

    PlotClusters(x, run_clusters, cluster_probs=cluster_probs, outFilename=outDir + "cluster_probs.png",
                 run_resp=run_resp, run_logL=run_logL, genes=genes, start=start)

    ## Violin plot
    #Ntot = len(labels["refs"])
    #full_resp = numpy.zeros((Ntot,run_resp.shape[1]))
    #for i,bitString in enumerate(labels["bitStrings"][0]):
    #    bitRow = numpy.where(x.bitStrings == bitString)[0][0]
    #    full_resp[i,:] = run_resp[bitRow,:]
#
    #combined = defaultdict(list)
    #prob_label = "predicted membership probability"
    #for j in range(len(run_clusters)):
    #    combined["allele"].extend(labels["refs"][0])
    #    combined[prob_label].extend(full_resp[:,j])
    #    combined["cluster"].extend([str(j+1)]*Ntot)
    #combined = pandas.DataFrame(combined)
    #try:
    #    seaborn.violinplot(y=prob_label, x="cluster", hue="allele", split=True, inner="quartile",data=combined)
    #except ValueError:
    #    seaborn.violinplot(y=prob_label, x="cluster", hue="allele", split=False, inner="quartile",data=combined)
    #pylab.title("Predicted cluster memberships (%d reads, loglik = %.2f)" % (Ntot,run_logL))
    #pylab.savefig(outDir + "/cluster_memberships.png")
    #pylab.close()

def PlotTimes(times, outFilename):
    E_times = [times[2*i + 1] - times[2*i] for i in range(len(times)//2 - 1)]
    M_times = [times[2*i + 2] - times[2*i + 1] for i in range(len(times)//2 - 2)]
    mean_e = numpy.mean(E_times)
    mean_m = numpy.mean(M_times)
    x_1 = [a for a in range(len(E_times))]
    x_2 = [a for a in range(len(M_times))]
    plt.plot(x_1, E_times, color='r')
    plt.plot(x_2, M_times, color='b')
    plt.title('Timestamps: Mean E step: %s, Mean M step: %s.' % (mean_e, mean_m))
    plt.savefig(outFilename[:])
    plt.clf()
    plt.cla()

def PlotLLs(run_lls, outFilename):
    x = [a for a in range(len(run_lls))]
    y = run_lls
    plt.plot(x[2:], y[2:])
    plt.savefig(outFilename[:])
    plt.clf()
    plt.cla()

def PlotCPs(cps, outFilename):
    cpt = cps.T
    ind = numpy.arange(len(cps))    # the x locations for the groups
    width = 0.8       # the width of the bars

    colors = ['r', 'b', 'y']
    for i in range(len(cpt)):
        plt.bar(ind, cpt[i], width, bottom=numpy.sum(cpt[:i], axis=0), color=colors[i % 3])

    plt.ylabel('Cluster Probs')
    plt.title('Cluster Probs Progression')
    plt.savefig(outFilename[:])
    plt.clf()
    plt.cla()

def PlotClusters(x, run_clusters, cluster_probs=None, xlim=None, outFilename=None, run_resp=None, run_logL=0.0,
                 genes=None, start=0):
    '''
        clusters can be either a pandas DataFrame, or the name of a file containing the clusters.
    '''
    if type(run_clusters) == str: run_clusters = files.LoadClusters(run_clusters)
    d = len(run_clusters)
    if cluster_probs is not None:
        pass
    elif run_resp is not None:
        cluster_probs = numpy.sum(run_resp.T*x.N, axis=1)/x.N.sum()
    else:
        cluster_probs = ['NA'] * d

    # Plot the mutational probabilities for each cluster:
    xvalues = [start + i for i in range(1,len(run_clusters[0])+1)]
    for i in range(d):
        pylab.plot(xvalues,run_clusters[i], 'o--', label="cluster %d ($p_%d=%s)"
                                                         % (i+1, i+1, str(cluster_probs[i])[0:5]))
    pylab.ylim(0.0,0.3)
    if type(xlim)==list and len(xlim)==2: pylab.xlim(xlim[0],xlim[1])
    #pylab.ylabel("Normalized mutation probability")
    pylab.ylabel("Mutation probability")
    pylab.xlabel("position (bp)")
    pylab.legend(loc="upper left")

    if type(genes) == list or type(genes) == str:
        if type(genes) == list and len(genes) == 1:
            genes = genes[0]
        elif type(genes) == list and len(genes) > 1:
            genes = ', '.join(genes[:-1]) + ' and ' + genes[-1]
        pylab.title("Clusters on {0} with loglik = {1:.2f}".format(genes,run_logL))
    else:
        pylab.title("loglik = %.2f" % run_logL)
    
    if type(outFilename) is str: 
        pylab.savefig(outFilename)
        pylab.close()
    else:
        pylab.show()

def PlotJackPotting(jackpot, plot_name):
    """
    
    :param jackpot: object of class jackPotting
    :param filename: string, file name without extension that will be used as prefix for the various plots
    :return: 
    """

    pylab.figure()
    pylab.scatter(range(len(jackpot.N_dist)), jackpot.N_dist)
    pylab.title('Deviation from the expected number of occurences')
    pylab.ylabel('Deviation (standard deviations)')
    pylab.xlabel('Reads (ordered by expected number of occurences)')
    if plot_name is not None:
        pylab.savefig(plot_name + '_sds.png')
        pylab.close()
    pylab.figure()
    pylab.scatter(numpy.log10(jackpot.N_theory), jackpot.N_dist)
    pylab.title('Deviation from the expected number of occurences vs expected number of occurences')
    pylab.ylabel('Deviation (standard deviations)')
    pylab.xlabel('log10 N_expected')
    if plot_name is not None:
        pylab.savefig(plot_name + '_sds_vs_N.png')
        pylab.close()
    pylab.figure()
    pylab.scatter(numpy.log10(jackpot.N_theory), numpy.log10(jackpot.Ns))
    pylab.title('N observed vs expected number of occurences')
    pylab.ylabel('log10 N_observed')
    pylab.xlabel('log10 N_expected')
    if plot_name is not None:
        pylab.savefig(plot_name + '_Nobs_vs_Nexp.png')
        pylab.close()














































































