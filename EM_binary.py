import numpy
import scipy.misc, scipy.stats
import multiprocessing
try:
    # website imports. IMPT: Include all directory imports in the try/except pattern to not break website code.
    from . import EM_files as files
    from . import EM_plots as plots
    from . import denominator as denom
    from .EM_lib import calculate_cluster_probs
except:
    import EM_files as files
    import EM_plots as plots
    import denominator as denom
    from EM_lib import calculate_cluster_probs

from scipy.optimize import newton_krylov
import time


# TODO: cluster_probs that is recorded is not the normalised. Doesn't affect the clustering but should be changed

def runEM(x, K, labels, threading=True, init_clusters=None, init_cluster_probs=None, init_plot=False, max_its=1000, seed=None,
          workDir=None, dmuts=3):
    """
    Runs the EM algorithm on the x dataset, separating it into d clusters. Creates a subdirectory containing:
    - cluster probabilities distribution text file
    - static plot of cluster probabilities
    - cluster responsibilities text file
    - loglik text file 

    @param: x - class containing 3 parameter: 'bitMatrix': n x l pandas dataFrame, where n is the number of independent reads, and l is the length/dimension of a bitvector, 'N': numpy array of size n, it shows the number of times this read occurs; 'bitStrings': numpy array of size n.
    @param: K - number of clusters into which the data will be separated
    @param: labels - Dictionary containing 2 Ntot x 1 pandas dataFrame, where Ntot is the number of reads: 'refs' and 'bitstring'. labels['refs'] contains the reference associated with each read in the bit-vector file, labels['bitStrings'] contains the bit- strings associated with each read in the bit-vector file.
    @param: init_clusters - prior distribution over clusters
    @param: init_clusters_probs - prior distribution over cluster relative abundance
    @param: max_its - maximum number of iterations for which a single run of EM is allowed to run.
    @param: a - 
    @param: b - 
    @param: workDir - 
    @param: seed - 
    @param: mask_threshold - mask positions with mutations above this decimal threshold; probably endogeneous modifications. default 0.25
    @param: ref_coords - dictionary from reference name: [start, stop] for relative plotting
    @param: dmuts - minimum distance between 2 mutations on a same read

    @return: loglik - the log likelihoods of this clustering run
    @return: clusters - numpy array of the final cluster probabilities 
    """
    start = time.time()
    print(x.bitMatrix[0])
    print(x.N[0])
    print(x.bitStrings[0])
    ref_coords = {}
    for ref in x.refs_info.refs.keys():
        ref_coords[ref] = [x.refs_info.refs[ref]['start'], x.refs_info.refs[ref]['start']]

    N = x.bitMatrix.shape[0]
    if not threading:
        print("Not using threading.")
        nthreads = 1
    else:
        nthreads = min(max(1, N / 40000), 18)
        nthreads = min(nthreads / K, multiprocessing.cpu_count() / K)
        nthreads = int(max(1, nthreads))
        nthreads = 18
        print("Using threading: %d threads." % nthreads)

    print("Using %d thread(s)." % nthreads)

    timestamps = []
    lls = []
    cps = []
    I = x.bitMatrix.shape[1]

    if init_clusters is None:
        if seed is not None:
            numpy.random.seed(seed)
        else:
            pass
        mu = numpy.asarray([scipy.stats.beta.rvs(1.5, 20, size=I) for k in range(K)]) # initize clusters with random values
    else:
        # Transform into floats in case the init_cluster is made of string
        mu = []
        for MPD in init_clusters:
            mu.append(numpy.asarray(map(float, MPD)))
        mu = numpy.asarray(mu)
    if init_cluster_probs is None:
        cluster_probs = numpy.asarray([1.0 / K] * K)  # intitialize uniform distribution cluster probs
    else:
        # Transform into floats in case the init_cluster is made of string
        cluster_probs = map(float, init_cluster_probs)

    timestamps.append(time.time())
    print("Starting Expectation step #0.")
    responsibilities, last_loglik = Expectation_step(x, mu, numpy.log(cluster_probs), nthreads=nthreads)

    lls.append(last_loglik)
    timestamps.append(time.time())

    if init_plot and type(workDir) == str:
        plots.PlotClusters(x, mu, cluster_probs=cluster_probs, outFilename=workDir + "/init_cluster_probs.png",
                           genes=ref_coords.keys(), start=ref_coords.values()[-1][0], run_logL=last_loglik)
        print("Initial likelihood =", last_loglik)

    # ----------------------------------------------- START OF ITERATIONS ----------------------------------------
    for its in range(max_its):
        print("Iteration: %d" % its)
        ## re-estimate cluster parameters
        #print("Starting Maximisation step #%d." % its)
        mu = Maximisation_step(x, mu, responsibilities, dmuts=dmuts)
        timestamps.append(time.time())

        real_cluster_probs = calculate_cluster_probs(responsibilities, x, mu, dmuts=dmuts, normalise=False)
        real_log_prob = numpy.log(real_cluster_probs)
        real_cluster_probs = real_cluster_probs / sum(real_cluster_probs)
        cps.append(real_cluster_probs)

        # Find the responsibilities for each read
        #print("Starting Expectation step #%d." % (its + 1))
        responsibilities, loglik = Expectation_step(x, mu, real_log_prob, nthreads=nthreads)
        timestamps.append(time.time())
        lls.append(loglik)

        # if the likelihood has converged (locally), stop iterating
        if lls[-1] < last_loglik:
            print("ERROR: The log-likelihood is decreasing: %f -> %f" % (last_loglik, lls[-1]))
            #print("Mu is %s" % mu)
            #if lls[-1] - last_loglik < -0.001 and last_loglik - lls[-3] < -0.001 or (K==1 and lls[-1] - last_loglik < 0.001):
            #if K == 1 and lls[-1] - last_loglik < 0.001:
                #print("Breaking on iteration: %d" % its)
                #break
        # if (numpy.abs(lls[-1] - last_loglik) / numpy.abs(lls[-1]) < 0.000) and its >= 100:
        #     break
        last_loglik = lls[-1]

    # ---------------------------------------------- END OF ITERATIONS -------------------------------------------
    end = time.time()
    print("Execution time: %f" % (end - start))
    cluster_probs = calculate_cluster_probs(responsibilities, x, mu, dmuts=dmuts)  # final (real) cluster probs
    print("Final real cluster probs: %s" % ", ".join([str(x) for x in cluster_probs]))
    if type(workDir) == str:
        files.SaveRun(x, lls, responsibilities, mu, labels["refs"], workDir)
        plots.PlotRun(x, last_loglik, lls, cps, responsibilities, mu, labels, workDir, ref_coords, timestamps, cluster_probs=cluster_probs)

    return lls, mu, cps


def Expectation_step(x, mu, real_log_prob, nthreads=1):
    N = x.bitMatrix.shape[0]
    I = x.bitMatrix.shape[1]
    K = len(real_log_prob)
    Px = numpy.zeros((N, I, K))  # Log of probability of a base in vector Xn given mus in cluster k
    bitMatrix = x.bitMatrix.values

    if nthreads == 1:
        for k in range(K):
            Px[:, :, k] = scipy.stats.bernoulli.logpmf(bitMatrix, mu[k])
    else:
        processes = []
        starts = []
        #print("Start multithreading E step with %d parallel processes" % (nthreads * K))
        Epool = multiprocessing.Pool(processes=nthreads)
        start = 0
        for i in range(nthreads):
            starts.append(start)
            end = start + int(N/nthreads)
            if i == nthreads - 1:
                end = N
            for k in range(K):
                processes.append(Epool.apply_async(logpmf_wrapper, args=(bitMatrix[start:end], mu[k])))
            start = end
        Epool.close()
        Epool.join()
        for i, process in enumerate(processes):
            k = i % K
            start = starts[int(i/K)]
            res = process.get()
            Px[start:start+len(res), :, k] = res
        #print("End of E step multithreading")

    Px = numpy.nansum(Px, axis=1)  # For each read, calculate P(x | mu_k) = prod_i [mu_ik^x_ni * (1 - mu_ik^(1 - x_ni))]
    responsibilities = numpy.add(Px, real_log_prob)  # responsibilities = log (P(x | mu_k)*pi_k)
    temp = scipy.misc.logsumexp(responsibilities, axis=1)  # temp = sum_k P(x| mu_k)*pi_k
    responsibilities = numpy.subtract(responsibilities.T, temp).T  # Like dividing, because in log
    responsibilities = numpy.exp(responsibilities)  # THIS IS THE FINAL VALUE
    loglik = numpy.dot(temp, x.N)
    return responsibilities, loglik


def logpmf_wrapper(x, mu):
    return scipy.stats.bernoulli.logpmf(x, mu)


def logmpf_wrapper(x, mu, start, k):
    return scipy.stats.bernoulli.logpmf(x, mu), start, k


def Maximisation_step(x, mu, responsibilities, dmuts=3):
    K = responsibilities.shape[1]
    for k in range(K):
        mu[k] = M_helper(x, responsibilities, k, mu[k], dmuts=dmuts)
    return mu


def M_helper(x, responsibilities, k, init_mu, dmuts=3):
    I = x.bitMatrix.shape[1]
    mup_k = numpy.nansum(responsibilities[:, k] * x.N * x.bitMatrix.values.T, axis=1) / \
            numpy.nansum(responsibilities[:, k] * x.N * numpy.isfinite(x.bitMatrix.values).T, axis=1)
    mu = newton_krylov(lambda mu_k: Der_fun(mu_k, mup_k, dmuts=dmuts),init_mu)
    zeros = numpy.where(mu <= 0.0)[0]
    if len(zeros) > 0:
        #print "MU has negative values: %s" % mu
        mu = numpy.maximum(mu, numpy.zeros(I))
    return mu


def Der_fun(mu_k, mup_k, dmuts=3):
    I = len(mu_k)
    mu_k_rev = mu_k[::-1]
    den, dens_k = denom.smart_denominator(mu_k, dmuts=dmuts)
    den, dens_k_rev = denom.smart_denominator(mu_k_rev)
    Der_ll = lambda alpha: denom.den_alpha(mu_k, mu_k_rev, alpha, dens_k, dens_k_rev, dmuts=dmuts) / den - mup_k[alpha]
    return numpy.asarray([Der_ll(alpha) for alpha in range(I)])


def print_everything(clusterprobs, loglik, mu, L_corr_clust):
    print( "MU is %s" % mu)
    print( "Cluster Probs is %s" % clusterprobs)
    print( "Loglik is %s " % loglik)
    print( L_corr_clust)

def fix(bitvector, dmuts=3):
    I = len(bitvector[0])
    for i in range(len(bitvector)):
        for j in numpy.where(bitvector[i] == 1)[0]:
            bitvector[i][max(0,j-dmuts):min(I,j+dmuts+1)] = [0]*(min(I,j+dmuts+1) - max(0,j-dmuts))
            bitvector[i][j] = 1
