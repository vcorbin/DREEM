import argparse
import numpy

def parse_args(convert_args=True):

    '''
    List of arguments and corresponding variables

    Mandatory:
                input_filename --> input_filename

    Optional:
                -d, --ncluster --> ncluster
                -n, --EM_its --> EM_its
                -i, --max_its --> max_its
                -r, --nreads --> nreads
                -m, --muts_threshold --> muts_threshold
    '''

    parser = argparse.ArgumentParser(description='Run a EM to cluster reads from a bam/sam file or a binary_vector file')

    parser.add_argument('input_filename',type=str,
        help='String. Preferably a binary_vector file (.txt) created by sam_to_binary(). If a sam/bam is given, a binary_vector file of the same name will be created (this part is not implemented yet.')

    parser.add_argument('-d','--ncluster',type=int,dest='ncluster',default=2,
        help='Integer, default: 2. Number of clusters used to run the EM.')

    parser.add_argument('-n','--EM_its',type=int,dest='EM_its',default=10,
        help='Integer, default: 10. Number of EM runs.')

    parser.add_argument('-i','--max_its',type=int,dest='max_its',default=300,
        help='Integer, default: 300. Number of iterations allowed in an EM run before it is stopped. The run will stop prior to this if convergence is obtained.')

    parser.add_argument('-r','--nreads',type=float,dest='nreads',default=1000000,
        help='Integer, default: 1 million. Number of reads to be loaded from the bit-vector file.')

    parser.add_argument('-m','--muts_threshold',type=int,dest='muts_threshold',default=0,
        help='Integer, default: 0. Minimum number of mutations for a read to be used in the clustering.')

    parser.add_argument('-M','--mask_threshold',type=float,dest='mask_threshold',default=25.0,
        help='Float, default: 25.0. Maximum threshold of mutation percent, over which positions are considered endogeneous and are masked')

    parser.add_argument('-I','--mask_index',type=str,dest='mask_index',default='no',
        help='String, default: "no". Index of the bases to be masked comma separated ("0,1,3").')

    parser.add_argument('-q','--q_threshold',type=float,dest='na_threshold',default=0.2,
        help='Float, default: 0.2. Threshold of the ratio of "?", "." and "N" a read from bit-vector file is allowed to have to be loaded.')

    parser.add_argument('-P', '--max_mut_prob', type=float, dest='max_mut_prob', default=-1.0,
        help='float, default: -1.0. Maximum number of mutations in a read. If -1.0, the mutation distribution is going to be used.')

    parser.add_argument('-W', '--workDir', type=str, dest='workDir', default=None,
        help='String, defautl: None. Directory to record the runs.')

    parser.add_argument('-s', '--seed', type=int, dest='seed', default=None,
        help='Integer, default: None. Random number seed.')

    parser.add_argument('-del', '--deletions', type=str, dest='deletions', default='0',
        help='String, default: 1. 0, ignore deletions. 1 if using deletions as mutation, 2 if delelting reads with deletions.')

    parser.add_argument('--no_threading',dest='threading',action='store_false',
        help='Default: False. Specify if clustering should not be threaded.')

    parser.add_argument('-j', '--remove_jackpotting', dest='remove_jackpotting', action='store_true',
        help='Default: False. Specify if the jackpotting needs to be corrected for.')

    parser.add_argument('-queue', '--queue', type=str, dest='queue', default='normal',
        help='String, default: <normal>. The queue to submit the job.')

    parser.add_argument('-recur', '--recur', type=str, dest='recur', default='no',
        help='String, default: "no". The queue to submit the job.')

    parser.add_argument('-trimUp', '--trimUp', type=int, dest='trimUp', default=400,
        help='Integer, default: 0. The number of bases to keep before the reference')

    parser.add_argument('-trimDown', '--trimDown', type=int, dest='trimDown', default=400,
        help='Integer, default: 0. The number of bases to keep after the reference')

    parser.add_argument('-dmuts', '--dmuts', type=int, dest='dmuts', default=3,
                        help='Integer, default: 3. The minimum number of bases between 2 mutations on a read.')


    args = parser.parse_args()

    if convert_args:
        if args.mask_index == 'no':
            args.mask_index = []
        else:
            args.mask_index = [int(x) for x in args.mask_index.split(',')]

    return args


