#!/usr/bin/python
import pysam
import matplotlib
import argparse
from scipy.spatial.distance import hamming
matplotlib.use('agg')

def reduce_by_barcode(inFile, min_Ms, expected_start, prefix, len_barcode, barcode_start, collapse):
    sam = pysam.AlignmentFile(inFile, 'r')

    valid_queries, reverse_reads, g_map, bad_reads = set_valid_queries(sam, min_Ms, expected_start, prefix, len_barcode, barcode_start)
    sam.close()

    graph(inFile, g_map)

    forward_reads = [None]*len(reverse_reads)

    #Fill forward reads with pairs of reverse reads
    # Maybe there is a way to refresh sam.fetch() without reopening file?
    if collapse:
        forward_all = {}

        sam = pysam.AlignmentFile(inFile, 'r')
        for read in sam.fetch():
            if not read.is_reverse:
                forward_all[read.query_name] = read
            if (not read.is_reverse) and (read.query_name in valid_queries):
                forward_reads[valid_queries[read.query_name]] = read

        num_forward = 0
        new_file = pysam.AlignmentFile("collapsed_%s" % inFile, "wh", template=sam)
        for i in range(len(reverse_reads)):
            new_file.write(reverse_reads[i])
            if forward_reads[i]:
                num_forward += 1
                new_file.write(forward_reads[i])
        
        print( "Reverse Reads in output: %s" % len(reverse_reads))
        print( "Forward Reads in output: %s" % num_forward)
        new_file.close()

        bad_reads_file = pysam.AlignmentFile("removed_%s" % inFile, "wh", template=sam)
        for read in bad_reads:
            bad_reads_file.write(read)
            if read.query_name in forward_all:
                bad_reads_file.write(forward_all[read.query_name])
        bad_reads_file.close()
        sam.close()

def set_valid_queries(sam, min_Ms, expected_start, prefix, len_barcode, barcode_start):
    read_map = {}
    valid_queries = {}
    g_map = {}
    reverse_reads = []
    total = 0
    num_decent = 0
    bad_ref_starts = 0
    no_barcode = 0
    for read in sam.fetch():
        total += 1
        decent, bad_ref_starts = is_decent(read, min_Ms, expected_start, bad_ref_starts)
        if decent:
            num_decent += 1
            loc, barcode, gs = get_barcode(read, prefix, len_barcode, barcode_start)
            if barcode:
                if barcode in read_map:
                    read_map[barcode].append(read)
                else:
                    read_map[barcode] = [read]
                if gs in g_map:
                    g_map[gs] += 1
                else:
                    g_map[gs] = 1
            else:
                no_barcode += 1
    print "Total Reads in input: %s" % total
    print "Total Reversed Bad Ref Starts: %s" % bad_ref_starts
    print "Total Decent Reads in input (reversed, > %s matches and reasonable reference start): %s" % (min_Ms, num_decent)
    print "Total Decent with no barcode identified: %s" % no_barcode
    bad_reads = []
    i = 0
    non_duplicated_barcodes = 0
    for barcode in read_map:
        if len(read_map[barcode]) == 1:
            non_duplicated_barcodes += 1
            reads = read_map[barcode]
        else:
            reads = choose_query(read_map[barcode], expected_start, bad_reads)
        for read in reads:
            valid_queries[read.query_name] = i
            reverse_reads.append(read)
            i += 1
    print("Nondubplicated Barcodes:", non_duplicated_barcodes)
    print("Collapsed Reversed Read:", len(reverse_reads))
    print("Total number of duplicate-barcode reads that are <=1 away from most common: %s" % (sum(g_map.values()) - non_duplicated_barcodes - len(bad_reads)))
    print("Total number of duplicate-barcode reads that are >1 away from most common: %s" % len(bad_reads))
    print("Averge number of duplicate-barcode reads that are >1 away from most common: %s" % (len(bad_reads)/(len(reverse_reads) - non_duplicated_barcodes)))
    return valid_queries, reverse_reads, g_map, bad_reads

def is_decent(read, min_Ms, expected_start, bad_ref_start):
    if not read.is_reverse:
        return False, bad_ref_start
    if abs(read.reference_start - expected_start) > 10 and read.reference_start != 0:
        return False, bad_ref_start + 1
    matches = 0
    clips = 0
    for (op, count) in read.cigartuples:
        if op == 4:
            clips += count
        if op == 0:
            matches += count
    if matches > min_Ms:
        return True, bad_ref_start
    return False, bad_ref_start

# Get location of barcode wrt reference file and the barcode itself
# If can't find it or is in bad position, return None, None
def get_barcode(read, prefix, len_barcode, barcode_start):
    loc = read.seq.find(prefix)
    if loc != -1 and abs(loc + len(prefix) + read.reference_start - barcode_start) < 6 :
        loc = loc + len(prefix)
        return loc + read.reference_start, read.seq[loc: loc+len_barcode],  read.seq[loc+len_barcode: loc+len_barcode + 10]
    return None, None, None

# Choose which ever read(s) from list of reads with same barcode
def choose_query(reads, expected_start, bad_reads):
    return [most_common(reads, expected_start, bad_reads)]

def most_common(reads, expected_start, bad_reads):
    length = 350 #len(reads[0].seq) + expected_start
    read_strings = {}
    str_to_query = {}
    for read in reads:
        st = read_to_string(read, length)
        str_to_query[st] = read
        if st in read_strings:
            read_strings[st] += 1
        else:
            read_strings[st] = 1
    best_str = max(read_strings, key=read_strings.get)
    for st in read_strings:
        if hamming(list(best_str), list(st))*len(st) > 1:
            bad_reads.append(str_to_query[st])
    return str_to_query[best_str]

# using read.get_aligned_pairs(), turn the read into its relevant sequence, corresponding to positions in the reference.
# This is useful as opposed to just looking at read.seq, because the strings can now be compared using hamming distance.
def read_to_string(read, length):
    pairs = read.get_aligned_pairs()
    # print pairs
    bases = ['N']*length
    for pair in pairs:
        if pair[0] and pair[1]:
            bases[pair[1]] = read.seq[pair[0]]
    return "".join(bases)

def graph(inFile, dic, min_percent = 0.005):
    total = sum(dic.values()) + 0.0
    print("Total decent reversed reads with barcode:", total)
    print("Max with same barcode:", max(dic.values()))
    import operator
    sorted_dic = sorted(list(dic.items()), key=operator.itemgetter(1))[::-1]
    for key in sorted_dic:
        val = dic[key[0]]/total
        if val >= min_percent:
            print("%s \t %s" % (key[0], val))
        else:
            break


if __name__ == "__main__":
    print("Starting")

    parser = argparse.ArgumentParser()
    parser.add_argument('input_filename', type=str,
                        help='Preferably a binary_vector file (.txt) created by sam_to_binary(). If a sam/bam is given, '
                             'a binary_vector file of the same name will be created.')
    parser.add_argument('-m', '--min_Ms', type=int, dest='min_Ms', default=200, help='Minimum number of matches or mutations in read.')
    parser.add_argument('-s', '--start', type=int, dest='expected_start', default=37, help='Expected reference start of reverse read')
    parser.add_argument('-p', '--prefix', type=str, dest='prefix', default="TGACTAGCGGAGGCTAGAAGGAGAGAGA" ,help='Prefix to barcode')
    parser.add_argument('-l', '--len_barcode', type=int, dest='len_barcode', default=12, help='Length of barcode')
    parser.add_argument('-b', '--barcode_start', type=int, dest='barcode_start', default=281, help='Expected reference start of barcode')
    parser.add_argument('-c', '--collapse', type=bool, dest='collapse', default=False, help='True if should collapse reads')

    args = parser.parse_args()

    import sys
    sys.stdout = open("log_G_" + args.input_filename, 'w')
    print( sys.argv)

    reduce_by_barcode(args.input_filename, args.min_Ms, args.expected_start, args.prefix, args.len_barcode, args.barcode_start, args.collapse)

