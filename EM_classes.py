import pandas
import numpy
import random
try:
    from . import EM_lib as lib
    from . import EM_files as files
except:
    import EM_lib as lib
    import EM_files as files
import os
from EM_lib import bitVect_to_bitString
from EM_lib import bitString_to_bitVect


class single_bitvector(object):
    """
    Object containing the information found in a line of a bitvector file.
    """
    def __init__(self, bitstring='', query_name='-', reference="Ref", index=""):
        self.bitstring = bitstring
        self.index = index
        self.query_name = query_name
        self.reference = reference
        self.N_deletions = bitstring.count('D')
        self.coverage = len(bitstring) - bitstring.count('.')
        self.bitvect, self.N_mutation = lib.bitString_to_bitVect(bitstring, 0)
        self.dict = {'Query_name': self.query_name, 'Bases_vector': self.bitstring, 'N_mutations': self.N_mutation,
                     'N_deletions': self.N_deletions, 'Coverage': self.coverage, 'Reference': self.reference,
                     'Index': self.index}
        self.header = ['Query_name', 'Bases_vector', 'N_mutations', 'N_deletions', 'Coverage', 'Reference', 'Index']


    def set_query_name(self, query_name):
        assert type(query_name) is str, "Error: query_name must be a string, not a %s" % type(query_name)
        self.query_name = query_name
        self.dict['Query_name'] = self.query_name


    def set_reference(self, reference):
        assert type(reference) is str, "Error: query_name must be a string, not a %s" % type(reference)
        self.reference = reference
        self.dict['Reference'] = self.reference


    def get_string(self):
        line_array = [str(self.dict[key]) for key in self.header]
        return self.__format_line(line_array)


    def get_header_string(self):
        return self.__format_line(self.header)


    def __format_line(self, x): # could assert they are all strings, but can't be bothered
        return "\t".join(x) + "\n"


class ind_bit_vect(object):
    '''
    Attributes:
    bitStrings - if bitStrings is specified it will be used to order self.N and self.bitMatrix
    N - array containing the number of occurrences of each read_csv
    bitMatrix - array of individual reads, which counts in the corresponding position of N
     
    '''
    def __init__(self, bit_vectors, N_occur, bitStrings=None, refs_info=None, filename=None):

        assert type(bit_vectors) in [dict, pandas.core.frame.DataFrame], "bit_vectors must be either a dict or a " \
                                                                         "pandas.core.frame.DataFrame."
        assert type(N_occur) in [dict, numpy.ndarray], "N_occur must be either a dict or a numpy.ndarray."

        if type(N_occur) == dict or type(bit_vectors) == dict:
            assert type(N_occur) == type(bit_vectors), "If one of N_occur or bit_vector is a dict the other must be too"


        assert type(refs_info) is comments or refs_info is None, "refs_info must be an object of class 'comments' or " \
                                                                 "must be 'None'"

        if bitStrings is None:
            assert type(N_occur) == dict, "N_occur must be a dict since bitStrings is not provided."
            keys = N_occur.keys()  # Get the keys from N_occur dictionary to make sure we preserve the same order.
            self.bitStrings = numpy.asarray(list(keys))
        else:
            assert type(bitStrings) == numpy.ndarray, "bitStrings must be a numpy.ndarray."
            self.bitStrings = bitStrings
            if type(N_occur) == dict:
                assert set(N_occur.keys()) == set(bitStrings), "All the keys of N_occur must be present in bitString."
                assert len(N_occur.keys()) == len(bitStrings), "bitString can't contain the same string twice."
                keys = bitStrings

        if type(bit_vectors) == dict:
            self.bitMatrix = pandas.DataFrame(numpy.asarray([bit_vectors[key] for key in keys]))
        else:
            self.bitMatrix = bit_vectors

        if type(N_occur) == dict:
            self.N = numpy.asarray([N_occur[key] for key in keys])
        else:
            self.N = N_occur

        assert self.bitMatrix.shape[0] == self.N.shape[0], "Error: Number of reads in bitMatrix is different that " \
                                                           "the number of reads in N (%d vs %d)." \
                                                           % (self.bitMatrix.shape[0], self.N.shape[0])
        assert self.bitMatrix.shape[0] == self.bitStrings.shape[0], "Error: Number of bitStrings differ from the " \
                                                                    "number of read in bitMatrix (%d vs %d)." % \
                                                                    (self.bitMatrix.shape[0], self.bitStrings.shape[0])

        if type(refs_info) is not comments:
            self.refs_info = comments(None)
        else:
            self.refs_info = refs_info
        self.filename = filename

    def get_pop_average(self):
        pop_average = (self.bitMatrix.T * self.N).sum(axis=1) / self.get_N_bases()
        return pop_average.values

    def get_std(self):
        std = numpy.sqrt((numpy.square(self.bitMatrix - self.get_pop_average()).T * self.N).sum(axis=1)) /\
                                                                                        numpy.sqrt(self.get_N_bases())
        return std.values

    def get_Ntotal(self):
        return self.N.sum()

    def get_N_bases(self):
        """
        Calculate the total number of reads for each base (na are skipped). If no na are present, it will equal Ntotal 
        for each base.
        :return: numpy.ndarray of integer which length is the number of bases in the sequence
        """
        return (self.N * numpy.isfinite(self.bitMatrix.values).T).sum(axis=1)

    def copy(self):
        return ind_bit_vect(self.bitMatrix, self.N, self.bitStrings)

    def collapse_similar_reads(self):
        """
        When the bitMatrix contain reads that are the same (the mutation are in the same location but the bases 
        are different) they are combined and N is modified accordingly
        :return: a new ind_bit_vect object
        """
        bitdic = {}
        for i, read in enumerate(self.bitMatrix.values):
            read = bitVect_to_bitString(read)
            try:
                bitdic[read] += self.N[i]
            except KeyError:
                bitdic[read] = self.N[i]
        reads, N = zip(*bitdic.items())
        bitvect = numpy.asarray([bitString_to_bitVect(read)[0] for read in reads])
        return ind_bit_vect(pandas.DataFrame(bitvect), numpy.asarray(N), bitStrings=numpy.asarray(reads),
                            refs_info=self.refs_info, filename=self.filename)

#
#    def save(self, filename, labels, start=0):
#        Nbitvect = x.N.sum()
#        Query_names = map(str,range(Nbitvect))
#        starts = map(str,[start] * Nbitvect)
#        N_muts = self.bitMatrix.sum(axis=1).astype(int).astype(str)
#
#        pandas.DataFrame({'Query_name':Query_names,'Binary_vector':self.bitStrings,'N_mutations':N_muts, \
#                          'Reference_name':labels['refs'],'Start_position':starts}).to_csv(filename,sep='\n')
#
#    def regenerate_bitStrings(self):
#        new_bitStrings = self.bitMatrix.replace([0.0, 1.0, numpy.nan], ['0', '1', '?']).values
#        self.bitStrings = map(''.join, new_bitStrings)


class bitvect_resp(object):
    """
    Load a responsibility file or dataframe into an object which can then be separated probabilistically into the
    individual clusters and transformed into ind_bit_vect objects.
    """
    def __init__(self, resp_df):
        if type(resp_df) == str:
            self.filename = resp_df
            resp_df = pandas.read_csv(resp_df, sep="\t", index_col=0)
        else:
            self.filename = None

        ncol = resp_df.shape[1]
        self.cluster_names = list(resp_df)[:ncol-2]
        self.Ncluster = ncol - 2
        self.responsibilities = numpy.array(resp_df.values[:,0:ncol-2],dtype=float)
        self.N = resp_df.values[:,ncol-2].astype(numpy.int64, copy=False)
        self.bitstr = resp_df.values[:,ncol-1]


    def split_probabilistically(self):
        """
        Split into bitvect_resp objects according to clusters. The reads are assigned to the clusters probabilistically.
        create a bitvector file for each cluster
        :return: an array of bitvect_resp object. One for each clusters.
        """
        if self.filename is not None:
            prefix = split_filename(self.filename).noext
        else:
            prefix = ''
        files = [open(prefix + '_' + name + '_bitvector.txt', 'w') for name in self.cluster_names]
        # write the header
        header = single_bitvector().get_header_string()
        for file in files:
            file.write(header)
        for i, cluster_resps in enumerate(self.responsibilities):
            #print(self.bitstr[i])
            bitvect_line = single_bitvector(self.bitstr[i])
            for j in range(self.N[i]):
                bitvect_line.set_query_name(str(i) + '_' + str(j))
                outstr = bitvect_line.get_string()
                prob = random.random()
                cum_prob = 0
                for j, resp in enumerate(cluster_resps):
                    cum_prob += resp
                    if prob <= cum_prob:
                        files[j].write(outstr)
                        break
                    else:
                        pass
        for file in files:
            file.close()


    def to_indbitvect(self, deletions='1'):
        bit_vectors = {}
        N_occur = {}
        for i, bitString in enumerate(self.bitstr):
            bitVect = lib.bitString_to_bitVect(bitString, deletions)[0]
#            if bitString in bit_vectors.keys():
#                print "Warning: %s has already occured in resp. N=%s N_new=%s" % (bitString, N_occur[bitString],
#                                                                                  self.N[i])
            bit_vectors[bitString] = bitVect
            N_occur[bitString] = self.N[i]
        return ind_bit_vect(bit_vectors, N_occur, bitStrings=self.bitstr)


class jackPotting(object):
    def __init__(self, N_theory, Ns, nmuts, bitStrings):
        '''
        What is N_theory?
        '''
        ordind = N_theory.argsort()
        self.N_theory = N_theory[ordind]
        self.Ns = Ns[ordind]
        self.nmuts = nmuts[ordind]
        self.bitStrings = bitStrings[ordind]
        minisone = self.N_theory.copy()
        minisone[minisone < 1.0] = 1.0
        self.N_dist = (self.Ns - minisone) / numpy.sqrt(minisone)

    def to_csv(self,filename):
        pandas.DataFrame(list(zip(self.N_theory,self.Ns,self.N_dist,self.nmuts,self.bitStrings))).to_csv(
            filename, sep="\t", header=['N_theory', 'N_obs', 'N_sd', 'nmut', 'bitString'], index_label=False
        )


class comments(object):
    def __init__(self, fileObj):
        """
        Create an object that contain the comments from a bitvector file (and subsequent cluster file) in the form of a
        a string and dictionary containing the useful information. The file object is moved to the header.
        :param fileObj: file object from which the comments are going to be extracted. None will create an object with
        self.txt = '' and self.refs = {'Ref': {'start':0, 'end':0, 'seq'=''}}, self.header = '' and sel.seq_length = 0.
        """
        self.text = ''
        if fileObj is None:
            self.header = ''
            info = {}
        else:
            assert fileObj.tell() == 0, "Error: The file object pointer does not point to the beginning of the file."
            info = {}
            line = fileObj.readline()
            while line[0] == '@':
                self.text += line
                try:
                    key, ref, value = line.strip('@\n').split("\t")
                except ValueError:
                    info = {}
                    break

                if ref not in info.keys():
                    info[ref] = {}
                else:
                    pass
                info[ref][key] = value
                line = fileObj.readline()
            self.header = line
        self.refs, self.seq_length = self._info_to_refs(info)

    def _info_to_refs(self, info):
        """
        Convert the info dictionary scrapped from the comment section into a dictionary that contains only field of
        interest (ref_start, ref_end, ref_sequence)
        :return: a dictionary containing a dictionary for each reference. the sub-dictionary contain the fields 'start',
        'end', and 'seq', and the length of the sequences (all the length should be the same).
        """
        refs = {}
        if info == {}:
            length = 0
            refs['Ref'] = {'start': 0,'end': 0 + length,'seq': ''}
        else:
            lengths = []
            for ref in info.keys():
                refs[ref] = {'start': int(info[ref]['start']),
                             'end': int(info[ref]['start']) + int(info[ref]['length']),
                             'seq': info[ref]['seq']}
                lengths.append(int(info[ref]['length']))
            length = lengths[0]
            assert all(x == length for x in lengths), "The lengths of the refs sequences differ"
        return refs, length

    def get_coordinates(self):
        """
        :return: a list of 3 elements: a string corresponding to the names of the gene(s) and  2 int corresponding to 
        the coordinates of the reference sequence (starting at 0 in case there are more than 1 ref).  
        """

        if self.refs is None:
            gene = 'ref'
            start = 0
        elif len(self.refs.keys()) > 1:
            gene = '_'.join(self.refs.keys())
            start = 0
        else:
            gene = list(self.refs.keys())[0]
            start = self.refs[gene]['start']
        return gene, start

    def get_sequence(self):
        """
        :return: a sequence if all sequences are the same, otherwise return an error.
        """
        seq = []
        for ref in self.refs.keys():
            seq.append(self.refs[ref]['seq'])
        assert(len(set(seq)) == 1), "There is more than 1 sequence in the header."
        return seq[0]


class clusterGroup(object):
    def __init__(self, filename):
        inF, header, self.cluster_info = files.load_header(filename)  # self.cluster_info is of class comments
        inF.close()
        self.cluster_names = header.strip().split("\t")[1:]
        self.df = pandas.read_csv(filename, sep="\t", index_col=0,comment="@")
        self.df.index = numpy.array(self.df.index).astype(int)
        self.filename = filename

    def aslist(self):
        """
        Create a list of lists out of the dataFrame
        :param self: Object of class clusterGroup
        :return: List of lists
        """
        return self.df.T.values

    def get_sequence(self):
        return self.cluster_info.get_sequence()


class singleCluster(object):
    def __init__(self, clusters, cluster_name):
        """
        Create an object containing the information for a single cluster
        :param clusters: either a clusterGroup object or a filename
        :param cluster_name: name of the cluster we want to load (name in the header in the cluster file)
        """
        assert type(clusters) is str or type(clusters) is clusterGroup, "Error: clusters must either be a filename" \
                                                                        " or a clusterGroup object"
        assert type(cluster_name) is str, "Error: cluster_name must be a string"
        if type(clusters) is str:
            clusters = clusterGroup(clusters)
        assert cluster_name in clusters.cluster_names, "Error: the cluster '%s' is not present in the cluster file '%s'"\
                                                       %(cluster_name, clusters.filename)

        self.cluster_info = clusters.cluster_info
        self.name = cluster_name
        self.filename = clusters.filename
        self.probs = clusters.df[cluster_name]

    def get_sequence(self):
        return self.cluster_info.get_sequence()


class split_filename(object):
    def __init__(self, filename):
        """
        From a filename, create an object which contain the path, the basename and the extension of the filename
        :param filename: str, name of the file of interest
        """
        self.fullname = filename
        self.basename = os.path.splitext(os.path.basename(filename))[0]  # no extension and no path. differs from
                                                                         # os.path.basename
        self.noext, self.ext = os.path.splitext(filename)
        self.nopath = self.basename + self.ext
        self.path = os.path.dirname(filename)
        self.abspath = os.path.abspath(filename)

















































