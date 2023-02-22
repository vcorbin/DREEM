#!/usr/bin/python

import pandas as pd
import sys
import os
from coordinates import get_coordinates
from make_bitvector import make_bitvector_filename
from EM_classes import split_filename
import argparse


def main(config_file, date=None, sample_ID=None, reference=None, max_its=300, queue='normal', nreads=1000000,
         na_threshold=0.2):
    """
    Run all the clusterings specified in the optiond
    :param config_file: required. Name of the config file. Needs to be an excel spreadsheet
    :param date: str (year-month-day). Default: None. Date of the sample to be clustered. Multiple dates can be passed
                 as a list. If None, no filter is applied.
    :param sample_ID: str. Default: None. Name of the sample to cluster. The name of the sample correspond to a fastq
                      file. Multiple dates can be passed as a list. If None no filter is applied.
    :param reference: str. Default: None. Filter the clustering by name of the reference. It correspond to a reference
                      file. Multiple references can be passed as a list. If None no filter is applied. 
    :return: Void
    """
    # files paths
    workdir = "/lab/solexa_rouskin/projects/HIV"
    ref_path = workdir + "/refs"
    fastq_path = workdir + "/data/fastq"
    bam_path = workdir + "/data/bam"
    bitvector_path = workdir + "/results/clustering"

    config_header = ['Date', 'Seq ID', 'Reference', 'Forward P', 'Reverse P']
    new_header = ['Date', 'Seq ID', 'Reference', 'Forward P', 'Reverse P']

    df = pd.read_excel(config_file)
    assert set(config_header).issubset(df.columns), \
        "Error: the config file doesn't contain the right columns: %s" % ", ".join(config_header)
    df = df[config_header]
    df.columns = new_header  # This is so that in case the names of the columns change, we don't need to change the code

    # Remove the line without reference, sample ID or primers
    df = df[df['Seq ID'].notnull()]  # sample ID
    df = df[df['Reference'].notnull()]  # reference
    df = df[df['Forward P'].notnull()]  # forward primer
    df = df[df['Reverse P'].notnull()]  # reverse primer

    # Filter by date
    df = __df_filter(df, 'Date', date)
    # Filter by sample_ID
    df = __df_filter(df, 'Seq ID', sample_ID)
    # Filter by reference
    df = __df_filter(df, 'Reference', reference)

    # Transform into a dictionary
    df_dict = __df_to_dic(df)

    # Loop over the sample_names,references and align
    # then loop over the primers and run make_bitvector.py and parallel.py
    for sample_ref in df_dict:
        sample_ID, ref = sample_ref.split(',')
        bowtie2_command, sam_file, ref_file = command_alignment(fastq_path, sample_ID, ref_path, ref, bam_path)
        bam_file = split_filename(sam_file).noext + ".bam"

        if os.path.isfile(bam_file):
            print("%s already exists, skipping the alignment." % bam_file)
            print_command = 'bsub -J ' + bam_file + ' -q ' + queue + ' echo "%s already exists, skipping the alignment."' % bam_file
            os.system(print_command)
        else:
            print("Aligning %s using %s as reference." % (sample_ID, ref))
            bowtie2_command = 'bsub -J ' + sam_file + ' -q ' + queue + ' ' + bowtie2_command
            print("aligning")
            os.system(bowtie2_command)
            bam_command = "bsub -w 'ended(\"" + sam_file + "\")' -J " + bam_file + " -q " + queue + " samtools view -bS -o " \
                          + bam_file + " " + sam_file
            os.system(bam_command)
            remove_sam_command = "bsub -w 'ended(\"" + bam_file + "\")' rm " + sam_file
            os.system(remove_sam_command)

        for primers in df_dict[sample_ref]:
            primers_str = ",".join(primers)
            # get coordinates form primers and reference
            # make the bitvector file name
            bit_command, temp_bit_file = command_makebitvector(bam_file, ref_file, bitvector_path, primers_str)
            bit_file = split_filename(temp_bit_file)
            bit_file = bit_file.path + "/" + bit_file.nopath.split("temp_")[1]
            if os.path.isfile(bit_file):
                print("%s already exists, skipping the bit-vectors making." % bit_file)
                print_command = 'bsub -J ' + bit_file + ' -q ' + queue + ' echo "%s already exists, skipping the making of bitvectors."'\
                                % bit_file
                os.system(print_command)
            else:
                bit_command = "bsub -w 'ended(\"" + bam_file + "\")' -q " + queue + " -J " + temp_bit_file + " " + bit_command
                os.system(bit_command)
                mv_bit_command = "bsub -w 'ended(\"" + temp_bit_file + "\")' -q " + queue + " -J " + bit_file + " mv " + temp_bit_file \
                                 + " " + bit_file
                os.system(mv_bit_command)

            # clustering
            cluster_command = "bsub -w 'ended(\"" + bit_file + "\")' parallel.py " + bit_file + " --recur yes -i " + \
                              str(max_its) + " --queue " + queue + " -r " + str(nreads) + " -q " + str(na_threshold)
            os.system(cluster_command)


def command_alignment(fastq_path, sample_ID, ref_path, ref, bam_path):
    """
    Create the command for the alignment of the fastq file
    :param fastq_path: str, path of the fastq file
    :param sample_ID: str, name of the sample, which indicates what the name of the fastq file will be
    :param ref_path: str, path of the reference file
    :param ref: str, name of the reference. Must correspond to the name of the index file (without the numeral and 
    extension) 
    :param bam_path: str, path of the location where the bam file will be created
    :return: a bowtie2 command, the name of the bam_file and the name of the reference file
    """
    bowtie2_command = ""
    ref_file_noext = ref_path + "/" + ref
    fastq_file_1 = fastq_path + "/" + sample_ID + '_1_sequence.fastq'
    fastq_file_2 = fastq_path + "/" + sample_ID + '_2_sequence.fastq'
    with open(fastq_file_1,'r') as f:
        line_1 = f.readline()
        line_2 = f.readline()
        line_3 = f.readline()
        line_4 = f.readline()  # This is the line for quality scores
    f.close()
    if any(letter.islower() for letter in line_4):
        qual_encoding = '--phred64'
        ascii_encoding = 'Phred+64'
    else:
        qual_encoding = ''
        ascii_encoding = 'Phred+33'
    print("Input qualities are encoded with: " + ascii_encoding)
    sam_file = bam_path + '/' + sample_ID + '_' + ref + '.sam'
    bowtie2_command = 'bowtie2 -L 15 --local --no-unal -x ' + ref_file_noext + ' -1 ' + fastq_file_1 + ' -2 ' \
                      + fastq_file_2 + ' -S ' + sam_file + ' ' + qual_encoding

    return bowtie2_command, sam_file, ref_file_noext + '.fasta'


def command_makebitvector(bam_file, ref_file, bitvector_path, primers, nreads=1000000):
    """
    Create the command for the creation of the bitvector file.
    :param bam_file: str, name of the bam file
    :param ref_file: str, name of the reference file
    :param bitvector_path: str, path of the location where the bitvector is going to be created
    :param primers: str, the forward and reverse primers delimited by ','
    :param nreads: int, maximum number of reads to analyze
    :return: a make_bitvector command, the name of the bit_vector file
    """
    primer_1_seq, primer_2_seq = primers.split(",")
    ref_coords, ref_seq, ref_name = get_coordinates(primer_1_seq, primer_2_seq, ref_file)
    bit_file = make_bitvector_filename(bam_file, ref_coords)
    bit_file = bitvector_path + '/temp_' + split_filename(bit_file).nopath
    command = 'make_bitvector.py -f ' + bam_file + ' -n ' + str(nreads) + ' -r ' + ref_file + ' -c ' + ",".join([str(x) for x in ref_coords])\
              + ' -o ' + bit_file
    return command, bit_file


def __df_filter(df, column, filter):
    """
    filter the dataframe according to values specified by filter
    :param df: pandas.DataFrame
    :param column: str, name of the column to be used for filtering
    :param filter: string or list of string
    :return: filtered dataframe
    """
    if filter is not None:
        assert column in df.columns, "Error: there is no column in the dataframe to be filtered that is named '%s'. " \
                                     "The dataframe columns are namesd: %s" % (column, ", ".join(df.columns))
        assert type(filter) is list or type(filter) is str, "Error: the filter for '%s' must be either a list or a " \
                                                            "string." % column
        if type(filter) is str:
            filter = [filter]
        else:
            pass
        df = df[df[column].isin(filter)]
        if df.empty:
            print("The values specified in filter are not present in the column '%s' of the dataframe." % column)
            sys.exit()
    else:
        pass
    return df


def __df_to_dic(df):
    """
    turn the dataframe into a dictionary to be looped over for the alignment, make_bitvector and clustering
    :param df: pandas.DataFrame
    :param field_names: list of the name of the column to be included in the dictionary: sample_ID, reference, 
    forward_primer, reverse_primer
    :return: a dictionary of dictionary with the sample IDs as keys and the reference as subkeys
    """
    df = df.reset_index()
    df_dict = {}
    for i in range(df.shape[0]):
        sample_ID = df['Seq ID'][i]
        reference = df['Reference'][i]
        for_p = df['Forward P'][i]
        rev_p = df['Reverse P'][i]
        try:
            df_dict[sample_ID + ',' + reference].append([for_p, rev_p])
        except KeyError:
            df_dict[sample_ID + ',' + reference] = [[for_p, rev_p]]
    return df_dict


def __convert_filter(filter):
    """
    Take a filter passed in the command line and parse it.
    :param filter: str (or None). series of filters delimited by ','
    :return: a list, containing the filters (or None)
    """
    if filter is not None:
        return filter.split(",")
    else:
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='From a specified config file, align fastq file and cluster. Lines in '
                                                 'the config file that are missing information are ignored.')

    parser.add_argument('inFile', type=str,
                        help='String, required. Config file containing the information for each sample')

    parser.add_argument('-i', '--max_its', type=int, dest='max_its', default=3000,
                        help='Integer, default: 3000. Number of iterations allowed in an EM run before it is stopped. '
                             'The run will stop prior to this if convergence is obtained.')

    parser.add_argument('--date', type=str, dest='date', default=None,
                        help='str. Filter the samples to be analysed by their dates(year-month-day). '
                             'If multiple dates, separate them with a comma, e.g.:"2018-02-08,2017-11-03". '
                             'If not specified, dates will not be used as a filter.')

    parser.add_argument('--sample_ID', type=str, dest='sample_ID', default=None,
                        help='str. filter the samples to be analysed by their names. If multiple dates, separate them '
                             'with a comma, e.g.:"180206Rou_D18-1064,180206Rou_D18-1065". '
                             'If not specified, the samples will not be filtered by their names.')

    parser.add_argument('--ref', type=str, dest='ref', default=None,
                        help='str. Filter the samples to be analysed by their alignment reference. '
                             'If multiple references, separate them with a comma, e.g.:"RREshort,NL43rna".'
                             'If not specified, references will not be used as a filter')

    parser.add_argument('-queue', '--queue', type=str, dest='queue', default='normal',
                        help='String, default: <normal>. The queue to submit the job.')

    parser.add_argument('-n', '--nreads', type=int, dest='nreads', default=1000000,
                        help='Integer, default: 1 million. Maximum number of reads to analyse.')

    parser.add_argument('-q','--q_threshold',type=float,dest='na_threshold',default=0.2,
                        help='Float, default: 0.2. Threshold of the ratio of "?", "." and "N" a read from bit-vector '
                             'file is allowed to have to be loaded.')

    args = parser.parse_args()

    # Convert the filters into lists:
    config_file = args.inFile
    max_its = args.max_its
    queue = args.queue
    nreads = args.nreads
    na_threshold = args.na_threshold
    date = __convert_filter(args.date)
    sample_ID = __convert_filter(args.sample_ID)
    reference = __convert_filter(args.ref)


    main(config_file, date=date, sample_ID=sample_ID, reference=reference, max_its=max_its, queue=queue, nreads=nreads,
         na_threshold=na_threshold)


# Tests:
def test___df_filter():
    df = pd.DataFrame({'a': ['1', '1', '2', '3', '3'], 'b': ['a', 'b', 'd', 'c', 'c']})
    test_df = __df_filter(df, 'a', '3').reset_index(drop=True)
    assert test_df.equals(pd.DataFrame({'a': ['3', '3'], 'b': ['c', 'c']})), "Error: df_filter() test1 failed."
    test_df = __df_filter(df, 'a', ['3']).reset_index(drop=True)
    assert test_df.equals(pd.DataFrame({'a': ['3', '3'], 'b': ['c', 'c']})), "Error: df_filter() test2 failed."
    test_df = __df_filter(df, 'b', ['a', 'c']).reset_index(drop=True)
    assert test_df.equals(pd.DataFrame({'a': ['1', '3', '3'], 'b': ['a', 'c', 'c']})), \
        "Error: df_filter() test3 failed."
    print("__df_filter() has been successfully tested.")


def test_command_alignment():
    fastq_path = 'fastq_path_test/test1'
    sample_ID = 'sample_test'
    ref_path = 'ref_test/test2/test3/'
    ref = 'reference'
    bam_path = '/users/wherever/'
    bowtie2_command = 'bowtie2 -L 12 --local --no -unal -x ref_test/test2/test3//reference ' \
                      '-1 fastq_path_test/test1/sample_test_1_sequence.fastq ' \
                      '-2 fastq_path_test/test1/sample_test_2_sequence.fastq ' \
                      '| samtools view -bS -o /users/wherever//sample_test_reference.bam -'
    bam_file = '/users/wherever//sample_test_reference.bam'
    ref_file = 'ref_test/test2/test3//reference.fasta'
    assert command_alignment(fastq_path, sample_ID, ref_path, ref, bam_path) == (bowtie2_command, bam_file, ref_file), \
        'Error: command_alignment() test failed.'
    print ("command_alignment() has been successfully tested.")


def test_command_makebitvector():
    bam_file = 'bampath/bamfile.bam'
    ref_file = 'refpath/reffile.fasta'
    bitvector_path = 'bitpath'
    primers = 'sfsd,erer'













































