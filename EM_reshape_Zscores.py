#!/usr/bin/python
try:
    from . import EM_lib as lib
except:
    import EM_lib as lib
import numpy
import pandas


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Reshape the Z-scores matrix by adding 0s in front and behind the '
                                                 'current Z-scores matrix')
    parser.add_argument('input_filename', type=str, help='File containing the Z-scores matrix')
    parser.add_argument('-o', '--output', type=str, dest='output_filename', default=None,
                        help='File containing the new Z-scores matrix. By default, the same as input_filename.')
    parser.add_argument('-e', '--expand', type=str, dest='expand', default='0,0',
                        help='2 integer separated by a comma. Number of 0s to put in front and behind the correlation '
                             'matrices.')

    args = parser.parse_args()
    expand = numpy.array(args.expand.split(","), dtype=int)
    if args.output_filename is None:
        args.output_filename = args.input_filename

    df = pandas.read_csv(args.input_filename, sep=" ", header=None, dtype=float)
    new_df = lib.expand_df(df, expand[0], expand[1])
    new_df.to_csv(args.output_filename, sep=" ", header=False, index=False)

