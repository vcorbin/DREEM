#!/usr/bin/python3
import sys
try:
    from . import EM_classes as EMClass
except:
    import EM_classes as EMClass


#TODO fix the future warning that shows up everytime
def split_reponsibilities(filename):
    '''
    Given a file containing the responsibilities of each cluster for a bitvector, splits this into
    separate files containing only the reads that belongs to each cluster

    @param: filename = path to the input file containing the reponsibilities
    '''

    resp = EMClass.bitvect_resp(filename)
    resp.split_probabilistically()

if __name__ == "__main__":
    filename = sys.argv[1]
    split_reponsibilities(filename)
