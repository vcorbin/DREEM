#!/usr/bin/python3

from EM_files import find_nread
import sys

filename = sys.argv[1]
count = find_nread(filename)
print(count)

