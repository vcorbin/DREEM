#!/usr/bin/python
import sys


if __name__ == "__main__":
    print(sys.argv)
    try:
        dir_i = sys.argv.index("-W")
        sys.argv[dir_i + 1] = 'test'
    except ValueError:
        sys.argv += ["-W","test2"]

    print(" ".join(sys.argv[1:]))

