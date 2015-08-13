#!/usr/bin/env python

from sys import argv

filew = open(argv[1] + ".plop", 'w')
with open(argv[1], 'r') as filer:
    for line in filer:
        if line.startswith(">"):
            filew.write(line)
        else:
            filew.write(line.strip())

filew.close()
