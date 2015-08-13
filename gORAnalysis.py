#!/usr/bin/env python

from sys import argv
from sys import exit

# CONSTANTS
REF1 = "HVRI.fasta"
REF2 = "HVRII.fasta"
CONFIG = "config.txt"

def main():
    ref1 = readInRef(REF1)
    ref2 = readInRef(REF2)
    config = readInConfig()

    print config


# Is given a file name and reads in the nucleotide seq as a string
# returns this string. 
def readInRef(fileName):
    toReturn = ''
    with open(fileName, 'r') as filer:
        for line in filer:
            if line.startswith(">"):
                continue
            else:
                toReturn += line.strip()

    return toReturn


# Reads in the config file into a dictionary so that
# the first entry is the key and the second entry is
# the value on each line. Returns a dictionary of this form
def readInConfig():
    dictio = {}

    with open(CONFIG, 'r') as filer:
        for line in filer:
            listl = line.split(",")
            listl = [x.strip() for x in listl]
            dictio[listl[0]] = listl[1]

    return dictio

        
if __name__ == '__main__':
    main()
