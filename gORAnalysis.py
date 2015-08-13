#!/usr/bin/env python

from sys import argv
from sys import exit
import os

# CONSTANTS
REF1 = "HVRI.fasta"
REF2 = "HVRII.fasta"
CONFIG = "config.txt"

DATA_LOC = "../GeneticsOfRace/"

def main():
    ref1 = readInRef(REF1)
    ref2 = readInRef(REF2)
    pairs = readInConfig()

    # Generate the FASTA files
    path = DATA_LOC
    fastaDict = makeFASTAs(path, pairs, ref1, ref2)


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


# path is the path to the directory containing to VCF files
# Returns a dictionary of {sampleId: newRef} where the newRef is
# the reference with the changes incorporated from the VCF file and
# the sampleId is the Id listed in the config file
def makeFASTAs(path, pairs, ref1, ref2):
    toReturn = {}

    # Get VCF files
    vcfFiles = [x for x in os.listdir(path) if x.endswith(".vcf")]

    # Get lists of ref1 and ref2 to be used
    firsts = pairs.keys()
    seconds = pairs.values()

    for vcf in vcfFiles:
        nub = vcf.split(".")[0]

        # Figure out which region the VCF file comes from
        if nub in firsts:
            newRef = makeChanges(ref1, vcf)
        elif nub in seconds:
            newRef = makeChanges(ref2, vcf)
        else:
            continue  # vcf file not in config file, skip it

        toReturn[nub] = newRef  

    return toReturn


def makeChanges(seq, vcf):
    print "MAKEING CHANGES"
        

        
if __name__ == '__main__':
    main()
