#!/usr/bin/env python

from sys import argv
from sys import exit
import os

#############  CONSTANTS  #################
# Files
REF1 = "HVRI.fasta"
REF2 = "HVRII.fasta"
CONFIG = "config3.txt"
DATA_LOC = "../GeneticsOfRace/"

# Indexes
POS_IDX = 1
REF_IDX = 3
ALT_IDX = 4
FILTER_IDX = 6

# Strings
PASS = "PASS"

# Numerical Constants
HVRI_OFFSET = 15951
##########################################


# CUSTOM EXCEPTIONS
class UserException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return str(self.msg)

class InternalException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)


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
            newRef = makeChanges(ref1, vcf, path, True)
        elif nub in seconds:
            newRef = makeChanges(ref2, vcf, path, False)
        else:
            continue  # vcf file not in config file, skip it


        #print vcf
        toReturn[nub] = newRef  
        #print newRef
        #print

    return toReturn


# If there are multiple alleles for the same variant, it makes the greedy
# choice and choses the first one by default
# When the variant file is from HVRI, need to throw out records below a
# threshhold
# SNPs prioritized over INDELs
def makeChanges(seq, vcf, path, isHVRI):
    indelStack = []
    validSnpPoses = []
    nucList = [x for x in seq]

    with open(os.path.join(path, vcf), 'r') as filer:
        for line in filer:
            if line.startswith("#"):  # Skip comment lines
                continue

            listL = line.split("\t")
            # make everything upper case so no case clashes
            listL = [x.upper() for x in listL]

            # position should be an int not a string
            listL[POS_IDX] = int(listL[POS_IDX])

            if listL[FILTER_IDX] != PASS:  # Skip vars that don't pass filter
                continue

            pos = listL[POS_IDX]
            refNuc = oneInd(nucList, pos)
            vcfOldNuc = listL[REF_IDX]
            vcfNewNuc = listL[ALT_IDX]

            # More than one allele, just grab first one
            if "," in vcfNewNuc:
                vcfNewNuc = vcfNewNuc.split(",")[0]
                listL[ALT_IDX] = vcfNewNuc


            # These two lines protect against variants from other regions
            # that end up in these files for some reason 
            if isHVRI:   # Adjust for offset within file
                pos = pos - HVRI_OFFSET
                refNuc = oneInd(nucList, pos)  # Recalc refNuc
                listL[POS_IDX] = pos
                if pos < 1:
                    continue
                

            if refNuc is None:
                continue

            ########

            # Save Indels for later
            if len(vcfOldNuc) != 1 or len(vcfNewNuc) != 1:
                indelStack.append(listL) 
                continue

            # VCF doesn't match reference, something wrong
            if vcfOldNuc != refNuc:
                raise UserException("VCF REF entry doesn't match" +\
                        " NCBI reference sequence. \nVCF: " + listL[REF_IDX] +\
                        "\nREF: " + oneInd(nucList, pos) + "\nVCF_POS: " +\
                        str(pos) + "\nSee VCF_File: " + vcf)

            # Replace SNP
            nucList[pos - 1] = vcfNewNuc
            validSnpPoses.append(pos)

    return addIndels(nucList, indelStack, validSnpPoses, vcf)



# Returns the entry in the list at position pos - 1. It is used
# to adjust between 1-indexed systems and 0-indexed systems
def oneInd(ls, pos):
    try:
        return ls[pos - 1]
    except IndexError:
        return None


# Inputs a list of characters that represenets a nucleotide sequence
# a list of Indel intries and a list of positions not to add Indels in
# as they were already included as SNPs.
# Assumes position in each indel list is already an int (not a string)
# Assumes all multiple alleles have been replaced by a single one (in
# alternate)
def addIndels(nucList, indels, exceptions, fileName):
    isDel = False
    # sort indels by position (should already be done, but to be safe)
    indels.sort(key=lambda x: x[POS_IDX])
    indels = [x for x in indels if x[POS_IDX] not in exceptions]

    for variant in indels:
        # insertion
        if len(variant[REF_IDX]) <= len(variant[ALT_IDX]):
            nucList[variant[POS_IDX] - 1] = variant[ALT_IDX]
            print "INSERTION"
        else:  # deletion
            # no. nucleotides deleted
            print "DELETION", fileName
            isDel = True
            idx = variant[POS_IDX] - 1
            diff = len(variant[REF_IDX]) - len(variant[ALT_IDX])
            for i in range(idx + 1, idx + diff + 1):
                nucList[i] = '-'  # Dashes mark a deletion

    nucString = "".join(nucList)
    nucString = nucString.replace('-', '')

    return nucString

        
if __name__ == '__main__':
    main()
