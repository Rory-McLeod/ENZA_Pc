#!/usr/bin/env python

import os.path
import gzip
import sys
from time import gmtime, strftime


class Main:

    fastQFileList = []
    refGenomeList = []
    mapperClass = []
    bamClass = []

    def __init__(self):
        return

    def openFastQFiles(self, fastQlist):
        # divide files to seperate fastQ files
        Main.printer("Checking fastQ files, please wait")
        run = 0
        for fastQ in fastQlist:
            if ".gz" in fastQ and ".fastq" in fastQ:
                if run == 0:
                    print "GZip file(s) found, trying to unzip it, please wait"
                    run = 1
                outfilename = fastQ.replace(".gz", "")
                try:

                    inF = gzip.open(fastQ, 'rb')
                    outF = open(outfilename, 'wb')
                    outF.write(inF.read())
                    inF.close()
                    outF.close()
                    self.fastQFileList.append(outfilename)
                except IOError:
                    print "Error unzipping file. please do this manually and rerun the program"
                    sys.exit(1)
            elif ".fastq" in fastQ:
                self.fastQFileList.append(fastQ)
            else:
                print "File " + str(fastQ) + " does not have the correct extension. please check if file is fastq" \
                                             " then rename extension to .fastq"
                sys.exit(1)
        return

    def openRefGenomes(self, refGenomeFileList):
        # select the reference genomes
        Main.printer("Checking genome files, please wait")
        run = 0
        for refGenome in refGenomeFileList:
            if ".fa" in refGenome or ".fsa_nt" in refGenome:
                if ".gz" in refGenome:
                    if run == 0:
                        print "GZip file(s) found, trying to unzip it, please wait"
                        run = 1
                    try:
                        outfilename = refGenome.replace(".gz", "")
                        inF = gzip.open(refGenome, 'rb')
                        outF = open(outfilename, 'wb')
                        outF.write(inF.read())
                        inF.close()
                        outF.close()
                        self.refGenomeList.append(outfilename)
                    except IOError:
                        print "Error unzipping file. please do this manually and rerun the program"
                        sys.exit(1)
                else:
                    self.refGenomeList.append(refGenome)
            else:
                print "File " + str(refGenome) + " does not have the correct extension. please check if file is fasta" \
                                                 " then rename extension to .fa or other fasta extensions"
                sys.exit(1)
        return

    def makeDirectory(self, directoryName):
        Main.printer("Creating directory, please wait")
        if not os.path.exists(directoryName):
            try:
                os.makedirs(directoryName)
            except OSError:
                print "Could not generate directory " + directoryName + " please check if you are allowed to generate" \
                                                                        "directorys"
                sys.exit(1)
        return

    @staticmethod
    def printer(line):
        timeSpot = strftime("%Y-%m-%d %H:%M:%S")
        print "[" + timeSpot + "]" + line
        return
