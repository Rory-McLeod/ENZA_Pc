#!/usr/bin/env python
"""
Created on: 29-08-2016
@author: H.J.C. Cornelisse
Class is used to do main functions like checking the inputs and creating directories
Todo:
- add excecute for common use
"""
import os.path
import gzip
import sys
from time import strftime
import subprocess
from os import listdir
import logging

class Main:
    fastQFileList = []
    refGenomeList = []
    mapperClass = []
    assemblerClass = []
    bamClass = []
    primerDesClass = []
    visualisationClass = []
    genesFile = ""
    gffFile = ""
    workDir = ""
    resultDir = ""
    genomeAdd = ""
    fastQAdd = ""
    threadList = list()
    logger = logging.getLogger('WorkTitleHere')

    def __init__(self):
        """
        Todo:
        - Might have to remove this... don't see the use of this
        """
        hdlr = logging.FileHandler(Main.workDir + '/logging.log')
        formatter = logging.Formatter('%(asctime)s : %(threadName)s - %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        Main.logger.addHandler(hdlr)
        Main.logger.setLevel(logging.DEBUG)
        return

    def openFastQFiles(self, fastQlist):
        """
        Reads the list of fastQ files, and check if the input is correct. if needed, it unzips it.
        If good, saved to self.fastQFileList (accessible from all modules)
        Args:
            fastQlist: input list from user containing all fastQ files

        Returns:
            None: returns to the place of calling
        Todo:
        - Check if file exist in current directory
        """
        # divide files to seperate fastQ files
        Main.logger.debug("Checking fastQ files, please wait")
        fromMap = False
        if fastQlist is None:
            fromMap = True
            fastQlist = [f for f in listdir("Reads") if os.path.isfile(os.path.join("Reads", f))]
            fastQlist.sort()
        run = 0
        if fromMap:
            Main.fastQAdd = "Reads/"
        for fastQ in fastQlist:
            Main.logger.debug("Running "+fastQ)
            if ".gz" in fastQ and ".fastq" in fastQ:
                if run == 0:
                    Main.logger.info("GZip file(s) found, trying to unzip it, please wait")
                    run = 1
                outfilename = fastQ.replace(".gz", "")
                try:
                    Main.logger.debug("unzipping "+fastQ)
                    inF = gzip.open(Main.fastQAdd + fastQ, 'rb')
                    outF = open(Main.fastQAdd + outfilename, 'wb')
                    outF.write(inF.read())
                    inF.close()
                    outF.close()
                    self.fastQFileList.append(outfilename)
                except IOError:
                    Main.logger.error("Error unzipping file. please do this manually and rerun the program")
                    sys.exit(1)
            elif ".fastq" in fastQ or ".fq" in fastQ:
                self.fastQFileList.append(fastQ)
            else:
                Main.logger.error("File " + str(fastQ) + " does not have the correct extension. "
                                  "please check if file is fastq then rename extension to .fastq")
                sys.exit(1)
        return

    def openRefGenomes(self, refGenomeFileList):
        """
        Reads the list of genome files, and check if the input is correct. if needed, it unzips it.
        If good, saved to self.refGenomeList (accessible from all modules)
        Args:
            refGenomeFileList: input list from user containing all genome files

        Returns:
            None: returns to the place of calling
        """
        # select the reference genomes
        Main.logger.debug("Checking genome files, please wait")
        fromMap = False
        if refGenomeFileList is None:
            fromMap = True
            refGenomeFileList = [f for f in listdir("Genomes") if os.path.isfile(os.path.join("Genomes", f))]
            refGenomeFileList.sort()
        run = 0
        if fromMap:
            Main.genomeAdd = "Genomes/"
        for refGenome in refGenomeFileList:
            Main.logger.debug("Running "+refGenome)
            if ".fa" in refGenome or ".fsa_nt" in refGenome:
                if ".gz" in refGenome:
                    if run == 0:
                        Main.logger.info("GZip file(s) found, trying to unzip it, please wait")
                        run = 1
                    try:
                        Main.logger.debug("unzipping "+refGenome)
                        outfilename = refGenome.replace(".gz", "")
                        inF = gzip.open(Main.genomeAdd + refGenome, 'rb')
                        outF = open(Main.genomeAdd + outfilename, 'wb')
                        outF.write(inF.read())
                        inF.close()
                        outF.close()
                        self.refGenomeList.append(outfilename)
                    except IOError:
                        Main.logger.error("Error unzipping file. please do this manually and rerun the program")
                        sys.exit(1)
                else:
                    self.refGenomeList.append(refGenome)
            else:
                Main.logger.error("File " + str(refGenome) + " does not have the correct extension. please check if file is fasta" \
                                  " then rename extension to .fa or other fasta extensions")
                sys.exit(1)
        return

    @staticmethod
    def makeDirectory(directoryName):
        """
        Makes the temporary folders for workfiles and resultfiles.
        Args:
            directoryName: name of the directory that needs to be created

        Returns:
            None: returns to the place of calling
        """
        Main.logger.debug("Creating directory, please wait")
        if not os.path.exists(directoryName):
            try:
                os.makedirs(directoryName)
            except OSError:
                Main.logger.error("Could not generate directory " + directoryName +
                                  " please check if you are allowed to generate directorys")
                sys.exit(1)
        return

    @staticmethod
    def printer(line):
        """
        Ooh you gotta love the timestamps on your output!
        Args:
            line: line to be printed

        Returns:
            None: returns to the place of calling
        """
        timeSpot = strftime("%Y-%m-%d %H:%M:%S")
        print "[" + timeSpot + "]" + line
        return

    @staticmethod
    def execute(cmd, worktext="working, please wait"):
        Main.logger.debug("Running command: " + cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        Main.printer(worktext)
        jobNr, stderr = p.communicate()
        if p.returncode == 0:
            Main.logger.info("Completed command: " + cmd)
        else:
            Main.logger.error("Error running command: " + cmd)
            Main.logger.error(stderr)
        return

