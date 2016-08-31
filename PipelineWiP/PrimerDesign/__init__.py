#!/usr/bin/env python

"""
Created on: 29-08-2016
@author: H.J.C. Cornelisse
Class is used to generate points of interessed from the mapped reads.
Todo:
- Testing POI generation (currently not used in this forServerRun.py)
- Check the generated POI, are they covered in all three reads?
"""

import subprocess
import threading
import sys
from Main import Main
import math
from collections import defaultdict
import numpy

class PrimerDesign(threading.Thread):

    def __init__(self, bamFile):
        """
        Method for initiating the primer design.
        threading.Thread is called for threaded use of this class
        Args:
            bamFile: raw file name to correspond to the correct files
        """
        threading.Thread.__init__(self)
        self.bamFile = bamFile
        return

    def run(self):
        """
        Method used for threaded use. the method mentioned below are run in a single thread.
        Returns:
            None: returns to the place of calling (start())
        """
        self.readsOnGenes()
        self.getTotalReads()
        self.getReadsPerGene()
        return

    @staticmethod
    def execute(cmd, worktext="Primerdesign, please wait"):
        """
        Method to execute functions on command line.
        Args:
            cmd: the command line command
            worktext: the printed text to the user

        Returns:
            None: returns to the place of calling
        """
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        Main.printer(worktext)
        jobNr, stderr = p.communicate()
        if p.returncode == 0:
            print "Done! "
        else:
            print "Error"
            for line in stderr:
                print line
            sys.exit(1)
        return

    @staticmethod
    def readRefGenome(referenceDB):
        """
        Reads the whole genome, and saved it for fast access
        Args:
            referenceDB: the name of the genome to which is mapped

        Returns:
            None: returns to the place of calling
        """
        Main.printer("read reference Genome")
        Genome = {}
        referenceFile = file(referenceDB)
        for line in referenceFile:
            line = line.rstrip()
            if ">" in line:
                chromosomName = line.replace(">", "")
                Genome[chromosomName] = list()
            else:
                for seqChar in line:
                    Genome[chromosomName].append(seqChar)
        return Genome

    @staticmethod
    def generatePrimer3Input(outputFile, POI):
        """
        Transfers the POI to a boulder IO file, used by primer3
        Args:
            outputFile: name for saved boulder IO file
            POI: a list with name and sequences of the POI's

        Returns:
            None: returns to the place of calling
        """
        outputFile = file(outputFile, mode='w')
        for key, value in POI.iteritems():
            outputFile.write("PRIMER_SEQUENCE_ID=" + key + "\n")
            outputFile.write("SEQUENCE=" + value + "\n")
            outputFile.write("PRIMER_MIN_SIZE=15\n"
                             "PRIMER_MAX_SIZE=21\n"
                             "PRIMER_NUM_NS_ACCEPTED=0\n"
                             "PRIMER_PRODUCT_SIZE_RANGE=200-500\n"  # change this to change the wanted product size
                             "PRIMER_PRODUCT_OPT_SIZE=300\n"
                             "PRIMER_FILE_FLAG=1\n"
                             "PRIMER_PICK_INTERNAL_OLIGO=0\n"
                             "PRIMER_EXPLAIN_FLAG=1")
        outputFile.close()
        return

    def readsOnGenes(self):
        """
        Makes a subset of the sorted mapped reads to only cover the genes
        Returns:
            None: returns to the place of calling
        """
        newBamFile = self.bamFile.replace(".sorted", "Genes.sorted")
        workLine = "samtools view -hL " + Main.genesFile + " " + self.bamFile + " | " \
                   "samtools view -bS | samtools sort -o " + newBamFile
        PrimerDesign.execute(workLine, "Generating gene reads only")
        self.bamFile = newBamFile
        return

    def getTotalReads(self):
        """
        Gets the total reads in a sorted bam file, and saves it into a totReads file
        Returns:
            None: returns to the place of calling
        """
        workLine = "samtools view -F 0x4 " + self.bamFile + " | cut -f 1 | sort | uniq | wc -l > " \
                                                       + self.bamFile.replace(".sorted", ".totReads")
        PrimerDesign.execute(workLine, "Getting total reads in genes")
        resultFile = file(self.bamFile.replace(".sorted", ".totReads"), mode='r')
        self.totalReads = int(resultFile.readline())
        return

    def getReadsPerGene(self):
        """
        Generate a number for the coverage, normalized against length of the gene and the total reads in the file
        Returns:
            None: returns to the place of calling
        """
        normRedPerLengtDic = defaultdict(dict)
        chrStartNormDic = defaultdict(dict)
        workLine = "coverageBed -abam " + self.bamFile + " -b " + Main.genesFile + " -s > " + self.bamFile.replace(".sorted", ".tsv")
        PrimerDesign.execute(workLine, "Generate reads per gene")
        outputFile = file(self.bamFile.replace(".sorted", ".tsv"), mode='r')
        for line in outputFile:
            line = line.split("\t")
            if line[1] == '-1':
                continue
            chromosom = line[0]
            reads = float(line[4])
            geneLength = float(line[14])
            startPos = int(line[1])
            endPos = int(line[2])
            if line[5] == '-' and chromosom in chrStartNormDic and startPos in chrStartNormDic[chromosom]:
                continue
            normReadsPerLengthPerGene = ((reads/geneLength)/self.totalReads) * math.pow(10, 9)
            if chromosom not in normRedPerLengtDic[normReadsPerLengthPerGene]:
                normRedPerLengtDic[normReadsPerLengthPerGene][chromosom] = list()
            normRedPerLengtDic[normReadsPerLengthPerGene][chromosom].append([startPos, endPos])
            chrStartNormDic[chromosom][startPos] = normReadsPerLengthPerGene
        self.readsPerGene = [normRedPerLengtDic, chrStartNormDic]
        return

    @staticmethod
    def getTreshold(*args):
        """
        Calculates the threshold used for selecting the POI. the threshold is calculated by taking the average of
        the coverage of a gene between n number of files, tkae the average of all thouse averages, and measure the
        standard deviation of the same list of averages
        Args:
            *args: multiple dictionaries orderded by chromosom, start position of gen, coverage of the gene
        Returns:
            treshold: float which is the average - standard deviation
        """
        nrFiles = len(args)
        chrom = set()
        startPos = defaultdict(set)
        sumRPKM = list()
        for valueList in args:
            chrom |= set(valueList.keys())
        chrom = sorted(chrom)
        for valueList in args:
            for key in chrom:
                if key in valueList:
                    startPos[key] |= set(valueList[key].keys())
        for chromkey, value in startPos.iteritems():
            for startPosKey in value:
                sum = 0
                counter = 0
                for i in xrange(nrFiles):
                    if chromkey in args[i] and startPosKey in args[i][chromkey]:
                        counter += 1
                        sum += int(args[i][chromkey][startPosKey])
                sumRPKM.append(float(sum)/float(counter))
        average = numpy.mean(sumRPKM)
        stdDev = numpy.std(sumRPKM)
        treshold = average - stdDev
        return treshold

    @staticmethod
    def getPOI(treshold, *args):
        """
        Generates a dictionary containing the start and end positions of the genes that have a higher or equal coverage
        rate compared to the threshold.
        Args:
            treshold: the set threshold generated by getThreshold
            *args: a dictionary with information coverage rate, chromosom, start - end position of gen

        Returns:
            perChrom: POI dictionary divided on chromosom, used to generate a fasta file
        """
        perChrom = defaultdict(list)
        for dictPerFile in args:
            for key, valueDict in dictPerFile.iteritems():
                if key >= treshold:
                    for chr, listStartEnd in valueDict.iteritems():
                        start = set()
                        end = set()
                        for startEnd in listStartEnd:
                            start.add(startEnd[0])
                            end.add(startEnd[1])
                        perChrom[chr] = [start, end]
        return perChrom
