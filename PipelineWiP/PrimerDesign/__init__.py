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
from NUCmer import NUCmer
import itertools

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


class PrimerDesignByDenovo:

    def __init__(self):
        self.coordsInfo = list()
        self.geneInfo = dict()
        self.hitList = dict()
        return

    def readCoords(self, coordsFile):
        coordsHits = dict()
        coordsFile = open(coordsFile, mode="r")
        for line in coordsFile:
            line = line.split()
            if len(line) == 19:
                if line[17] not in coordsHits:
                    coordsHits[line[17]] = NUCmer(line[17])
                coordsHits[line[17]].hit(int(line[0]), int(line[1]), line[18])
        self.coordsInfo.append(coordsHits)
        coordsFile.close()
        return

    def readGenes(self, geneFile=Main.gffFile):
        geneFile = open(geneFile, mode="r")
        for line in geneFile:
            if "exon" in line:
                line = line.split()
                if line[0] not in self.geneInfo:
                    self.geneInfo[line[0]] = NUCmer(line[0])
                self.geneInfo[line[0]].hit(int(line[3]), int(line[4]), line[8])
        geneFile.close()
        return

    def generateHitList(self):
        self.HitList = dict()
        allScaffold = set()
        i = 0
        for coordsDict in self.coordsInfo:
            allScaffold.update(set(coordsDict.keys()))
        for scaffold in allScaffold:
            HitList = list()
            if all(scaffold in key for key in self.coordsInfo):
                StartList = list()
                EndList = list()
                for coordsDict in self.coordsInfo:
                    coordsDict[scaffold].combineStartEnds()
                    StartList.append(coordsDict[scaffold].starts)
                    EndList.append(coordsDict[scaffold].ends)
                for starts in StartList:
                    for start in starts:
                        trueFalseCheck = list()
                        i += 1
                        y = 0
                        end = EndList[StartList.index(starts)][starts.index(start)]
                        for compareStarts in StartList:
                            y += 1
                            if i == y:
                                continue
                            checkStarts = [i for i in compareStarts if i <= start]
                            checkEnds = [i for i in EndList[StartList.index(compareStarts)] if i > start]
                            indexSetStart = set()
                            indexSetEnd = set()
                            for item in checkStarts:
                                indexSetStart.add(compareStarts.index(item))
                            for item in checkEnds:
                                indexSetEnd.add(EndList[StartList.index(compareStarts)].index(item))
                            indexSet = indexSetStart.intersection(indexSetEnd)
                            if len(indexSet) > 0:
                                trueFalseCheck.append(True)
                            else:
                                trueFalseCheck.append(False)
                            for index in indexSet:
                                if end > EndList[StartList.index(compareStarts)][index]:
                                    end = EndList[StartList.index(compareStarts)][index]
                        if all(check for check in trueFalseCheck):
                            HitList.append([start, end])
            HitList.sort()
            HitList = list(HitList for HitList,_ in itertools.groupby(HitList))
            print scaffold + " " + str(len(HitList))
            self.HitList[scaffold] = HitList

    def generateGeneSpecificHitList(self):
        outputFile = open(Main.workDir+"/denovoPoI.gff", mode="w")
        delimiter = "\t"
        endLine = "\n"
        intersectScaffold = set.intersection(set(self.HitList.keys()), set(self.geneInfo.keys()))
        for scaffold in intersectScaffold:
            Hitlist = list()
            self.geneInfo[scaffold].combineStartEnds()
            geneStart = self.geneInfo[scaffold].starts
            geneEnd = self.geneInfo[scaffold].ends
            CoverStart = list()
            CoverEnd = list()
            for item in self.HitList[scaffold]:
                start = item[0]
                end = item[1]
                GeneHit = False
                CoverStart.append(start)
                CoverEnd.append(end)
                Gstart = [i for i in geneStart if i <= start]
                Gend = [i for i in geneEnd if i > start]
                GindexSetStart = set()
                GindexSetEnd = set()
                for Gitem in Gstart:
                    GindexSetStart.add(geneStart.index(Gitem))
                for Gitem in Gend:
                    GindexSetEnd.add(geneEnd.index(Gitem))
                GindexSet = GindexSetStart.intersection(GindexSetEnd)
                if len(GindexSet) > 0 :
                    GeneHit = True
                for index in GindexSet:
                    if end > geneEnd[index]:
                        end = geneEnd[index]
                if GeneHit is True:
                    Hitlist.append([start, end])

            for Gstart in geneStart:
                start = Gstart
                end = geneEnd[geneStart.index(Gstart)]
                CoverHit = False
                Cstart = [i for i in CoverStart if i <= start]
                Cend = [i for i in CoverEnd if i > start]
                CindexSetStart = set()
                CindexSetEnd = set()
                for Citem in Cstart:
                    CindexSetStart.add(CoverStart.index(Citem))
                for Citem in Cend:
                    CindexSetEnd.add(CoverEnd.index(Citem))
                CindexSet = CindexSetStart.intersection(CindexSetEnd)
                if len(CindexSet) > 0:
                    CoverHit = True
                for index in CindexSet:
                    if end > CoverEnd[index]:
                        end = CoverEnd[index]
                if CoverHit is True:
                    Hitlist.append([start, end])
            Hitlist.sort()
            Hitlist = list(Hitlist for Hitlist, _ in itertools.groupby(Hitlist))
            print scaffold + " " + str(len(Hitlist))
            starts = set()
            for item in Hitlist:
                if item[0] not in starts:
                    outputFile.write(scaffold+delimiter+"Dundee"+delimiter+
                                     "CDS"+delimiter+item[0]+delimiter+item[1])
                    print str(item[0]) + " : " + str(item[1])
                    starts.add(item[0])
        return

