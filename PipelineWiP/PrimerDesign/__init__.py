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
from NUCmer import NUCmer
import itertools
import VisualisationTools

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
                line = line.replace(">", "")
                line = line.split(" ")
                chromosomName = line[0]
                Genome[chromosomName] = list()
            else:
                for seqChar in line:
                    Genome[chromosomName].append(seqChar)
        return Genome

    @staticmethod
    def readGFF(gffFile, genome):
        POI = dict()
        gffFile = open(gffFile, mode='r')
        for line in gffFile:
            line = line.split("\t")
            scaffold = line[0]
            start = int(line[3]) - 1
            end = int(line[4]) - 1
            id = line[8]
            id = id.replace("ID=", "")
            id = id.strip()
            sequence = ""
            while start < end:
                sequence += genome[scaffold][start]
                start += 1
            POI[id] = sequence
        return POI

    @staticmethod
    def saveFasta(outputFile, POI):
        outputFile = open(Main.workDir + "/" + outputFile, mode="w")
        for key, value in POI.iteritems():
            outputFile.write(">" + str(key) + "\n" + value + "\n")
        outputFile.close()

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
            outputFile.write("SEQUENCE_ID=POI" + key + "\n")
            outputFile.write("SEQUENCE_TEMPLATE=" + value + "\n")
            outputFile.write("PRIMER_MIN_SIZE=15\n"
                             "PRIMER_MAX_SIZE=21\n"
                             "PRIMER_PICK_LEFT_PRIMER=1\n"
                             "PRIMER_PICK_RIGHT_PRIMER=1\n"
                             "PRIMER_NUM_NS_ACCEPTED=0\n"
                             "PRIMER_PRODUCT_SIZE_RANGE=200-500\n"  # change this to change the wanted product size
                             "PRIMER_PRODUCT_OPT_SIZE=300\n"
                             "P3_FILE_FLAG=1\n"
                             "PRIMER_PICK_INTERNAL_OLIGO=0\n"
                             "PRIMER_EXPLAIN_FLAG=1\n"
                             "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/mnt/apps/primer3-2.3.0/src/primer3_config/\n"
                             "=\n")
        outputFile.close()
        workLine = "/mnt/apps/primer3-2.3.0/src/primer3_core -format_output -output=PrimerList.txt " + outputFile.name
        PrimerDesign.execute(workLine, "Running primer3")
        return

class PrimerDesignByDenovo:

    def __init__(self):
        self.coordsInfo = list()
        self.geneInfo = dict()
        self.hitList = dict()
        return

    def readCoords(self, coordsFile):
        coordsHits = dict()
        coordsFile = open(coordsFile+".coords", mode="r")
        for line in coordsFile:
            line = line.split()
            if len(line) == 16:
                if line[14] not in coordsHits:
                    coordsHits[line[14]] = NUCmer(line[14])
                coordsHits[line[14]].hit(int(line[0]), int(line[1]), line[15])
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
        i = 0
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
                    outputFile.write(scaffold+delimiter+"Dundee"+delimiter+"CDS"+delimiter+str(item[0])+
                                     delimiter+str(item[1])+delimiter+"."+delimiter+"+"+delimiter+"0"+delimiter+
                                     "ID="+str(i)+endLine)
                    i += 1
                    starts.add(item[0])
        outputFile.close()
        return


class PrimerDesignByMapping:

    y = 0

    def __init__(self):
        self.coordsInfo = list()
        self.HitList = dict()
        self.geneInfo = dict()
        return

    def generateCoords(self, depthPerPos, depthLimit=12, counterLimit=28):
        print depthLimit
        print counterLimit
        print len(depthPerPos)
        coordsHits = dict()
        allscaffold = list(depthPerPos.keys())
        for scaffold in allscaffold:
            counter = 0
            pos = 0
            for depth in depthPerPos[scaffold]:
                pos += 1
                if int(depth) >= depthLimit:
                    counter += 1
                elif counter >= counterLimit:
                    if scaffold not in coordsHits:
                        coordsHits[scaffold] = VisualisationTools.Mapping(scaffold)
                    PrimerDesignByMapping.y += 1
                    coordsHits[scaffold].hit(pos-counter, pos-1, PrimerDesignByMapping.y)
                    counter = 0
                else:
                    counter = 0
        self.coordsInfo.append(coordsHits)
        return

    def readGenes(self, geneFile=Main.gffFile):
        geneFile = open(geneFile, mode="r")
        for line in geneFile:
            if "exon" in line:
                line = line.split()
                if line[0] not in self.geneInfo:
                    self.geneInfo[line[0]] = VisualisationTools.Mapping(line[0])
                self.geneInfo[line[0]].hit(int(line[3]), int(line[4]), line[8])
        geneFile.close()
        return

    ## TODO: optimise this part. why does it take longer now? why does it not skip if i == y?
    def generateHitList(self):
        self.HitList = dict()
        allScaffold = set()
        i = 0
        for coordsDict in self.coordsInfo:
            allScaffold.update(set(coordsDict.keys()))
        for scaffold in allScaffold:
            HitList = list()
            if all(scaffold in key for key in self.coordsInfo):
                startList = list()
                endList = list()
                for coordsDict in self.coordsInfo:
                    startList.append(coordsDict[scaffold].starts)
                    endList.append(coordsDict[scaffold].ends)
                for starts in startList:
                    for start in starts:
                        trueFalseCheck = list()
                        i += 1
                        y = 0
                        end = endList[startList.index(starts)][starts.index(start)]
                        for compareStarts in startList:
                            y += 1
                            if i == y:
                                continue
                            checkStarts = [i for i in compareStarts if i <= start]
                            checkEnds = [i for i in endList[startList.index(compareStarts)] if i > start]
                            indexSetStart = set()
                            indexSetEnd = set()
                            for item in checkStarts:
                                indexSetStart.add(compareStarts.index(item))
                            for item in checkEnds:
                                indexSetEnd.add(endList[startList.index(compareStarts)].index(item))
                            indexSet = indexSetStart.intersection(indexSetEnd)
                            if len(indexSet) > 0:
                                trueFalseCheck.append(True)
                            else:
                                trueFalseCheck.append(False)
                            for index in indexSet:
                                if end > endList[startList.index(compareStarts)][index]:
                                    end = endList[startList.index(compareStarts)][index]
                        if all(check for check in trueFalseCheck):
                            HitList.append([start, end])
            HitList.sort()
            HitList = list(HitList for HitList,_ in itertools.groupby(HitList))
            print len(HitList)
            self.HitList[scaffold] = HitList

    def generateGeneSpecificHitList(self):
        i = 0
        outputFile = open(Main.workDir + "/MapperPoI.gff", mode="w")
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
                if len(GindexSet) > 0:
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
            starts = set()
            for item in Hitlist:
                if item[0] not in starts:
                    outputFile.write(scaffold + delimiter + "Dundee" + delimiter + "CDS" + delimiter + str(item[0]) +
                                     delimiter + str(item[1]) + delimiter + "." + delimiter + "+" + delimiter + "0" +
                                     delimiter + "ID=" + str(i) + endLine)
                    i += 1
                    starts.add(item[0])
        outputFile.close()
        return