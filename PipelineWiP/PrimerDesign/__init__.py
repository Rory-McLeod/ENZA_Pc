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
    @staticmethod
    def runIntersect(coordsFile, outputName):
        numb = len(coordsFile)
        while len(coordsFile) > 1:
            coordsFileList = list()
            workLine = ""
            workList = list()
            for i, coordsFile in enumerate(coordsFile):
                if i != 0:
                    outPutDir = Main.workDir + "/" + str(numb) + ".bed"
                    numb += 1
                    workList.append(outPutDir)
                    workLine += " -b " + coordsFile + " > " + outPutDir
                    coordsFileList.append(outPutDir)
                    PrimerDesign.execute(workLine, "Running intersect, please wait")
                workLine = "intersectBed -a " + coordsFile
            coordsFile = coordsFileList
        fromBED = open(coordsFile[0], mode="r")
        toGFFFileName = Main.workDir + outputName
        toGFF = open(toGFFFileName, mode="w")
        for i, line in enumerate(fromBED):
            line = line.split("\t")
            toGFF.write(line[0] + "\t" + "." + "\t" + "CDS" + "\t" + line[1] + "\t" + line[
                2].rstrip() + "\t" + "." + "\t" + "+" + "\t" + "0" + "\t"
                        + "ID=" + str(i) + "\n")
        return toGFFFileName

    @staticmethod
    def readGenes(y, geneFile=Main.gffFile):
        print geneFile
        geneInfo = dict()
        geneFile = open(geneFile, mode="r")
        for line in geneFile:
            if "exon" in line:
                line = line.split()
                if line[0] not in geneInfo:
                    geneInfo[line[0]] = VisualisationTools.Mapping(line[0])
                geneInfo[line[0]].hit(int(line[3]), int(line[4]), line[8])
        geneFile.close()
        outputFile = Main.workDir + "/" + y + ".bed"
        bedFile = open(outputFile, mode="w")
        for scaffold, mapper in geneInfo.iteritems():
            for starts, ends in zip(mapper.starts, mapper.ends):
                bedFile.write(scaffold + "\t" + str(starts) + "\t" + str(ends) + "\n")
        bedFile.close()
        return outputFile


class PrimerDesignByDenovo:

    def __init__(self):
        self.y = 0
        self.coordsFile = list()
        self.coordsFile.append(PrimerDesign.readGenes("D" + str(self.y), Main.gffFile))
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
        coordsFile.close()
        outputFile = Main.workDir + "/D" + str(self.y) + ".bed"
        bedFile = open(outputFile, mode="w")
        for scaffold, mapper in coordsHits.iteritems():
            for starts, ends in zip(mapper.starts, mapper.ends):
                bedFile.write(scaffold + "\t" + str(starts) + "\t" + str(ends) + "\n")
        bedFile.close()
        self.coordsFile.append(outputFile)
        return


class PrimerDesignByMapping:

    def __init__(self):
        self.HitList = dict()
        self.coordsFile = list()
        self.y = 0
        self.coordsFile.append(PrimerDesign.readGenes("M"+str(self.y), Main.gffFile))
        return

    def generateCoords(self, depthPerPos, depthLimit=12, counterLimit=28):
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
        outputFile = Main.workDir + "/M" + str(self.y) + ".bed"
        bedFile = open(outputFile, mode="w")
        for scaffold, mapper in coordsHits.iteritems():
            for starts, ends in zip(mapper.starts, mapper.ends):
                bedFile.write(scaffold + "\t" + str(starts) + "\t" + str(ends) + "\n")
        bedFile.close()
        self.coordsFile.append(outputFile)
        self.y += 1
        return


