#!/usr/bin/env python

"""
Created on: 29-08-2016
@author: H.J.C. Cornelisse
Class is used to generate points of interessed from the mapped reads.
Dependencies:
- primer3 2.3.0 (tested)
- bedtools 2.23.0 (tested)
"""

from Main import Main
from NUCmer import NUCmer
import VisualisationTools


class PrimerDesign:
    """
    overall primer design. all function in this module are directly accessible.
    """
    def __init__(self):
        Main.logger.debug("PrimerDesign:")
        return

    @staticmethod
    def readRefGenome(referenceDB):
        """
        Reads the whole genome, and saved it for fast access
        :param referenceDB: the name of the genome to which is mapped
        :return: returns to the place of calling
        """
        Main.logger.debug("readRefGenome: r." + referenceDB)
        Main.logger.info("read reference Genome")
        Genome = {}
        referenceFile = file(referenceDB)
        chromosomName = ""
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
    def removeDuplicate(gffFile):
        """
        removes duplicates in the gff file (from the regions subsets
        :param gffFile: location to the gff file
        :return:
        """
        Main.logger.debug("removeDuplicate: g." + gffFile)
        gffFile = open(gffFile, mode='rw+')
        items = list()
        startPrev = 0
        endPrev = 0
        scaffoldPrev = ""
        for line in gffFile:
            line = line.split("\t")
            scaffold = line[0]
            start = int(line[3])
            end = int(line[4])
            if start != startPrev and end != endPrev and scaffold == scaffoldPrev:
                items.append([scaffold, start, end])
            elif scaffold != scaffoldPrev:
                items.append([scaffold, start, end])
            scaffoldPrev = scaffold
            startPrev = start
            endPrev = end
        gffFile.seek(0)
        for i, line in enumerate(items):
            gffFile.write(line[0] + "\t" + "." + "\t" + "CDS" + "\t" + str(line[1]) + "\t" + str(line[2]) +
                          "\t" + "." + "\t" + "+" + "\t" + "0" + "\t" + "ID=" + str(i) + "\n")
        gffFile.truncate()
        gffFile.close()

    @staticmethod
    def readGFF(gffFile, genome):
        """
        reads the gff file to the start and the end of the region, and saves the sequence to POI
        :param gffFile: location to gff subset
        :param genome: genome loaded in memory
        :return: POI dictorionary with id and sequence as value
        """
        Main.logger.debug("readGFF: gff." + gffFile + " g." + str(len(genome)))
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
        """
        saves the POI dictionary to a fasta file
        :param outputFile: location to save fasta file to
        :param POI: dictionary with id and sequence
        :return: -
        """
        Main.logger.debug("saveFasta: o." + outputFile + " p." + str(len(POI)))
        outputFile = open(outputFile, mode="w")
        for key, value in POI.iteritems():
            outputFile.write(">" + str(key) + "\n" + value + "\n")
        outputFile.close()

    @staticmethod
    def generatePrimer3Input(outputFile, POI):
        """
        generates the primers depending on the POI (id and sequence)
        :param outputFile: position to save the primer design settings and inputs
        :param POI: dictionary with ID and sequence
        :return:
        """
        Main.logger.debug("generatePrimer3Input: o." + outputFile + " p." + str(len(POI)))
        outputFilePart = outputFile.split(".")
        outputDir = outputFilePart[0]
        Main.makeDirectory(outputDir)
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
                             "P3_FILE_FLAG=0\n"
                             "PRIMER_PICK_INTERNAL_OLIGO=0\n"
                             "PRIMER_EXPLAIN_FLAG=1\n"
                             "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/mnt/apps/primer3-2.3.0/src/primer3_config/\n"
                             "=\n")
        outputFile.close()
        workLine = "/mnt/apps/primer3-2.3.0/src/primer3_core -format_output -output=" + outputDir +\
                   "/PrimerList.txt " + outputFile.name
        Main.execute(workLine, "Running primer3")
        return

    @staticmethod
    def runIntersect(coordsFile, outputName):
        """
        runs bedtools intersect to find the intersection between read information
        :param coordsFile: file with the coordinations of the alignment
        :param outputName: location of the gff output file
        :return: gff file name
        """
        Main.logger.debug("runIntersect: c." + str(len(coordsFile)) + " o." + outputName)
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
                    Main.execute(workLine, "Running intersect, please wait")
                workLine = "intersectBed -a " + coordsFile
            coordsFile = coordsFileList
        fromBED = open(coordsFile[0], mode="r")
        toGFFFileName = Main.workDir + outputName
        toGFF = open(toGFFFileName, mode="w")
        for i, line in enumerate(fromBED):
            line = line.split("\t")
            toGFF.write(line[0] + "\t" + "." + "\t" + "CDS" + "\t" + line[1] + "\t" + line[
                2].rstrip() + "\t" + "." + "\t" + "+" + "\t" + "0" + "\t" +
                        "ID=" + str(i) + "\n")
        return toGFFFileName

    @staticmethod
    def runMethodIntersect(coordsFile, outputName):
        """
        Tool to find the intersection between the tools
        :param coordsFile: files with the coordinations
        :param outputName: location of the intersection file
        :return: location to the intersection file
        """
        Main.logger.debug("runMethodIntersect: c." + str(len(coordsFile)) + " o." + outputName)
        workLine = "intersectBed -a " + coordsFile[0] + " -b " + coordsFile[1] + " > " + outputName
        Main.execute(workLine, "Running method intersect, please wait")
        return outputName

    @staticmethod
    def runMethodSubstract(coordsFile, outputName):
        """
        removes the intersection from the subsets, and makes a new subset with only unique regions for that method
        :param coordsFile: coordsfile like mapper and intersection
        :param outputName: outputfile like mapperunique
        :return: file location
        """
        Main.logger.debug("runMethodSubstract: c." + str(len(coordsFile)) + " o." + outputName)
        workLine = "subtractBed -a " + coordsFile[0] + " -b " + coordsFile[1] + " > " + outputName
        Main.execute(workLine, "Running substraction, please wait")
        return outputName

    @staticmethod
    def readGenes(y, geneFile=Main.gffFile):
        """
        reads the gene file, and saves the gene information as bed file
        :param y: number send by caller
        :param geneFile: location of the gene file
        :return: location of the output file
        """
        Main.logger.debug("readGenes: y." + str(y) + " g." + geneFile)
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
    """
    all steps used for the denovo method
    """

    def __init__(self):
        """
        initiation step for designing regions for the denovo method
        """
        Main.logger.debug("PrimerDesignByDenovo:")
        self.y = 0
        self.coordsFile = list()
        if Main.gffFile is not None:
            self.coordsFile.append(PrimerDesign.readGenes("D" + str(self.y), Main.gffFile))
            self.y += 1
        self.geneInfo = dict()
        self.hitList = dict()
        return

    def readCoords(self, coordsFile):
        """
        reads the coordinates from the results of nucmer
        :param coordsFile: location to the nucmer file
        :return:
        """
        Main.logger.debug("readCoords: c." + coordsFile)
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
            mapper.combineStartEnds()
            for starts, ends in zip(mapper.starts, mapper.ends):
                bedFile.write(scaffold + "\t" + str(starts) + "\t" + str(ends) + "\n")
        bedFile.close()
        self.coordsFile.append(outputFile)
        self.y += 1
        return


class PrimerDesignByMapping:
    """
    all steps used for the mapping method
    """

    y = 0

    def __init__(self):
        """
        initiation step for designing regions for the mapping method
        """
        Main.logger.debug("PrimerDesignByMapping:")
        self.HitList = dict()
        self.coordsFile = list()
        self.y = 0
        if Main.gffFile is not None:
            self.coordsFile.append(PrimerDesign.readGenes("M"+str(self.y), Main.gffFile))
            self.y += 1
        return

    def generateCoords(self, depthPerPos, depthLimit=12, counterLimit=28):
        """
        generates coord depending on the settings. the coverage is received from bowtie
        :param depthPerPos: information of the depth per position
        :param depthLimit: miminum depth that needs to be reached
        :param counterLimit: miminum window size
        :return:
        """
        Main.logger.debug("PrimerDesignByMapping generateCoords: df." + str(len(depthPerPos)) +
                          " d." + str(depthLimit) + " c." + str(counterLimit))
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


