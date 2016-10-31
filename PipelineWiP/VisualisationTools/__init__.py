#!/usr/bin/env python

from __future__ import division
import threading
from Main import Main

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class VisualisationTools(threading.Thread):
    """
    Currently, this module is not yet implemented in the whole system due to testing
    Used to generate plots and statistics regarding the outcomes of the other modules
    """

    inputDir = ""
    outputDir = ""
    depthOccurance = {}
    depthPerPos = {}
    allList = []
    totalSize = 0

    def __init__(self, inputDir):
        Main.logger.debug("VisualisationTools: i." + inputDir)
        """
        Method for initiating the visualisation tool.
        threading.Thread is called for threaded use of this class
        :param inputDir: the location of the output of divers systems (default: workDir)
        :param outputDir: the location to save the output of these methods to (default: resultDir)
        """
        threading.Thread.__init__(self)
        self.inputDir = inputDir
        return

    def run(self):
        Main.logger.debug("VisualisationTools run:")
        """
        Method called upon by start(), specific for threated runs.
        Please call upon this after finishing the mapping and transfer in ReadAligner
        :returns: None (Null): returns to the place of calling
        """
        self.readBedToLocal()
        self.readCovBedtoLocal()
        oneDepth = self.getCoveragePercentage()
        fiveDepth = self.getCoveragePercentage(5)
        tenDepth = self.getCoveragePercentage(10)
        fifteenDepth = self.getCoveragePercentage(15)
        twentyDepth = self.getCoveragePercentage(20)
        outputFile = file(self.inputDir+".CovRes", mode='w')
        outputFile.write("1 depth:\t" + str(oneDepth) +
                         "\n5 depth:\t" + str(fiveDepth) +
                         "\n10 depth:\t" + str(tenDepth) +
                         "\n15 depth:\t" + str(fifteenDepth) +
                         "\n20 depth:\t" + str(twentyDepth))
        return

    def readBedToLocal(self):
        Main.logger.debug("VisualisationTools readBedToLocal:")
        """
        Reads the BED file into the memory for fast access.
        The BED file here is a BED file giving read depth per position
        Returns:
        """
        bedFile = file(self.inputDir+".bed", mode='r')
        bedList = {}
        scaffoldList = []
        scaffold = ""
        self.allList = []
        for line in bedFile:
            line = line.split('\t')
            if scaffold != line[0]:
                if scaffold != "":
                    bedList[scaffold] = scaffoldList
                scaffold = line[0]
                scaffoldList = []
            scaffoldList.append(line[2])
            self.allList.append(line[2])
        bedList[scaffold] = scaffoldList
        self.depthPerPos = bedList
        return

    def readCovBedtoLocal(self):
        Main.logger.debug("VisualisationTools readCovBedtoLocal:")
        """
        Reads the occurance of a read depth over the whole genome.
        Returns:
        None (Null): returns to the place of calling
        """
        covFile = file(self.inputDir + ".CovBed", mode='r')
        self.depthOccurance = {}
        for line in covFile:
            line = line.split('\t')
            if "genome" in line[0]:
                if line[1] in self.depthOccurance:
                    self.depthOccurance[line[1]] += int(line[2])
                else:
                    self.depthOccurance[line[1]] = int(line[2])
        self.totalSize = int(line[3])
        return

    def getCoveragePercentage(self, depth=1):
        Main.logger.debug("VisualisationTools getCoveragePercentage: d." + str(depth))
        """
        Returns the percentage of coverage over the whole genome by an set read depth
        Input:
        depth (int): depth that you want to get the coverage rate from
        Returns:
        percentageMore (float): the percentage of coverage at the given depth
        """
        lessVal = 0
        for key, value in self.depthOccurance.iteritems():
            key = int(key)
            if key < depth:
                lessVal += value
        percentageLess = lessVal/self.totalSize
        percentageMore = 1.0 - percentageLess
        return percentageMore

class Mapping:

    def __init__(self, scaffold):
        Main.logger.debug("Mapping: s." + scaffold)
        self.hitList = list()
        self.starts = list()
        self.ends = list()
        self.startEnd = list()
        self.scaffold = scaffold
        return

    def hit(self, start, end, ContigID):
        Main.logger.debug("Mapping hit: s." + str(start) + " e." + str(end) + " c." + str(ContigID))
        hitItem = [start, end, ContigID]
        self.startEnd.append([start, end])
        self.starts.append(start)
        self.ends.append(end)
        self.hitList.append(hitItem)
        return

    def getStartEnds(self, startEnd=0):
        Main.logger.debug("Mapping getStartEnds: s." + str(startEnd))
        returnList = list()
        for hitItem in self.hitList:
            returnList.append(hitItem[startEnd])
        return returnList


class scaffoldClass:
    colorList = list()
    isolates = list()
    ylabel = ""
    xlabel = ""

    def __init__(self, scaffold):
        self.scaffold = scaffold
        self.hitDict = dict()
        return

    def hit(self, isolat, start, end):
        if isolat not in self.hitDict:
            self.hitDict[isolat] = list()
        self.hitDict[isolat].append([start, end])
        return

    def makePlot(self, number):
        Main.logger.debug(self.scaffold)
        plt.ioff()
        plt.figure(number)
        for isolat, value in iter(self.hitDict.items()):
            i = scaffoldClass.isolats.index(isolat)
            for startEnd in value:
                line = plt.plot([startEnd[0], startEnd[1]], [i + 1, i + 1])
                plt.setp(line, color=scaffoldClass.colorList[i])
        plt.ylim(0, len(scaffoldClass.isolates)+1)
        plt.ylabel(scaffoldClass.ylabel)
        plt.xlabel(scaffoldClass.xlabel)
        plt.grid(False)
        plt.title(str(self.scaffold))
        place = Main.resultDir + "/" + str(self.scaffold) + ".svg"
        plt.savefig(place, format='svg', dpi=1200)
        plt.clf()
        return
