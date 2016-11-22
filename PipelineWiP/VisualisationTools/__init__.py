#!/usr/bin/env python
"""
Created on: 29-08-2016
@author: H.J.C. Cornelisse
Class is used to analyse the mapped reads.
Dependencies:
- matplotlib (not tested)
"""
from __future__ import division
import threading
from Main import Main

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class VisualisationTools(threading.Thread):
    """
    class that shows the different coverage rates according to the results of the read alignment
    """

    inputDir = ""
    outputDir = ""
    depthOccurance = {}
    depthPerPos = {}
    allList = []
    totalSize = 0

    def __init__(self, inputDir):
        """
        Method for initiating the visualisation tool.
        threading.Thread is called for threaded use of this class
        :param inputDir: the location of the output of divers systems (default: workDir)
        """
        Main.logger.debug("VisualisationTools: i." + inputDir)
        threading.Thread.__init__(self)
        self.inputDir = inputDir
        return

    def run(self):
        """
        Method called upon by start(), specific for threated runs.
        Please call upon this after finishing the mapping and transfer in ReadAligner
        :return:
        """
        Main.logger.debug("VisualisationTools run:")
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
        """
        Reads the BED file into the memory for fast access.
        The BED file here is a BED file giving read depth per position
        :return:
        """
        Main.logger.debug("VisualisationTools readBedToLocal:")
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
        """
        Reads the occurance of a read depth over the whole genome.
        :return:
        """
        Main.logger.debug("VisualisationTools readCovBedtoLocal:")
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
        """
        Returns the percentage of coverage over the whole genome by an set read depth
        :param depth: depth that you want to get the coverage rate from
        :return: percentageMore (float): the percentage of coverage at the given depth
        """
        Main.logger.debug("VisualisationTools getCoveragePercentage: d." + str(depth))
        lessVal = 0
        for key, value in self.depthOccurance.iteritems():
            key = int(key)
            if key < depth:
                lessVal += value
        percentageLess = lessVal/self.totalSize
        percentageMore = 1.0 - percentageLess
        return percentageMore

class Mapping:
    """
    Saving information to work with the primer design in combination with the direct mapping
    """

    def __init__(self, scaffold):
        """
        start with scaffold
        :param scaffold:
        """
        Main.logger.debug("Mapping: s." + scaffold)
        self.hitList = list()
        self.starts = list()
        self.ends = list()
        self.startEnd = list()
        self.scaffold = scaffold
        return

    def hit(self, start, end, ContigID):
        """
        Information that saved the start, end and the origin of the coverage
        :param start: start position on the scaffold
        :param end: end position on the scaffold
        :param ContigID: name of the contig
        :return:
        """
        Main.logger.debug("Mapping hit: s." + str(start) + " e." + str(end) + " c." + str(ContigID))
        hitItem = [start, end, ContigID]
        self.startEnd.append([start, end])
        self.starts.append(start)
        self.ends.append(end)
        self.hitList.append(hitItem)
        return

    def getStartEnds(self, startEnd=0):
        """
        returns a list of all the starts or ends of the scaffold
        :param startEnd: 0 for starts, 1 for ends
        :return: list of all the starts or ends in the scaffold
        """
        Main.logger.debug("Mapping getStartEnds: s." + str(startEnd))
        returnList = list()
        for hitItem in self.hitList:
            returnList.append(hitItem[startEnd])
        return returnList


class scaffoldClass:
    """
    Not yet in use
    """
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
