#!/usr/bin/env python
#.. todo:: test this shit out!
""
from __future__ import division
import threading

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

    def __init__(self, inputDir, outputDir):
        """
        Method for initiating the visualisation tool.
        threading.Thread is called for threaded use of this class
        :param inputDir: the location of the output of divers systems (default: workDir)
        :param outputDir: the location to save the output of these methods to (default: resultDir)
        """
        threading.Thread.__init__(self)
        self.inputDir = inputDir
        self.outputDir = outputDir
        return

    def run(self):
        """
        Method called upon by start(), specific for threated runs.
        Please call upon this after finishing the mapping and transfer in ReadAligner
        :returns: None (Null): returns to the place of calling
        """
        self.readBedToLocal()
        # self.readCovBedtoLocal()
        # oneDepth = self.getCoveragePercentage()
        # fiveDepth = self.getCoveragePercentage(5)
        # tenDepth = self.getCoveragePercentage(10)
        # fifteenDepth = self.getCoveragePercentage(15)
        # twentyDepth = self.getCoveragePercentage(20)
        # outputFile = file(self.inputDir+".CovRes", mode='w')
        # outputFile.write("1 depth:\t" + str(oneDepth) +
        #                  "\n5 depth:\t" + str(fiveDepth) +
        #                  "\n10 depth:\t" + str(tenDepth) +
        #                  "\n15 depth:\t" + str(fifteenDepth) +
        #                  "\n20 depth:\t" + str(twentyDepth))
        return

    def readBedToLocal(self):
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

    def saveAsCSV(self):
        """
        Saves the information from the BED depth per position to a CSV file for visualisation
        important: whole genome is big! think before visualising!
        Returns:
        None (Null): returns to the place of calling
        """
        for key, value in self.depthPerPos.iteritems():
            i = 1
            outputFile = file(self.outputDir + "/" + key + ".csv", mode='w')
            outputFile.write("pos,depth\n")
            for depth in value:
                outputFile.write(str(i)+","+str(depth))
                i += 1
            outputFile.close()
        return

    def getCoveragePercentage(self, depth=1):
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
        print self.totalSize
        percentageLess = lessVal/self.totalSize
        percentageMore = 1.0 - percentageLess
        return percentageMore
