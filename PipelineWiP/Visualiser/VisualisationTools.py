#!/usr/bin/env python

from __future__ import division
import threading

class VisualisationTools(threading.Thread):

    inputDir = ""
    outputDir = ""
    depthOccurance = {}
    depthPerPos = {}
    allList = []
    totalSize = 0

    def __init__(self, inputDir, outputDir):
        threading.Thread.__init__(self)
        self.inputDir = inputDir
        self.outputDir = outputDir
        return

    def run(self):
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
        lessVal = 0
        for key, value in self.depthOccurance.iteritems():
            key = int(key)
            if key < depth:
                lessVal += value
        print self.totalSize
        percentageLess = lessVal/self.totalSize
        percentageMore = 1.0 - percentageLess
        return percentageMore
