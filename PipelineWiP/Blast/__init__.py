#!/usr/bin/env python

import threading
from Main import Main
import copy

class Blast(threading.Thread):
    y = 0

    def __init__(self, query, genomeList):
        Main.logger.debug("Blast: q." + query + " g." + str(len(genomeList)))
        Blast.y += 1
        self.y = copy.deepcopy(Blast.y)
        threading.Thread.__init__(self)
        self.genomeList = genomeList
        self.query = query
        return

    def run(self):
        Main.logger.debug("Blast run:")
        allAlias = ""
        for genome in self.genomeList:
            self.makeDatabase(str(genome), Main.workDir)
            allAlias += Main.workDir + "/" + genome + " "
        allAlias = allAlias.rstrip()
        self.aliasTool(allAlias, Main.workDir)
        self.doBlast(Main.workDir)
        self.interpertBlast(Main.workDir)
        self.removeHits()
        return

    def makeDatabase(self, genomeFile, outputDir):
        Main.logger.debug("Blast makeDatabase: g." + genomeFile + " o." + outputDir)
        workline = "makeblastdb -in " + Main.genomeAdd + genomeFile + " -parse_seqids -dbtype nucl -out " + outputDir +\
                   "/" + genomeFile
        Main.execute(workline)
        return

    def aliasTool(self, allAlias, outputDir=Main.workDir):
        Main.logger.debug("Blast aliasTool: a." + allAlias + " o." + outputDir)
        output = outputDir+"/others"
        workline = "blastdb_aliastool -dblist \"" + allAlias + "\" -dbtype nucl -out " + output + " -title others"
        Main.execute(workline, "Making one database")
        return

    def doBlast(self, outputDir=Main.workDir):
        Main.logger.debug("Blast doBlast: o." + outputDir)
        workline = "blastn -db " + outputDir + "/others -query " + self.query + " -out " + outputDir +\
                   "/blastResult"+str(self.y)+".csv -outfmt \"10 std\""
        Main.execute(workline, "Running blast")
        return

    def interpertBlast(self, workDir=Main.workDir):
        Main.logger.debug("Blast interpertBlast: w." + workDir)
        blastFile = open(workDir+"/blastResult" + str(self.y) + ".csv")
        self.hitSet = set()
        for line in blastFile:
            line = line.split(",")
            try:
                self.hitSet.add(int(line[0]))
            except ValueError:
                Main.logger.error("Something went wrong here: "+line[0])
        return

    def removeHits(self):
        Main.logger.debug("Blast removeHits:")
        gffFile = self.query[:-2]
        inputFile = open(gffFile + "gff", 'r')
        outputFile = open(gffFile + "unique.gff", 'w')
        for line in inputFile:
            lineSplit = line.split("\t")
            id = lineSplit[8]
            id = int(id.replace("ID=", ""))
            if id not in self.hitSet:
                outputFile.write(line)
        inputFile.close()
        outputFile.close()
        return