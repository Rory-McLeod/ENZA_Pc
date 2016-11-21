#!/usr/bin/env python
"""
Created on: 29-08-2016
latest update: 01-11-2016
@author: H.J.C. Cornelisse
Module is used to run all parts used for the blasting against other spp
Dependencies:
- Blast 2.2.28+ (tested)
Todo:
- tests
- add comments
"""
import threading
from Main import Main
import copy


class Blast(threading.Thread):
    """
    Blast method takes care of all functions regarding removal of regions that occur in other genomes.
    """
    y = 0

    def __init__(self, query, genomeList):
        """
        Initation step for running blast in threads.
        :param query: location to fasta file
        :param genomeList: list of genome locations
        """
        Main.logger.debug("Blast: q." + query + " g." + str(len(genomeList)))
        Blast.y += 1
        self.y = copy.deepcopy(Blast.y)
        threading.Thread.__init__(self)
        self.genomeList = genomeList
        self.query = query
        return

    def run(self):
        """
        This part is called upon by running start().
        Before this is called upon, a local database should be made and renamed to others
        :return:
        """
        Main.logger.debug("Blast run:")
        self.doBlast(Main.workDir)
        self.interpertBlast(Main.workDir)
        self.removeHits()
        return

    @staticmethod
    def makeDatabase(genomeFile, outputDir):
        """
        Makes a local database from genome file
        :param genomeFile: location to genome file
        :param outputDir: location to store the local database
        :return:
        """
        Main.logger.debug("Blast makeDatabase: g." + genomeFile + " o." + outputDir)
        workline = "makeblastdb -in " + Main.genomeAdd + genomeFile + " -parse_seqids -dbtype nucl -out " + outputDir +\
                   "/" + genomeFile
        Main.execute(workline)
        return

    @staticmethod
    def aliasTool(allAlias, outputDir=Main.workDir):
        """
        Combines all the local databases to one local database
        :param allAlias: names of all the local databases
        :param outputDir: location to save the combined database
        :return:
        """
        Main.logger.debug("Blast aliasTool: a." + allAlias + " o." + outputDir)
        output = outputDir+"/others"
        workline = "blastdb_aliastool -dblist \"" + allAlias + "\" -dbtype nucl -out " + output + " -title others"
        Main.execute(workline, "Making one database")
        return

    def doBlast(self, outputDir=Main.workDir):
        """
        Runs the blast based on the local database and the query.
        Before this, a local database should be made named others
        :param outputDir: location where the blast results are saved
        :return:
        """
        Main.logger.debug("Blast doBlast: o." + outputDir)
        workline = "blastn -db " + outputDir + "/others -query " + self.query + " -out " + outputDir +\
                   "/blastResult"+str(self.y)+".csv -outfmt \"10 std\""
        Main.execute(workline, "Running blast")
        return

    def interpertBlast(self, workDir=Main.workDir):
        """
        Opens the blast results, and finds all POI that are found in other genomes
        :param workDir: location where the results are saved
        :return:
        """
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
        """
        Opens the gff file, and removes all POI that are found in other genomes.
        :return:
        """
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