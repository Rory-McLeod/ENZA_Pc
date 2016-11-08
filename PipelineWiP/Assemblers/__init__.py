#!/usr/bin/env python
"""
Created on: 29-08-2016
latest update: 01-11-2016
@author: H.J.C. Cornelisse
Module is used to run all parts used for the de novo assembly
Dependencies:
- QUAST 2.3 (tested)
- SPAdes 3.1.1. (tested)
Todo:
- tests
- add comments
"""
import threading
from Main import Main
import os


class Assemblers(threading.Thread):
    # TODO: change gffChanger for output in Main.workdir instead of program place
    """
    Assembler class is used as a static class to execute the report of the de novo assembly (quast)
    """
    def __init__(self):
        Main.logger.debug("Assemblers:")
        threading.Thread.__init__(self)
        return


    @staticmethod
    def quast(contigs):
        """
        Quast is used to export the quality of the de novo assembly.
        Before running Quast, a gff file is transformed as a gene indicator for QUAST.
        The location of QUAST is hardcoded. if you want to use this module, update the location to your location

        :param contigs:
        :return:
        """
        Main.logger.debug("Quast: "+contigs)
        geneFile = Assemblers.gffChanger(Main.gffFile)
        workline = "/mnt/apps/quast/quast-2.3/quast.py -o " + Main.resultDir + " -R " + Main.genomeAdd +\
                   Main.refGenomeList[0] + " -G " + geneFile + " " + contigs
        Main.execute(workline, "Running quast, please wait!")
        return

    @staticmethod
    def gffChanger(gffFileName):
        """
        This function changes a gff3 file into a text file readable for QUAST.
        Only the exons are saved. this is hardcoded, but can be changed to any other type
        :param gffFileName:
        :return:
        """
        Main.logger.debug("gffChanger: "+gffFileName)
        i = 0
        gffFile = open(gffFileName, mode="r")
        quastFileName = gffFileName[:-4]+".txt"
        quastFile = open(quastFileName, mode="w")
        for line in gffFile:
            if "exon" in line:
                line = line.split("\t")
                if "-" in line[6]:
                    start = line[4]
                    end = line[3]
                else:
                    start = line[3]
                    end = line[4]
                quastFile.write(line[0] + "\t" + str(i) + "\t" + start + "\t" + end + "\n")
                i += 1
        gffFile.close()
        quastFile.close()
        return quastFileName

class Spades(Assemblers):
    """
    Spades is run in a thread. originally, it was part of the assembler, but this is changed over time.
    Todo:
    - remove from assemblers class and add to threading.Thread class
    """
    def __init__(self, fileForward, fileReversed, workDir):
        """
        Initiation function for spades. the inputs are simple parameters
        :param fileForward: location of the forward read file
        :param fileReversed: location of the reversed read file
        :param workDir: location of the output of spades
        """
        Main.logger.debug("Spades: f." + fileForward + " r." + fileReversed + " w." + workDir)
        threading.Thread.__init__(self)
        self.fileForward = fileForward
        self.fileReversed = fileReversed
        self.workDir = workDir

    def run(self):
        """
        Run is called upon by making an instance of the class spades, and then start(). this is a normal method
        for threaded classes.
        The function calls upon the executer, which runs a terminal command. the command is to make a de novo assembly
        by pairwise reads, and with the output in an unique map for each assembly.
        The location of SPAdes is hardcoded. if you want to use this module, update the location to your location

        :return:
        """
        Main.logger.debug("Spades run:")
        outputName = self.fileForward.split(".")
        outputName = outputName[0]
        self.outputDir = self.workDir+"/"+outputName
        workline = '/mnt/apps/SPAdes-3.1.1-Linux/bin/spades.py -1 ' + Main.fastQAdd + self.fileForward+' -2 ' \
                   + Main.fastQAdd + self.fileReversed+' -o ' + self.outputDir
        Main.execute(workline, "Running spades, please wait")
        os.rename(self.outputDir+"/contigs.fasta", self.outputDir+"/"+outputName+".fa")
        self.outputDir = self.outputDir+"/"+outputName+".fa"
        return
