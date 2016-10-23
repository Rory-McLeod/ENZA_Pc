#!/usr/bin/env python

import threading
from Main import Main
import os


class Assemblers(threading.Thread):
    # TODO: change gffChanger for output in Main.workdir instead of program place
    def __init__(self):
        Main.logger.debug("Assemblers:")
        threading.Thread.__init__(self)
        return

    @staticmethod
    def quast(contigs):
        Main.logger.debug("Quast: "+contigs)
        geneFile = Assemblers.gffChanger(Main.gffFile)
        workline = "/mnt/apps/quast/quast-2.3/quast.py -o " + Main.resultDir + " -R " + Main.genomeAdd +\
                   Main.refGenomeList[0] + " -G " + geneFile + " " + contigs
        Main.execute(workline, "Running quast, please wait!")
        return

    @staticmethod
    def gffChanger(gffFileName):
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
    def __init__(self, fileForward, fileReversed, workDir):
        Main.logger.debug("Spades: f." + fileForward + " r." + fileReversed + " w." + workDir)
        threading.Thread.__init__(self)
        self.fileForward = fileForward
        self.fileReversed = fileReversed
        self.workDir = workDir

    def run(self):
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
