#!/usr/bin/env python

import subprocess
import threading
from Main import Main
import sys
import os


class Assemblers(threading.Thread):
    # TODO: remove execute from Assembler to Main module
    # TODO: change gffChanger for output in Main.workdir instead of program place
    def __init__(self):
        threading.Thread.__init__(self)
        return

    @staticmethod
    def execute(cmd, worktext="Assembling, please wait"):
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
    def quast(contigs):
        geneFile = Assemblers.gffChanger(Main.gffFile)
        workline = "/mnt/apps/quast/quast-2.3/quast.py -o " + Main.resultDir + " -R " + Main.refGenomeList[0] +\
                   " -G " + geneFile + " " + contigs
        print workline
        Assemblers.execute(workline, "Running quast, please wait!")
        return

    @staticmethod
    def gffChanger(gffFileName):
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
    # TODO: change print to Main Printer
    # TODO: run Main execute instead of assembler
    def __init__(self, fileForward, fileReversed, workDir):
        threading.Thread.__init__(self)
        self.fileForward = fileForward
        self.fileReversed = fileReversed
        self.workDir = workDir

    def run(self):
        outputName = self.fileForward.split(".")
        outputName = outputName[0]
        self.outputDir = self.workDir+"/"+outputName
        workline = '/mnt/apps/SPAdes-3.1.1-Linux/bin/spades.py -1 '+self.fileForward+' -2 '+self.fileReversed+' -o ' + \
                   self.outputDir
        Assemblers.execute(workline, "Running spades, please wait")
        print "Changing name"
        os.rename(self.outputDir+"/contigs.fasta", self.outputDir+"/"+outputName+".fa")
        self.outputDir = self.outputDir+"/"+outputName+".fa"
        return
