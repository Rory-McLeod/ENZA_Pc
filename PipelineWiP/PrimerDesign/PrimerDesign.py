#!/usr/bin/env python

import subprocess
import threading
import sys
from Main.Main import Main

class PrimerDesign(threading.Thread):

    depthPerPos = []
    POI = []
    Genome = []
    referenceDB = ""

    def __init__(self, depthPerPos, referenceDB):
        threading.Thread.__init__(self)
        self.depthPerPos = depthPerPos
        self.referenceDB = referenceDB
        return

    def run(self):
        self.generatePOI()
        # self.readRefGenome()
        # self.savePOI()
        return

    def execute(cmd, worktext="Mapping, please wait"):
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

    def generatePOI(self, depthLimit=10, nucLimit=28):
        Main.printer("generate POI's")
        counter = 0
        beginPos = 0
        i = 0
        dPP = self.depthPerPos
        length = len(dPP)
        while i < length:
            if int(dPP[i]) >= depthLimit:
                if counter == 0:
                    beginPos = i
                counter += 1
            else:
                if counter > 0:
                    if counter >= nucLimit:
                        endPos = i
                        self.POI.append([beginPos, endPos])
                    counter = 0
            i += 1
        for listing in self.POI:
            Main.printer(str(listing[0]+1) + " - " + str(listing[1]+1) + ": " + dPP[listing[0]].rstrip() + " - " + dPP[listing[1]-1].rstrip())
        return

    def readRefGenome(self):
        Main.printer("read reference Genome")
        self.Genome = []
        referenceFile = file(self.referenceDB)
        for line in referenceFile:
            line = line.rstrip()
            if ">" not in line:
                for seqChar in line:
                    self.Genome.append(seqChar)
        return

    def savePOI(self):
        Main.printer("print POI sequences")
        counter = 0
        for positions in self.POI:
            sequence = ""
            beginPos = positions[0]
            endPos = positions[1]
            while beginPos < endPos:
                sequence += str(self.Genome[beginPos])
                beginPos += 1
            print str(counter) + ": sequence: "+sequence
            counter += 1
        return
