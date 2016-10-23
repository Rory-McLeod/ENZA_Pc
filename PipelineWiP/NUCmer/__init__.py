#!/usr/bin/env python

import threading
from Main import Main


class NUCmer:
    Contigset = set()
    DuplicateSet = set()

    def __init__(self, scaffold):
        Main.logger.debug("Nucmer:"+scaffold)
        self.hitList = list()
        self.scaffold = scaffold
        return

    def hit(self, start, end, ContigID):
        Main.logger.debug("Nucmer hit: S"+str(start) + " E" + str(end) + " " + ContigID)
        duplicate = False
        if ContigID in NUCmer.Contigset:
            duplicate = True
        hitItem = [start, end, ContigID, duplicate]
        if duplicate:
            NUCmer.DuplicateSet.add(ContigID)
        NUCmer.Contigset.add(ContigID)
        self.hitList.append(hitItem)
        return

    def getStartEnds(self, startEnd=0):
        Main.logger.debug("Nucmer getStartEnds: "+startEnd)
        returnList = list()
        for hitItem in self.hitList:
            returnList.append(hitItem[startEnd])
        return returnList

    def combineStartEnds(self):
        Main.logger.debug("Nucmer combineStartEnds:")
        start = 0
        end = 0
        self.startEnd = list()
        self.starts = list()
        self.ends = list()
        self.hitList.sort(key=lambda x: int(x[0]))
        for hitItem in self.hitList:
            if end > hitItem[0] >= start:
                if hitItem[1] > end:
                    end = hitItem[1]
            else:
                if start != 0:
                    self.starts.append(start)
                    self.ends.append(end)
                    self.startEnd.append([start, end])
                start = hitItem[0]
                end = hitItem[1]
        self.startEnd.append([start, end])
        self.starts.append(start)
        self.ends.append(end)

class NUCmerRun(threading.Thread):

    def __init__(self, contigs):
        Main.logger.debug("NucmerRun: "+contigs)
        threading.Thread.__init__(self)
        self.contigs = contigs
        return

    def run(self):
        Main.logger.debug("NucmerRun run:")
        workline = "nucmer -maxmatch -p " + self.contigs + " " + \
                   Main.genomeAdd + Main.refGenomeList[0] + " " + self.contigs
        Main.execute(workline, "Running NUCmer, please wait")
        workline = "show-coords -l " + self.contigs + ".delta > " + self.contigs + ".coords"
        Main.execute(workline, "Generating coordinate locations, please wait")
        workline = "show-snps -r " + self.contigs + ".delta > " + self.contigs + ".snps"
        Main.execute(workline, "Generation snp information, please wait")
        return