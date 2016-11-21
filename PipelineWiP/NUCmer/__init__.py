#!/usr/bin/env python
"""
Created on: 29-08-2016
@author: H.J.C. Cornelisse
Class is used to do functions regarding mapping of contigs and saving the locations
Dependencies:
- mummer 3.23 (tested)
"""

import threading
from Main import Main


class NUCmer:
    """
    class used to save coverage locations from primerdesign
    """
    Contigset = set()
    DuplicateSet = set()

    def __init__(self, scaffold):
        """
        coverage is saved by scaffold for easy access
        :param scaffold: scaffold name of the contig
        """
        Main.logger.debug("Nucmer:"+scaffold)
        self.hitList = list()
        self.scaffold = scaffold
        return

    def hit(self, start, end, ContigID):
        """
        information about the coverage, like the start and end position, and which contig the coverage originates from
        The hits are saved in sets
        :param start: start position of coverage
        :param end: end position of coverage
        :param ContigID: name of contig
        :return:
        """
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
        """
        Returns a list of starts or ends, depending on the entry
        :param startEnd: 0 for starts, 1 for ends
        :return: list of all starts or ends
        """
        Main.logger.debug("Nucmer getStartEnds: "+startEnd)
        returnList = list()
        for hitItem in self.hitList:
            returnList.append(hitItem[startEnd])
        return returnList

    def combineStartEnds(self):
        """
        combines overlapping regions that have the same overlap, but different starts or ends (or both)
        :return:
        """
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
    """
    Thread that runs the nucmer tool to align the contigs against the reference genome
    Todo:
    - mapview is not working at this moment due to the output problem
    """

    def __init__(self, contigs):
        """
        initiation step for running nucmer. the system saves the location of the contig that is used in this nucmer run
        :param contigs: location to contig
        """
        Main.logger.debug("NucmerRun: "+contigs)
        threading.Thread.__init__(self)
        self.contigs = contigs
        return

    def run(self):
        """
        running of the nucmer alignment.
        first, the contigs are aligned, then the coordinations are saved for snp calling.
        when the snps are called, a new coords file is generated for coverage and mapview
        :return:
        """
        Main.logger.debug("NucmerRun run:")
        workline = "nucmer -maxmatch -p " + self.contigs + " " + \
                   Main.genomeAdd + Main.refGenomeList[0] + " " + self.contigs
        Main.execute(workline, "Running NUCmer, please wait")
        workline = "show-coords -l " + self.contigs + ".delta > " + self.contigs + ".coords"
        Main.execute(workline, "Generating coordinate locations, please wait")
        workline = "delta-filter -r -q " + self.contigs + ".delta > " + self.contigs + ".filter"
        Main.execute(workline, "Generating coordinate locations, please wait")
        workline = "show-snps -r -I -T -H " + self.contigs + ".filter > " + self.contigs + ".snps"
        Main.execute(workline, "Generation snp information, please wait")
        workline = "show-coords -r -c -l " + self.contigs + ".delta > " + self.contigs + ".Mapcoords"
        Main.execute(workline, "Generating coordinate locations, please wait")
        workline = "mapview -f pdf -p " + self.contigs + " " + self.contigs + ".Mapcoords"
        Main.execute(workline, "Generating mapview in pdf, please wait")
        return