import subprocess
import threading
from Main import Main


class NUCmer:
    Contigset = set()
    DuplicateSet = set()

    def __init__(self, scaffold):
        self.hitList = list()
        self.scaffold = scaffold
        return

    def hit(self, start, end, ContigID):
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
        returnList = list()
        for hitItem in self.hitList:
            returnList.append(hitItem[startEnd])
        return returnList

    def combineStartEnds(self):
        start = 0
        end = 0
        self.startEnd = list()
        self.starts = list()
        self.ends = list()
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

    @staticmethod
    def execute(cmd, worktext="NUCmer, please wait"):
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        Main.printer(worktext)
        jobNr, stderr = p.communicate()
        if p.returncode == 0:
            print "Done! "
        else:
            print "Error"
            for line in stderr:
                print line
        return

    def __init__(self, contigs):
        threading.Thread.__init__(self)
        self.contigs = contigs
        return

    def run(self):
        workline = "nucmer -maxmatch -p " + self.contigs + " " + \
                   Main.refGenomeList[0] + " " + self.contigs
        NUCmerRun.execute(workline, "Running NUCmer, please wait")
        workline = "show-coords -l " + self.contigs + ".delta > " + self.contigs + ".coords"
        NUCmerRun.execute(workline, "Generating coordinate locations, please wait")
        return