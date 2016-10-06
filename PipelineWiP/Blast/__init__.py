import threading
import subprocess
from Main import Main
import sys

class Blast(threading.Thread):
    y = 0

    def __init__(self, query, genomeList):
        Blast.y += 1
        threading.Thread.__init__(self)
        self.genomeList = genomeList
        for genome in self.genomeList:
            print genome
        self.query = query
        return

    @staticmethod
    def execute(cmd, worktext="Primerdesign, please wait"):
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        Main.printer(worktext)
        jobNr, stderr = p.communicate()
        if p.returncode == 0:
            print "Done! "
        else:
            print "Error"
            print stderr
            sys.exit(1)
        return

    def run(self):
        allAlias = ""
        for genome in self.genomeList:
            self.makeDatabase(str(genome))
            allAlias += genome + " "
        allAlias = allAlias.rstrip()
        self.aliasTool(allAlias)
        self.doBlast()
        self.interpertBlast()
        return

    def makeDatabase(self, genomeFile):
        workline = "makeblastdb -in " + genomeFile + " -parse_seqids -dbtype nucl -out " + genomeFile
        Blast.execute(workline)
        return

    def aliasTool(self, allAlias, outputDir=Main.workDir):
        output = outputDir+"/others"
        workline = "blastdb_aliastool -dblist \"" + allAlias + "\" -dbtype nucl -out "+output + " -title others"
        print workline
        Blast.execute(workline, "Making one database")
        return

    def doBlast(self, outputDir=Main.workDir):
        workline = "blastn -db " + outputDir + "/others -query " + self.query + " -out " + outputDir +\
                   "/blastResult.csv -outfmt \"10 std\""
        Blast.execute(workline, "Running blast")
        return

    def interpertBlast(self, workDir=Main.workDir):
        blastFile = open(workDir+"/blastResult"+ Blast.y +".csv")
        self.hitSet = set()
        for line in blastFile:
            line = line.split(",")
            self.hitSet.add(line[0])
        return
