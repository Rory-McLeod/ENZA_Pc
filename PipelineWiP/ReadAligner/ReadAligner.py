#!/usr/bin/env python
import subprocess
import threading
import sys
from Main.Main import Main

class ReadAligner(threading.Thread):

    fastQFile1 = ""
    fastQFile2 = ""
    referenceDB = ""
    samFile = ""

    def __init__(self):
        threading.Thread.__init__(self)
        return

    @staticmethod
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


class Bowtie2(ReadAligner):

    referenceFileName = ""
    outputDir = ""

    def __init__(self, fastQFile1, fastQFile2, referenceFileName, outputDir):
        threading.Thread.__init__(self)
        self.fastQFile1 = fastQFile1
        self.fastQFile2 = fastQFile2
        self.referenceFileName = referenceFileName
        self.outputDir = outputDir
        return

    def run(self):
        Main.printer("Running thread to execute with more speed")
        self.bowTie2Index(self.referenceFileName, self.outputDir)
        self.bowTie2Map(self.outputDir)

    def bowTie2Index(self, referenceFileName, outputdir):
        Main.printer("Working on generating an index database for Bowtie2")
        referenceName = referenceFileName.split(".")
        referenceName = referenceName[0]
        workLine = "bowtie2-build " + referenceFileName + " " + outputdir + "/" + referenceName
        self.execute(workLine, "Making index, please wait")
        self.referenceDB = referenceName
        return

    def bowTie2Map(self, outputDir):
        fastQName = self.fastQFile1.split(".")[0]
        samFile = outputDir + "/Bow" + self.referenceDB + fastQName + ".sam"
        workLine = "bowtie2 -x " + outputDir + "/" + self.referenceDB + \
                   " -1 " + self.fastQFile1 + \
                   " -2 " + self.fastQFile2 + \
                   " -q -S " + samFile
        self.execute(workLine)
        self.samFile = samFile.replace(".sam", "")
        return


class BWA(ReadAligner):

    def BWA(self):
        return


class BamTools(ReadAligner):

    def __init__(self, samFile, referenceDB):
        threading.Thread.__init__(self)
        self.samFile = samFile
        self.referenceDB = referenceDB
        return

    def run(self):
        Main.printer("Running thread to change sam to usable files")
        self.referenceIndex()
        self.samToBam()
        self.bamToBai()
        self.bamToBed()
        self.bamToCoverageRate()
        # self.snpCalling()

    def samToBam(self):
        workLine = "samtools view -bS " + self.samFile + ".sam | samtools sort -o " + self.samFile + ".sorted"
        self.execute(workLine, "Creating sorted BAM file from SAM file")
        return

    def bamToBai(self):
        workLine = "samtools index " + self.samFile + ".sorted " + self.samFile + ".bai"
        self.execute(workLine, "Creating index from BAM file")
        return

    def bamToBed(self):
        workLine = "genomeCoverageBed " \
                   "-ibam " + self.samFile + ".sorted " \
                   "-g workDir/" + self.referenceDB + " -d >" + self.samFile + ".bed"
        self.execute(workLine, "Generating BED file from sorted BAM file")
        return

    def bamToCoverageRate(self):
        workLine = "genomeCoverageBed " \
                   "-ibam " + self.samFile + ".sorted " \
                   "-g workDir/" + self.referenceDB + " > " + self.samFile + ".CovBed"
        self.execute(workLine, "Generating depth info file from BAM")
        return

    def referenceIndex(self):
        workLine = "samtools faidx " + self.referenceDB
        self.execute(workLine, "Generating reference index for samtools")
        return

    def snpCalling(self):
        workLine = "samtools mpileup -ugf " + self.referenceDB + " " + self.samFile + ".sorted |" \
                   " bcftools call -vmO v -o " + self.samFile + ".vcf"
        self.execute(workLine, "Variance calling")
        # workLine = "bcftools view " + self.samFile + ".bcf | vcfutils.pl varFilter -D 8000 > " + self.samFile + ".vcf"
        # self.execute(workLine, "Variance calling, changes BCF to VCF")
        return

