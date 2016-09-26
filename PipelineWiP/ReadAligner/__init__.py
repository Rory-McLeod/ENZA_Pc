#!/usr/bin/env python
"""
Created on: 29-08-2016
@author: H.J.C. Cornelisse
Class is used to map the reads against the reference genome.
Todo:
- Write unit testing
- Write BWA aligner
"""

import subprocess
import threading
import sys
from Main import Main

class ReadAligner(threading.Thread):

    fastQFile1 = ""
    fastQFile2 = ""
    referenceDB = ""
    samFile = ""

    def __init__(self):
        """
        Method for initiating the read aligner.
        threading.Thread is called for threaded use of this class
        """
        threading.Thread.__init__(self)
        return

    @staticmethod
    def execute(cmd, worktext="Mapping, please wait"):
        """
        Method to execute functions on command line.
        Args:
            cmd: the command line command
            worktext: the printed text to the user

        Returns:
            None: returns to the place of calling
        """
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        Main.printer(worktext)
        stderr = p.communicate()[1]
        if p.returncode == 0:
            print "Done! "
        else:
            print "Error"
            for line in stderr:
                print line
            sys.exit(1)
        Main.printer(stderr)
        return


class Bowtie2(ReadAligner):

    referenceFileName = ""
    outputDir = ""

    def __init__(self, fastQFile1, fastQFile2, referenceFileName, outputDir):
        """
        Method for initiating the read aligner, specific Bowtie2.
        threading.Thread is called for threaded use of this class
        Args:
            fastQFile1: location of the fastq file (forward)
            fastQFile2: location of the fastq file (reversed)
            referenceFileName: location of the reference genome
            outputDir: location of the output of these systems (default = workDir)
        """
        threading.Thread.__init__(self)
        self.fastQFile1 = fastQFile1
        self.fastQFile2 = fastQFile2
        self.referenceFileName = referenceFileName
        self.outputDir = outputDir
        return

    def run(self):
        """
        Method called upon by start(), starts the threaded use.
        The method runs two methods in a single thread
        Returns:
            None: returns to the place of calling
        """
        Main.printer("Running thread to execute with more speed")
        self.bowTie2Index(self.referenceFileName, self.outputDir)
        self.bowTie2Map(self.outputDir)

    def bowTie2Index(self, referenceFileName, outputdir):
        """
        Generates an index database for mapping the reads against (for use in Bowtie2
        Args:
            referenceFileName: the reference genome full file name, input by user
            outputdir: the location to save the index database

        Returns:
            None: returns to the place of calling
        """
        Main.printer("Working on generating an index database for Bowtie2")
        referenceName = referenceFileName.split(".")
        referenceName = referenceName[0]
        workLine = "bowtie2-build " + referenceFileName + " " + outputdir + "/" + referenceName
        self.execute(workLine, "Making index, please wait")
        self.referenceDB = referenceName
        return

    def bowTie2Map(self, outputDir):
        """
        Method to map the reads against the generated index database by Bowtie2
        Args:
            outputDir: the location to save the SAM file

        Returns:
            None: returns to the place of calling
        """
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
    """
    Todo: work in progress
    """

    def BWA(self):
        return


class BamTools(ReadAligner):

    def __init__(self, samFile, referenceDB):
        """
        Method for initiating the read aligner, specific bam/sam/bed tools.
        threading.Thread is called for threaded use of this class
        Args:
            samFile: the workname of the output of bowtie2/bwa for further use
            referenceDB: the reference database used to map the reads against

        Returns:
            None: returns to the place of calling
        """
        threading.Thread.__init__(self)
        self.samFile = samFile
        self.referenceDB = referenceDB
        return

    def run(self):
        """
        Method called upon by start(), starts the threaded use.
        The method runs five methods in a single thread
        Returns:
            None: returns to the place of calling
        """
        Main.printer("Running thread to change sam to usable files")
        self.samToBam()
        self.bamToBai()
        self.bamToBed()
        self.bamToCoverageRate()
        self.snpCalling()

    def samToBam(self):
        """
        Transfers the SAM file to a sorted BAM file. The sorted BAM file is important to speed up other processes
        Returns:
            None: returns to the place of calling
        """
        workLine = "samtools view -bS " + self.samFile + ".sam | samtools sort -o " + self.samFile + ".sorted"
        self.execute(workLine, "Creating sorted BAM file from SAM file")
        return

    def bamToBai(self):
        """
        Generates a index file from a sorted BAM file for speed access to certain position in the genome
        Returns:
            None: returns to the place of calling
        """
        workLine = "samtools index " + self.samFile + ".sorted " + self.samFile + ".bai"
        self.execute(workLine, "Creating index from BAM file")
        return

    def bamToBed(self):
        """
        Generates a depth coverage per position (BED format) file from a sorted BAM file
        Returns:
            None: returns to the place of calling
        """
        workLine = "genomeCoverageBed " \
                   "-ibam " + self.samFile + ".sorted " \
                   "-g workDir/" + self.referenceDB + " -d >" + self.samFile + ".bed"
        self.execute(workLine, "Generating BED file from sorted BAM file")
        return

    def bamToCoverageRate(self):
        """
        Generates a file (CovBed) that shows the occurrence of a read depth per scaffold/chromosom and the whole genome
        Returns:
            None: returns to the place of calling
        """
        workLine = "genomeCoverageBed " \
                   "-ibam " + self.samFile + ".sorted " \
                   "-g workDir/" + self.referenceDB + " > " + self.samFile + ".CovBed"
        self.execute(workLine, "Generating depth info file from BAM")
        return

    def snpCalling(self):
        """
        Finds the SNP between the reference and the mapped reads.
        Returns:
            None: returns to the place of calling
        """
        workLine = "samtools mpileup -ugf " + Main.refGenomeList[0] + " " + self.samFile + ".sorted |" \
                   " bcftools call -vmO v -o " + self.samFile + ".vcf"
        self.execute(workLine, "Variance calling")
        return

