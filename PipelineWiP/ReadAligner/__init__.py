#!/usr/bin/env python
"""
Created on: 29-08-2016
@author: H.J.C. Cornelisse
Class is used to map the reads against the reference genome.
Dependencies:
- Bowtie2 2.2.1 (tested)
"""

import threading
from Main import Main

class ReadAligner(threading.Thread):
    """
    class to work with the reads in the pipeline for direct mapping method
    """

    fastQFile1 = ""
    fastQFile2 = ""
    referenceDB = ""
    samFile = ""

    def __init__(self):
        """
        Just the start of another day
        """
        Main.logger.debug("ReadAligner:")
        threading.Thread.__init__(self)
        return

class Bowtie2(ReadAligner):
    """
    class containing all steps needed for bowtie2
    """
    referenceFileName = ""
    outputDir = ""

    def __init__(self, fastQFile1, fastQFile2, referenceFileName, outputDir):
        """
        Method for initiating the read aligner, specific Bowtie2.
        threading.Thread is called for threaded use of this class
        :param fastQFile1: location of the fastq file (forward)
        :param fastQFile2: location of the fastq file (reversed)
        :param referenceFileName: location of the reference genome
        :param outputDir: location of the output of these systems (default = workDir)
        """
        Main.logger.debug("Bowtie2: f1." + fastQFile1 + " f2." + fastQFile2 +
                          " r." + referenceFileName + " o." + outputDir)
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
        :return:
        """
        Main.logger.debug("Bowtie2 run:")
        Main.logger.info("Running thread to execute with more speed")
        self.bowTie2Index(self.referenceFileName, self.outputDir+"/")
        self.bowTie2Map(self.outputDir+"/")

    def bowTie2Index(self, referenceFileName, outputdir):
        """
        Generates an index database for mapping the reads against (for use in Bowtie2)
        :param referenceFileName: the reference genome full file name, input by user
        :param outputdir: the location to save the index database
        :return:
        """
        Main.logger.debug("Bowtie2 bowTie2Index: r." + referenceFileName + " o." + outputdir)
        Main.logger.info("Working on generating an index database for Bowtie2")
        referenceName = referenceFileName.split(".")
        referenceName = referenceName[0]
        workLine = "bowtie2-build " + Main.genomeAdd + referenceFileName + " " + outputdir + referenceName
        Main.execute(workLine, "Making index, please wait")
        self.referenceDB = referenceName
        return

    def bowTie2Map(self, outputDir):
        """
        Method to map the reads against the generated index database by Bowtie2
        :param outputDir: the location to save the SAM file
        :return:
        """
        Main.logger.debug("Bowtie2 bowtie2map: o." + outputDir)
        fastQName = self.fastQFile1.split(".")[0]
        if outputDir == "":
            samFile = fastQName + ".sam"
        else:
            samFile = outputDir + "/Bow" + self.referenceDB + fastQName + ".sam"
        workLine = "bowtie2 -x " + outputDir + self.referenceDB + \
                   " -1 " + Main.fastQAdd + self.fastQFile1 + \
                   " -2 " + Main.fastQAdd + self.fastQFile2 + \
                   " -a --score-min \"L,-0.06,-0.06\" --phred33 " \
                   "--fr --very-sensitive --un-conc-gz " + fastQName + ".unalign.gz -S " + samFile
        Main.execute(workLine, "Running Bowtie2, please wait")
        self.samFile = samFile.replace(".sam", "")
        return


class BamTools(ReadAligner):
    """
    class used for analysing the read alignment
    """

    def __init__(self, samFile, referenceDB):
        """
        Method for initiating the read aligner, specific bam/sam/bed tools.
        threading.Thread is called for threaded use of this class
        :param samFile: the workname of the output of bowtie2/bwa for further use
        :param referenceDB: the reference database used to map the reads against
        """
        Main.logger.debug("BamTools: s." + samFile + " r." + referenceDB)
        threading.Thread.__init__(self)
        self.samFile = samFile
        self.referenceDB = referenceDB
        return

    def run(self):
        """
        Method called upon by start(), starts the threaded use.
        The method runs five methods in a single thread
        :return:
        """
        Main.logger.debug("BamTools run:")
        Main.logger.info("Running thread to change sam to usable files")
        self.samToBam()
        self.bamToBai()
        self.bamToBed()
        self.bamToCoverageRate()
        self.snpCalling()

    def samToBam(self):
        """
        Transfers the SAM file to a sorted BAM file. The sorted BAM file is important to speed up other processes
        :return:
        """
        Main.logger.debug("BamTools samToBam:")
        workLine = "samtools view -bS " + self.samFile + ".sam | samtools sort -o " + self.samFile + ".sorted"
        Main.execute(workLine, "Creating sorted BAM file from SAM file")
        return

    def bamToBai(self):
        """
        Generates a index file from a sorted BAM file for speed access to certain position in the genome
        :return:
        """
        Main.logger.debug("Bamtools bamToBai:")
        workLine = "samtools index " + self.samFile + ".sorted " + self.samFile + ".bai"
        Main.execute(workLine, "Creating index from BAM file")
        return

    def bamToBed(self):
        """
        Generates a depth coverage per position (BED format) file from a sorted BAM file
        :return:
        """
        Main.logger.debug("Bamtools bamToBed:")
        workLine = "genomeCoverageBed " \
                   "-ibam " + self.samFile + ".sorted " \
                   "-g workDir/" + self.referenceDB + " -d >" + self.samFile + ".bed"
        Main.execute(workLine, "Generating BED file from sorted BAM file")
        return

    def bamToCoverageRate(self):
        """
        Generates a file (CovBed) that shows the occurrence of a read depth per scaffold/chromosom and the whole genome
        :return:
        """
        Main.logger.debug("Bamtools bamToCoverageRate:")
        workLine = "genomeCoverageBed " \
                   "-ibam " + self.samFile + ".sorted " \
                   "-g workDir/" + self.referenceDB + " > " + self.samFile + ".CovBed"
        Main.execute(workLine, "Generating depth info file from BAM")
        return

    def snpCalling(self):
        """
        Finds the SNP between the reference and the mapped reads.
        :return:
        """
        Main.logger.debug("Bamtools snpCalling:")
        workLine = "samtools mpileup -ugf " + Main.genomeAdd +  Main.refGenomeList[0] + " " + self.samFile +\
                   ".sorted | bcftools call -vmO v -o " + self.samFile + ".vcf"
        Main.execute(workLine, "Variance calling")
        return

