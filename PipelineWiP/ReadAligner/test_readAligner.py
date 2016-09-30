from unittest import TestCase
import ReadAligner
import os
import os.path
import Main

class TestBowtie2(TestCase):

    def test_bowTie2Index(self):
        DIR = "test_data/reference"
        for name in os.listdir(DIR):
            if len(name) > 15:
                os.remove(os.path.join(DIR, name))
        bowtie2 = ReadAligner.Bowtie2("test_data/reads/reads_1.fq", "test_data/reads/reads_2.fq",
                                      "test_data/reference/lambda_virus.fa", "")
        bowtie2.bowTie2Index(bowtie2.referenceFileName, bowtie2.outputDir)
        self.assertEqual(len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]), 7)
        for name in os.listdir(DIR):
            if len(name) > 15:
                os.remove(os.path.join(DIR, name))

    def test_bowTie2Map(self):
        bowtie2 = ReadAligner.Bowtie2("test_data/reads/reads_1.fq", "test_data/reads/reads_2.fq",
                                      "test_data/reference/lambda_virus.fa", "")
        bowtie2.bowTie2Index(bowtie2.referenceFileName, bowtie2.outputDir)
        bowtie2.bowTie2Map(bowtie2.outputDir)
        self.assertEqual(bowtie2.samFile, "test_data/reads/reads_1")
        os.remove(bowtie2.samFile+".sam")


class TestBamTools(TestCase):

    def test_samToBam(self):
        bowtie2 = ReadAligner.Bowtie2("test_data/reads/reads_1.fq", "test_data/reads/reads_2.fq",
                                      "test_data/reference/lambda_virus.fa", "")
        bowtie2.bowTie2Index(bowtie2.referenceFileName, bowtie2.outputDir)
        bowtie2.bowTie2Map(bowtie2.outputDir)
        bamTools = ReadAligner.BamTools(bowtie2.samFile, bowtie2.referenceDB)
        bamTools.samToBam()
        self.assertTrue(os.path.isfile(bowtie2.samFile+".sorted"))
        os.remove(bowtie2.samFile+".sam")
        os.remove(bowtie2.samFile+".sorted")

    def test_bamToBai(self):
        bowtie2 = ReadAligner.Bowtie2("test_data/reads/reads_1.fq", "test_data/reads/reads_2.fq",
                                      "test_data/reference/lambda_virus.fa", "")
        bowtie2.bowTie2Index(bowtie2.referenceFileName, bowtie2.outputDir)
        bowtie2.bowTie2Map(bowtie2.outputDir)
        bamTools = ReadAligner.BamTools(bowtie2.samFile, bowtie2.referenceDB)
        bamTools.samToBam()
        bamTools.bamToBai()
        self.assertTrue(os.path.isfile(bowtie2.samFile + ".bai"))
        os.remove(bowtie2.samFile + ".sam")
        os.remove(bowtie2.samFile + ".sorted")
        os.remove(bowtie2.samFile + ".bai")

    def test_bamToBed(self):
        bowtie2 = ReadAligner.Bowtie2("test_data/reads/reads_1.fq", "test_data/reads/reads_2.fq",
                                      "test_data/reference/lambda_virus.fa", "")
        bowtie2.bowTie2Index(bowtie2.referenceFileName, bowtie2.outputDir)
        bowtie2.bowTie2Map(bowtie2.outputDir)
        bamTools = ReadAligner.BamTools(bowtie2.samFile, bowtie2.referenceDB)
        bamTools.samToBam()
        bamTools.bamToBed()
        self.assertTrue(os.path.isfile(bowtie2.samFile + ".bed"))
        os.remove(bowtie2.samFile + ".sam")
        os.remove(bowtie2.samFile + ".sorted")
        os.remove(bowtie2.samFile + ".bed")

    def test_bamToCoverageRate(self):
        bowtie2 = ReadAligner.Bowtie2("test_data/reads/reads_1.fq", "test_data/reads/reads_2.fq",
                                      "test_data/reference/lambda_virus.fa", "")
        bowtie2.bowTie2Index(bowtie2.referenceFileName, bowtie2.outputDir)
        bowtie2.bowTie2Map(bowtie2.outputDir)
        bamTools = ReadAligner.BamTools(bowtie2.samFile, bowtie2.referenceDB)
        bamTools.samToBam()
        bamTools.bamToCoverageRate()
        self.assertTrue(os.path.isfile(bowtie2.samFile + ".CovBed"))
        os.remove(bowtie2.samFile + ".sam")
        os.remove(bowtie2.samFile + ".sorted")
        os.remove(bowtie2.samFile + ".CovBed")

    def test_snpCalling(self):
        bowtie2 = ReadAligner.Bowtie2("test_data/reads/reads_1.fq", "test_data/reads/reads_2.fq",
                                      "test_data/reference/lambda_virus.fa", "")
        bowtie2.bowTie2Index(bowtie2.referenceFileName, bowtie2.outputDir)
        bowtie2.bowTie2Map(bowtie2.outputDir)
        bamTools = ReadAligner.BamTools(bowtie2.samFile, bowtie2.referenceDB)
        bamTools.samToBam()
        Main.Main.refGenomeList.append("test_data/reference/lambda_virus.fa")
        bamTools.snpCalling()
        self.assertTrue(os.path.isfile(bowtie2.samFile + ".vcf"))
        os.remove(bowtie2.samFile + ".sam")
        os.remove(bowtie2.samFile + ".sorted")
        os.remove(bowtie2.samFile + ".vcf")

