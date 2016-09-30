from unittest import TestCase
import PrimerDesign
import Main
import os
import VisualisationTools

# class TestPrimerDesign(TestCase):

    # def test_readRefGenome(self):
    #     genome = PrimerDesign.PrimerDesign.readRefGenome("test_data/reference/lambda_virus.fa")
    #     self.assertEqual(len(genome["gi|9626243|ref|NC_001416.1| Enterobacteria phage lambda, complete genome"]), 48502)

    # def test_readGFF(self):
    #     genome = PrimerDesign.PrimerDesign.readRefGenome("test_data/reference.fasta")
    #     POI = PrimerDesign.PrimerDesign.readGFF("test_data/genes.gff", genome)
    #     self.assertEqual(len(POI), 9)

    # def test_saveFasta(self):
    #     Main.Main.workDir = "test_data"
    #     outputProduct = "test.fasta"
    #     genome = PrimerDesign.PrimerDesign.readRefGenome("test_data/reference.fasta")
    #     POI = PrimerDesign.PrimerDesign.readGFF("test_data/genes.gff", genome)
    #     PrimerDesign.PrimerDesign.saveFasta(outputProduct, POI)
    #     self.assertTrue(os.path.isfile("test_data/test.fasta"))
    #     os.remove("test_data/test.fasta")

    # def test_generatePrimer3Input(self):
    #     self.fail()


class TestPrimerDesignByMapping(TestCase):
    # def test_generateCoords(self):
    #     testResult = True
    #     testAgainst = {"PHYCAscaffold_5": 3779, "PHYCAscaffold_4": 4699,
    #                    "PHYCAscaffold_6": 2697, "PHYCAscaffold_1": 5777,
    #                    "PHYCAscaffold_3": 4252, "PHYCAscaffold_2": 4932}
    #     visualisationTools = VisualisationTools.VisualisationTools("test_data/test", "test_data/")
    #     visualisationTools.readBedToLocal()
    #     primerMapping = PrimerDesign.PrimerDesignByMapping()
    #     primerMapping.generateCoords(visualisationTools.depthPerPos)
    #     for scaffoldDict in primerMapping.coordsInfo:
    #         for key, value in scaffoldDict.iteritems():
    #             if len(value.starts) != testAgainst[key]:
    #                 print key
    #                 testResult = False
    #     self.assertTrue(testResult)

    # def test_readGenes(self):
    #     testResult = True
    #     testAgainst = {"PHYCAscaffold_1": 2963, "PHYCAscaffold_2": 2673}
    #     primerMapping = PrimerDesign.PrimerDesignByMapping()
    #     primerMapping.readGenes("test_data/test.gff")
    #     for key, value in primerMapping.geneInfo.iteritems():
    #         if len(value.starts) != testAgainst[key]:
    #             print key
    #             testResult = False
    #     self.assertTrue(testResult)

    def test_generateHitList(self):
        testResult = True
        testAgainst = {"PHYCAscaffold_5": 3779, "PHYCAscaffold_4": 4699,
                       "PHYCAscaffold_6": 2697, "PHYCAscaffold_1": 5777,
                       "PHYCAscaffold_3": 4252, "PHYCAscaffold_2": 4932}
        visualisationTools = VisualisationTools.VisualisationTools("test_data/test", "test_data/")
        visualisationTools.readBedToLocal()
        primerMapping = PrimerDesign.PrimerDesignByMapping()
        primerMapping.generateCoords(visualisationTools.depthPerPos)
        primerMapping.generateHitList()
        for scaffold, value in primerMapping.HitList.iteritems():
            if len(value) != testAgainst[scaffold]:
                print scaffold
                testResult = False
        self.assertTrue(testResult)

#
#     def test_generateGeneSpecificHitList(self):
#         self.fail()
#
# class TestPrimerDesignByDenovo(TestCase):
#     def test_readCoords(self):
#         self.fail()
#
#     def test_readGenes(self):
#         self.fail()
#
#     def test_generateHitList(self):
#         self.fail()
#
#     def test_generateGeneSpecificHitList(self):
#         self.fail()
