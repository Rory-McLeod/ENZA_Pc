from unittest import TestCase
import NUCmer
import Main
import os

class TestNUCmer(TestCase):
    def test_hit(self):
        nucmer = NUCmer.NUCmer("testContig")
        nucmer.hit(10, 30, 1)
        nucmer.hit(15, 60, 2)
        nucmer.hit(80, 100, 3)
        self.assertEqual(len(nucmer.hitList), 3)

    def test_getStartEnds(self):
        nucmer = NUCmer.NUCmer("testContig")
        testResult = True
        TestList = [10, 15, 80]
        nucmer.hit(10, 30, 1)
        nucmer.hit(15, 60, 2)
        nucmer.hit(80, 100, 3)
        startList = nucmer.getStartEnds()
        for start in startList:
            if start not in TestList:
                testResult = False
        self.assertTrue(testResult)

    def test_combineStartEnds(self):
        nucmer = NUCmer.NUCmer("testContig")
        testResult = True
        TestList = [10, 80]
        nucmer.hit(10, 30, 1)
        nucmer.hit(15, 60, 2)
        nucmer.hit(80, 100, 3)
        nucmer.combineStartEnds()
        startList = nucmer.starts
        for start in startList:
            if start not in TestList:
                print start
                testResult = False
        self.assertTrue(testResult)

class TestNUCmerRun(TestCase):
    def test_run(self):
        nucmer = NUCmer.NUCmerRun("test_data/contigs_1.fasta")
        Main.Main.refGenomeList.append("test_data/reference.fasta")
        nucmer.run()
        self.assertTrue(os.path.isfile("test_data/contigs_1.fasta.coords"))
