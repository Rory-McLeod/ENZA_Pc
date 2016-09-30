from unittest import TestCase
import Main
import os

class TestMain(TestCase):
    def test_openFastQFiles(self):
        fastQList = ["test_data/reads/longreads.fq", "test_data/reads/longreads.fq.gz"]
        worker = Main.Main()
        worker.openFastQFiles(fastQList)
        self.assertEqual(len(worker.fastQFileList), 2)

    def test_openRefGenomes(self):
        refGenomeList = ["test_data/reference/lambda_virus.fa"]
        worker = Main.Main()
        worker.openRefGenomes(refGenomeList)
        self.assertEqual(len(worker.refGenomeList), 1)

    def test_makeDirectory(self):
        DIR = "Test_Dir_Make"
        worker = Main.Main()
        worker.makeDirectory(DIR)
        self.assertTrue(os.path.isdir(DIR))
        os.rmdir(DIR)
