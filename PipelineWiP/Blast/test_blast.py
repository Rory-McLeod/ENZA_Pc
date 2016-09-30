from unittest import TestCase
import subprocess
from Blast import Blast

class TestBlast(TestCase):

    def test_blastVersion(self):
        workline = "blastn -version"
        p = subprocess.Popen(workline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout, stderr = p.communicate()
            testResult = True
            if p.returncode == 0:
                if "2.2.28" not in stdout:
                    testResult = False
                self.assertTrue(testResult)
            else:
                raise EnvironmentError(stderr)
        except EnvironmentError as error:
            self.fail(error)


    def test_makeDatabase(self):
        tempBlast = Blast("", "")
        tempBlast.makeDatabase("test_data/reference.fasta")
        workline = "blastdbcmd -db test_data/reference.fasta -info"
        p = subprocess.Popen(workline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout, stderr = p.communicate()
            testResult = True
            if p.returncode == 0:
                if "1 sequences; 10,000 total bases" not in stdout:
                    testResult = False
                self.assertTrue(testResult)
            else:
                raise EnvironmentError(stderr)
        except EnvironmentError as error:
            self.fail(error)

    def test_aliasTool(self):
        tempBlast = Blast("", "")
        tempBlast.makeDatabase("test_data/reference.fasta")
        tempBlast.aliasTool("test_data/reference.fasta", "test_data")
        workline = "blastdbcmd -db test_data/others -info"
        p = subprocess.Popen(workline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout, stderr = p.communicate()
            testResult = True
            if p.returncode == 0:
                if "1 sequences; 10,000 total bases" not in stdout:
                    testResult = False
                self.assertTrue(testResult)
            else:
                raise EnvironmentError(stderr)
        except EnvironmentError as error:
            self.fail(error)

    def test_doBlast(self):
        tempBlast = Blast("", "")
        tempBlast.makeDatabase("test_data/reference.fasta")
        tempBlast.aliasTool("test_data/reference.fasta", "test_data")
        tempBlast.query = "test_data/contigs_1.fasta"
        tempBlast.doBlast("test_data")
        try:
            outputFile = open("test_data/blastResult.csv")
            testResult = True
            for line in outputFile:
                if "gi|49175990|ref|NC_000913.2|" not in line:
                    testResult = False
            self.assertTrue(testResult)
        except IOError as error:
            self.fail(error)

    def test_interpertBlast(self):
        tempBlast = Blast("", "")
        tempBlast.makeDatabase("test_data/reference.fasta")
        tempBlast.aliasTool("test_data/reference.fasta", "test_data")
        tempBlast.query = "test_data/contigs_1.fasta"
        tempBlast.doBlast("test_data")
        tempBlast.interpertBlast("test_data")
        self.assertEqual(len(tempBlast.hitSet), 3)
