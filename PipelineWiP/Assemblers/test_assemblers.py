"""
Created on: 29-08-2016
Latest update: 02-11-2016
@author: H.J.C. Cornelisse
Module is used to test the assembler method
Dependencies:
- QUAST 2.3
- SPAdes 3.1.1.
Todo:
- add more tests
- add comments
"""
import unittest
import subprocess
from Assemblers import Assemblers


class TestAssemblers(unittest.TestCase):

    def test_quast(self):
        """
        Calls upon the QUAST internal test module. Requires the test files of QUAST.
        Before testing, make sure that quast_test_output is empty.

        Todo:
        - empty quast_test_output before running test
        :return:
        """
        workline = "/mnt/apps/quast/quast-2.3/quast.py --test"
        p = subprocess.Popen(workline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout, stderr = p.communicate()
            if p.returncode == 0:
                logFile = open("quast_test_output/quast.log", mode='r')
                testResult = False
                for line in logFile:
                    if "TEST PASSED" in line:
                        testResult = True
                self.assertTrue(testResult)
            else:
                raise EnvironmentError(stderr)
        except IOError as error:
            self.fail(error)
        except EnvironmentError as error:
            self.fail(error)

    def test_gffChanger(self):
        """
        Test if gffChanger can open a file, read the content, and saves the file
        Before calling text, make sure that the gffchanger.gff file is in test_data and that the output does not exist.

        Todo:
        - remove gffchanger.txt before running test
        :return:
        """
        gffFile = "test_data/gffchanger.gff"
        Assemblers.gffChanger(gffFile)
        resultListStarts = [41144, 41689, 41195, 41687]
        try:
            outputFile = open("test_data/gffchanger.txt")
            testResult = True
            for line in outputFile:
                line = line.split("\t")
                if "PHYCAscaffold_1" not in line[0]:
                    testResult = False
                IDLine = line[1].split()
                idTag = int(IDLine[0].replace("ID=", ""))
                if int(line[2]) != resultListStarts[idTag]:
                    testResult = False
            self.assertTrue(testResult)
        except IOError as error:
            self.fail(error)


class TestSpades(unittest.TestCase):

    def test_run(self):
        """
        Calls upon SPAdes
        :return:
        """
        workline = "/mnt/apps/SPAdes-3.1.1-Linux/bin/spades.py --test"
        p = subprocess.Popen(workline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout, stderr = p.communicate()
            if p.returncode == 0:
                logFile = open("spades_test/spades.log", mode='r')
                testResult = False
                for line in logFile:
                    if "========= TEST PASSED CORRECTLY." in line:
                        testResult = True
                self.assertTrue(testResult)
            else:
                raise EnvironmentError(stderr)
        except IOError as error:
            self.fail(error)
        except EnvironmentError as error:
            self.fail(error)
