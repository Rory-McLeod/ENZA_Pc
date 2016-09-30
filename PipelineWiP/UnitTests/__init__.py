import unittest
from Assemblers import test_assemblers
from Blast import test_blast
from VisualisationTools import test_visualisationTools
from ReadAligner import test_readAligner
from Main import test_main
from NUCmer import test_NUCmer
from PrimerDesign import test_primerDesign

class UnitTests:

    def __init__(self):
        return

    def runAssemblerTest(self):
        print "*" * 70
        print "Running tests on Assembler module"
        print "=" * 70
        unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().loadTestsFromModule(test_assemblers))
        return

    def runBlastTest(self):
        print "*" * 70
        print "Running tests on Blast module"
        print "=" * 70
        unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().loadTestsFromModule(test_blast))

    def runVisualisationTools(self):
        print "*" * 70
        print "Running tests on VisualisationTools module"
        print "=" * 70
        unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().loadTestsFromModule(test_visualisationTools))

    def runReadAlignerTest(self):
        print "*" * 70
        print "Running tests on ReadAligner module"
        print "=" * 70
        unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().loadTestsFromModule(test_readAligner))

    def runMainTest(self):
        print "*" * 70
        print "Running tests on Main module"
        print "=" * 70
        unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().loadTestsFromModule(test_main))

    def runNUCmerTest(self):
        print "*" * 70
        print "Running tests on NUCmer module"
        print "=" * 70
        unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().loadTestsFromModule(test_NUCmer))

    def runPrimerDesignTest(self):
        print "*" * 70
        print "Running tests on PrimerDesign module"
        print "=" * 70
        unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().loadTestsFromModule(test_primerDesign))
