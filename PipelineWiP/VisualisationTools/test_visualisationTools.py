from unittest import TestCase
from VisualisationTools import VisualisationTools
from VisualisationTools import Mapping


class TestVisualisationTools(TestCase):
    def test_readBedToLocal(self):
        visualisationTools = VisualisationTools("test_data/test", "test_data/")
        visualisationTools.readBedToLocal()
        self.assertEqual(len(visualisationTools.depthPerPos), 6)

    def test_readCovBedtoLocal(self):
        visualisationTools = VisualisationTools("test_data/test", "test_data/")
        visualisationTools.readCovBedtoLocal()
        self.assertEqual(len(visualisationTools.depthOccurance), 6651)

    def test_getCoveragePercentage(self):
        visualisationTools = VisualisationTools("test_data/test", "test_data/")
        visualisationTools.readCovBedtoLocal()
        self.assertAlmostEqual(visualisationTools.getCoveragePercentage(5), 0.65295766, 4)

class TestMapping(TestCase):
    def test_hit(self):
        mapper = Mapping("test")
        mapper.hit(1, 10, 1)
        mapper.hit(2, 19, 2)
        mapper.hit(25, 30, 2)
        self.assertEqual(len(mapper.hitList), 3)

    def test_getStartEnds(self):
        mapper = Mapping("test")
        mapper.hit(1, 10, 1)
        mapper.hit(2, 19, 2)
        mapper.hit(25, 30, 2)
        returnList = mapper.getStartEnds()
        self.assertEqual(len(returnList), 3)
