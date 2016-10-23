#!/usr/bin/env python
"""
Todo:
- remove execute from classes and pass to Main class
- add check for GFF file
"""

import optparse
from VisualisationTools import VisualisationTools
import Assemblers
import Main
import ReadAligner
import threading
import PrimerDesign
from Blast import Blast
import copy
import NUCmer


def MapperPreRun():
    if len(worker.fastQFileList) > 2:
        fastQPairs = len(worker.fastQFileList) - 1
        i = 0
        while i < fastQPairs:
            mapper = ReadAligner.Bowtie2(worker.fastQFileList[i], worker.fastQFileList[i + 1],
                                         worker.refGenomeList[0], options.output_filepath)
            worker.mapperClass.append(mapper)
            mapper.start()
            i += 2
    else:
        mapper = ReadAligner.Bowtie2(worker.fastQFileList[0], worker.fastQFileList[1],
                                     worker.refGenomeList[0], options.output_filepath)
        worker.mapperClass.append(mapper)
        mapper.start()
    for mapper in worker.mapperClass:
        mapper.join()
    for mapper in worker.mapperClass:
        bamWorker = ReadAligner.BamTools(mapper.samFile, mapper.referenceDB)
        bamWorker.start()
        worker.bamClass.append(bamWorker)

    for bamWorker in worker.bamClass:
        bamWorker.join()

    for bamWorker in worker.bamClass:
        visualisationTool = VisualisationTools(bamWorker.samFile)
        visualisationTool.start()
        worker.visualisationClass.append(visualisationTool)

    for visualisationTool in worker.visualisationClass:
        visualisationTool.join()

    mapperPrimer = PrimerDesign.PrimerDesignByMapping()
    for visualisationTool in worker.visualisationClass:
        mapperPrimer.generateCoords(visualisationTool.depthPerPos)
    return PrimerDesign.PrimerDesign.runIntersect(mapperPrimer.coordsFile, "/MapperPoI.gff")

def DeNovoPreRun():
    if len(worker.fastQFileList) > 2:
        fastQPairs = len(worker.fastQFileList) - 1
        i = 0
        while i < fastQPairs:
            assembler = Assemblers.Spades(worker.fastQFileList[i], worker.fastQFileList[i + 1], options.output_filepath)
            worker.assemblerClass.append(assembler)
            assembler.start()
            i += 2
    else:
        assembler = Assemblers.Spades(worker.fastQFileList[0], worker.fastQFileList[1], options.output_filepath)
        worker.assemblerClass.append(assembler)
        assembler.start()
    contigs = ""
    for assembler in worker.assemblerClass:
        assembler.join()
        contigs += assembler.outputDir + " "
    thread = threading.Thread(Assemblers.Assemblers.quast(contigs))
    thread.start()
    Main.Main.threadList.append(thread)
    contigs = contigs.split(" ")
    nucmerList = list()
    for contig in contigs:
        if len(contig) > 0:
            nucmerRun = NUCmer.NUCmerRun(contig)
            nucmerRun.start()
            nucmerList.append(nucmerRun)
    for nucmerRun in nucmerList:
        nucmerRun.join()
    denovoPrimer = PrimerDesign.PrimerDesignByDenovo()
    for contig in contigs:
        print contig
        if len(contig) > 0:
            denovoPrimer.readCoords(contig)
    return PrimerDesign.PrimerDesign.runIntersect(denovoPrimer.coordsFile, "/denovoPoI.gff")

def intersect():

    return

fastQFileList = []
refGenomeList = []

parser = optparse.OptionParser()
parser.add_option('-o', '--output',
                  dest="output_filepath",
                  default="workDir",
                  help="Output directory for all the subfiles like sam, bam and fasta files"
                  )
parser.add_option('-r', '--resultOutput',
                  dest="result_filepath",
                  default="resultDir",
                  help="Result directory for all the results like coverage plots and primer outputs"
                  )
parser.add_option('-Q', '--fastQ',
                  dest="fastQFile",
                  action="append",
                  help="All the fastQ files (if pairwise data, input first the first file, and then the second file)"
                  )
parser.add_option('-q', '--fastQMethod',
                  dest="fastQMethod",
                  action="append",
                  default=0,
                  help="0 for pairwise fastQ data, other methods not yet made"
                  )
parser.add_option('-g', '--referenceGenome',
                  dest="refGenome",
                  action="append",
                  help="Reference genomes. each of this genome will be used in the program, only fasta format!"
                  )
parser.add_option('-m', '--mappingMethod',
                  dest="mapping",
                  default=0,
                  help="Mapping method, 0 for Bowtie2, other methods are coming"
                  )
parser.add_option('-p', '--primerDesign',
                  dest="primer",
                  default=0,
                  action="append",
                  help="attributes for primerDesign"
                  )
parser.add_option('-a', '--assemblerMethod',
                  dest="assembler",
                  default=0,
                  help="Assembler methods, 0 for SPADes, 1 for Velvet, 2 for MiRa"
                  )
parser.add_option('-G', '--GFFFile',
                  dest="gffFile",
                  help="GFF file. Make sure the chromosomes/scaffold have the same name as the reference genome!")

options, args = parser.parse_args()


Main.Main.makeDirectory(options.output_filepath)
Main.Main.makeDirectory(options.result_filepath)
Main.Main.gffFile = options.gffFile
Main.Main.workDir = options.output_filepath
Main.Main.resultDir = options.result_filepath
Main.Main.threadList = list()
worker = Main.Main()
worker.openRefGenomes(options.refGenome)
worker.openFastQFiles(options.fastQFile)
threadList = list()
for fastQFile in worker.fastQFileList:
    workLine = "fastqc " + Main.Main.fastQAdd + fastQFile + " -o " + Main.Main.resultDir + " -q --noextract"
    thread = threading.Thread(Main.Main.execute(workLine, "Generating fastQC reports in the background"))
    thread.start()
    Main.Main.threadList.append(thread)

preRun = list()
methodList = list()
thread = threading.Thread(methodList.append(MapperPreRun()))
thread.start()
preRun.append(thread)
thread = threading.Thread(methodList.append(DeNovoPreRun()))
thread.start()
preRun.append(thread)
for thread in preRun:
    thread.join()

methodList.append(PrimerDesign.PrimerDesign.runMethodIntersect(methodList, Main.Main.workDir+"/Intersect.gff"))
PrimerDesign.PrimerDesign.removeDuplicate(methodList[2])
methodList.append(PrimerDesign.PrimerDesign.runMethodSubstract([methodList[0], methodList[2]], Main.Main.workDir +
                                                               "/MapperUnique.gff"))
methodList.append(PrimerDesign.PrimerDesign.runMethodSubstract([methodList[1], methodList[2]], Main.Main.workDir +
                                                               "/DenovoUnique.gff"))
genome = PrimerDesign.PrimerDesign.readRefGenome(Main.Main.genomeAdd + Main.Main.refGenomeList[0])
blastList = list()
for PipMethod in methodList:
    fastaFile = PipMethod.rstrip()[:-3]+"fa"
    PrimerDesign.PrimerDesign.saveFasta(fastaFile, PrimerDesign.PrimerDesign.readGFF(PipMethod, genome))
    blastList.append(fastaFile)

otherGenomes = copy.copy(Main.Main.refGenomeList)
del otherGenomes[0]
blastThread = list()
for blastItem in blastList:
    print blastItem
    blastResult = Blast(blastItem, otherGenomes)
    blastResult.start()
    blastThread.append(blastResult)

for blastItem in blastThread:
    blastItem.join()

for k, blastItem in enumerate(blastList):
    item = blastItem.rstrip()[:-2]
    thread = threading.Thread(PrimerDesign.PrimerDesign.generatePrimer3Input(
        item+"primSets", PrimerDesign.PrimerDesign.readGFF(item+"unique.gff", genome)))
    thread.start()
    Main.Main.threadList.append(thread)

for thread in Main.Main.threadList:
    thread.join()
