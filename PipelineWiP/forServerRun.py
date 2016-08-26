#!/usr/bin/env python

#python forServerRun.py -Q AD84_S8_L001_R1_001P100.fastq -Q AD84_S8_L001_R2_001P100.fastq -g Phyca11_unmasked_genomic_scaffolds.fasta

import optparse
from Visualiser import VisualisationTools
from Main import Main
from ReadAligner import ReadAligner


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

options, args = parser.parse_args()

worker = Main.Main()
worker.makeDirectory(options.output_filepath)
worker.makeDirectory(options.result_filepath)
worker.openRefGenomes(options.refGenome)
worker.openFastQFiles(options.fastQFile)
if len(worker.fastQFileList) > 2:
    fastQPairs = len(worker.fastQFileList) - 1
    i = 0
    while i < fastQPairs:
        mapper = ReadAligner.Bowtie2(worker.fastQFileList[i], worker.fastQFileList[i+1],
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
    visualisationTool = VisualisationTools.VisualisationTools(bamWorker.samFile, options.result_filepath)
    visualisationTool.start()


