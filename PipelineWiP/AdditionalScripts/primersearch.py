# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 13:02:01 2016

@author: HC42845
"""

from tkinter import filedialog
import os

totalList = set()
truePos = set()
trueNeg = set()
falsePos = set()
falseNeg = set()

directory = filedialog.askdirectory()
for filename in os.listdir(directory):
    referenceFile = open(directory+"/"+filename)
    lengthList = list()
    POI = ""
    for line in referenceFile:
        if "Primer name" in line:
            line = line.split()
            POI = line[2].rstrip()
            totalList.add(POI)
            lengthList = list()
        elif "Amplimer length:" in line:
            line = line.split()
            lengthList.append(int(line[2]))
        elif line in ['\n', '\r\n'] and len(POI) > 0:
            lengthList = [x for x in lengthList if x >= 200 and x <= 500]
            if len(lengthList) > 0:
                truePos.add(POI)
            else:
                falseNeg.add(POI)
lengthList = [x for x in lengthList if x >= 200 and x <= 500]
if len(lengthList) > 0:
    truePos.add(POI)
else:
    falseNeg.add(POI)
            
directory = filedialog.askdirectory()
for filename in os.listdir(directory):
    otherOrg = open(directory+"/"+filename)
    lengthList = list()
    POI = ""
    for line in otherOrg:
        if "Primer name" in line:
            line = line.split()
            POI = line[2].rstrip()
            totalList.add(POI)
            lengthList = list()
        elif "Amplimer length:" in line:
            line = line.split()
            lengthList.append(int(line[2]))
        elif line in ['\n', '\r\n'] and len(POI) > 0:
            lengthList = [x for x in lengthList if x >= 50 and x <= 1000]
            if len(lengthList) > 0:
                falsePos.add(POI)
                if POI in trueNeg:
                    trueNeg.remove(POI)
            else:
                trueNeg.add(POI)
lengthList = [x for x in lengthList if x >= 50 and x <= 1000]
if len(lengthList) > 0:
    falsePos.add(POI)
    if POI in trueNeg:
        trueNeg.remove(POI)
else:
    trueNeg.add(POI)
print("Total: "+str(len(totalList)))
print("TruePos: "+str(len(truePos)))
print("TrueNeg: "+str(len(trueNeg)))
print("FalsPos: "+str(len(falsePos)))
print("FalseNeg: "+str(len(falseNeg)))