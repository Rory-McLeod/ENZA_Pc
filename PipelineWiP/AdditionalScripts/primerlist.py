# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 09:19:23 2016

@author: HC42845
"""

from tkinter import filedialog

primerFile = filedialog.askopenfile()
outputfile = filedialog.asksaveasfile()

getPrimers = False
counter = 100
resultList = list()
for line in primerFile:
    if "PRIMER PICKING RESULTS FOR" in line:
        line = line.split(" ")
        POI = line[4].rstrip()
    elif "No mispriming library specified" in line:
        counter = 0
        getPrimers = True
    elif getPrimers == True and counter == 3:
        if "LEFT PRIMER" in line:
            line = line.split()
            Lprimer = line[9].rstrip()
        else:
            getPrimers = False
    elif getPrimers == True and counter == 4:
        if "RIGHT PRIMER" in line:
            line = line.split()
            Rprimer = line[9].rstrip()
            outputfile.write(POI + "\t" + Lprimer + "\t" + Rprimer + "\n")
    counter += 1
primerFile.close()
outputfile.close()