import matplotlib.pyplot as plt
from tkinter import filedialog

class scaffoldClass:
    
    def __init__(self, scaffold, start, end):
        self.scaffold = scaffold
        self.hitDict = dict()
        self.start = start
        self.end = end
    
    def hit(self, isolat, start, end):
        if isolat not in self.hitDict:
            self.hitDict[isolat] = list()
        if start < self.start:
            start = self.start
        if end > self.end:
            end = self.end
        self.hitDict[isolat].append([start, end])
        return
    
    def makePlot(self, number):
        print(self.scaffold)
        plt.ioff()
        colorList = ['g', 'g', 'c', 'c', 'm', 'm', 'b', 'b', 'darkred', 'darkred', 'orange', 'orange']
        #colorList = ['red','green','blue','indigo','gold','darkseagreen','indigo','gold','darkseagreen','lime','darkgray','mistyrose']
        #isolats = ["Genes","Intersect Unique", "Intersect", "Y006D", "Q108D", "AD84D", "Y006M", "Q108M", "AD84M"]
        isolats = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
        i = 0
        plt.figure(number)
        for isolat, value in iter(self.hitDict.items()):
            i = isolats.index(isolat)
            for startEnd in value:
                line = plt.plot([startEnd[0], startEnd[1]],[i+1,i+1], linewidth=4.0)
                plt.setp(line, color=colorList[i])
        plt.ylim(0,len(isolats)+1)
        plt.xlim(self.start, self.end)
        plt.ylabel("Isolat and unique regions P.capsici RxLR/CRN")
        plt.xlabel("Position (bp)")
        plt.yticks(range(1, len(isolats)+1),isolats)
        plt.grid(True)
        plt.title("Protein id: " + str(self.scaffold))
        place = "C:\\Users\\HC42845\\Desktop\\Graphs\\"+str(self.scaffold)+".svg"
        plt.savefig(place, format='svg', dpi=1200)
        plt.clf()
        return plt
        
scaffoldDict = dict()
#isolats = ["Genes","Intersect Unique", "Intersect", "Y006D", "Q108D", "AD84D", "Y006M", "Q108M", "AD84M"]
isolats = ["RxLR", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
for isolat in isolats:
    print(isolat)
    inputFile = filedialog.askopenfile()
    for line in inputFile:
        line = line.split("\t")
        if isolat == "Genes":# or isolat == "Intersect Unique" or isolat == "Mapper" or isolat == "Denovo":
            start = int(line[3])
            end = int(line[4])
        elif isolat == "RxLR":
            start = int(line[1])
            end = int(line[2])
        elif isolat == "Intersect":
            try:
                start = int(line[9])
                end = int(line[10])
            except ValueError:
                start = 0
                end = 0
        else: # isolat == "Mapper" or isolat == "Denovo":
            try:
                start = int(line[7])
                end = int(line[8])
            except ValueError:
                start = 0
                end = 0
#        else:
#            start = int(line[5])
#            end = int(line[6])
        scaffold = line[3].rstrip()
        if scaffold not in scaffoldDict:
            scaffoldDict[scaffold] = scaffoldClass(scaffold, start, end)
        if isolat is not "RxLR":
            scaffoldDict[scaffold].hit(isolat, start, end)
    
for key, value in iter(scaffoldDict.items()):
    value.makePlot(key)