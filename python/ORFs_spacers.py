'''
First created Mar 22, 2016

Takes a genome and spacer as an input.
1. Look for a protospacer.
2. If found, store
start (S) and end (E) coordinates. If not found, reverse complement pattern
and try again (else not found).
3. From the S coordinate upstream store a sequence of length t in a list.
4. From the E coordinate store a sequence of length t in a list.

@author: Ville Hoikkala 2020
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import csv
from decimal import Decimal
from math import factorial

class spacer:
    start = 0
    stop = 0
    id = ""
    length = 0
    strand = ""
    mRNAcapable = 0
    ORF = ""
    ORF_description = ""
    intergenic = "yes"

class ORF:
    def __init__(self):
        self.start = 0
        self.stop = 0
        self.ID = ""
        self.strand = ""
        self.length= 0
        self.hits = 0
        self.mRNAhits = 0
        self.totalHPN = 0.0
        self.HPNF = 0.0
        self.HPNR = 0.0
        self.mRNAHPN = 0.0
        self.mRNAhitRatio = 0.0
        self.mRNAhitRatioUnique = 0.0
        self.annotation = ""
        self.P_value = 1
        self.uniqueFs = []
        self.uniqueRs = []
        self.uniquemRNAhits = []

def getSpacerPositions(file):
    parsed = list(csv.reader(open(file, 'r'), delimiter='\t')) #opens file as tab-delimited
    numberOfEntries = len(parsed)-3 #Count the number of entries. Sam files contain three rows of non-entry data
    print ("Number of spacers in SAM file: " + str(numberOfEntries))
    #spacerList = [spacer() for i in range(numberOfEntries)] #Create list to store spacers as objects
    spacerList = [] #Create list to store spacers as objects
    lineCount = 0 #for skipping the header lines of the SAM file
    spacerCount = 0 #keep count of spacers
    for line in parsed:
        if lineCount > 2: #Skip the headers
            strand = int(line[1])
            if ((strand == 16) or (strand == 0)): #if the spacers hits the genome
                ID = line[0] #name of query
                length = int(len((line[9]))) #spacer sequence is stored in column 9
                start = int(line[3]) #protospacer start is in column 3
                stop = int(line[3]) + length #we can deduce the stop by knowing the start and length
                spacerList.insert(spacerCount, spacer()) #insert a new spacer object into the list
                setattr(spacerList[spacerCount], 'start', int(start))
                setattr(spacerList[spacerCount], 'stop', int(stop))
                setattr(spacerList[spacerCount], 'id', ID)
                setattr(spacerList[spacerCount], 'length', length)
                if strand == 0:
                    strand = "f"
                if strand == 16:
                    strand = "r"
                setattr(spacerList[spacerCount], 'strand', strand)
                spacerCount = spacerCount + 1
        lineCount = lineCount + 1
    return spacerList

def getORFs(csvfile):
    file = csvfile
    rowCounter = 0
    ORFcount = 0
    ORFlist = []
    for row in file: #Create an ORF object from each row
        if rowCounter > 0: #Skip header row
            id = row[0]
            annotation = row[1]
            start = row[2]
            stop = row[3]
            length = row[4]
            direction = row[5]
            ORFlist.insert(ORFcount, ORF())
            setattr(ORFlist[ORFcount], "id", id)
            setattr(ORFlist[ORFcount], "start", int(start))
            setattr(ORFlist[ORFcount], "stop", int(stop))
            setattr(ORFlist[ORFcount], "strand", str(direction))
            setattr(ORFlist[ORFcount], "length", int(length))
            setattr(ORFlist[ORFcount], "annotation", annotation)
            ORFcount = ORFcount + 1
        rowCounter = rowCounter + 1
    print ("Number of ORFs: " + str(len(ORFlist)))
    return ORFlist


""" *** PROGRAM STARTS HERE *** """
print("-------------------")
print("Starting ORFmapper")
upStreamPAMsOrigin = [] #List for all spacers where the PAMs come from
downStreamPAMsOrigin = [] #And downstream
totalHits = 0

inputFile = sys.argv[1] #1st argument: Bowtie2 generated SAM-file
print (inputFile)
ORFfile = sys.argv[2] #.csv file that contains ORF id, start, stop and direction
outputFilename = sys.argv[3] #output file
upOutfileO = outputFilename+"_ORF_stats.csv" #actual output filename for ORF stats
upOutfileS = outputFilename+"_ORF_spacer_stats.csv" #actual output filename for spacer stats

with open(ORFfile, 'rt') as csvfile:
    ORFreader = csv.reader(csvfile, delimiter = ";")
    ORFs = list(ORFreader)

ORFobjects = getORFs(ORFs) #Make a list with all ORF objects


#Create spacer objects from each protospacer and make list
spacers = getSpacerPositions(inputFile) #This list contains spacer objects that hit the genome on any strand

for s in spacers: #Go through all spacers
    for o in ORFobjects: #and check each ORF against each spacer
        if (s.start >= o.start and s.stop <= o.stop): #if spacer is inside ORF
            if (s.strand == "f" and s.start not in o.uniqueFs):
                o.uniqueFs.append(s.start) #add starting positions to list of unique spacers of this ORF
                #print("Adding start pos " + str(s.start) + " to uniqueFs, length now " + str(len(o.uniqueFs)))
            if (s.strand == "r" and s.start not in o.uniqueRs):
                o.uniqueRs.append(s.start)
                #print("Adding to Rs, length now " + str(len(o.uniqueRs)))
            o.hits = o.hits + 1 #add one to ORF's 'hits' attribute
            setattr(s, "ORF", o.id) #add id of ORF to spacer's 'ORF' attribute
            setattr(s, "ORF_description", o.annotation) #add id of ORF to spacer's 'ORF' attribute
            setattr(s, "intergenic", "no")
            if (s.strand != o.strand): #if they are on separate strands, then the spacer is mRNA capable
                o.mRNAhits = o.mRNAhits + 1 #add one to ORF's 'mRNAhits' attribute
                setattr(s, "mRNAcapable", 1) #label the spacer mRNA capable
                if (s.start not in o.uniquemRNAhits):
                    o.uniquemRNAhits.append(s.start) #add starting positions to list of unique spacers of this ORF

#Calculate HPN (normalizing number of spacers hits against the length of the ORF)
for o in ORFobjects:
    totalUniqueHits = len(o.uniqueFs) + len(o.uniqueRs)
    hpn = round(Decimal(totalUniqueHits/o.length),4)
    setattr(o, "HPN", hpn)

#Calculate HPN for mRNA targeting spacers
for o in ORFobjects:
    hpnmRNA = round(Decimal(len(o.uniquemRNAhits)/o.length),4)
    setattr(o, "mRNAHPN", hpnmRNA)

for o in ORFobjects:
    HPNF = len(o.uniqueFs)/o.length
    setattr(o, "HPNF", HPNF)
    HPNR = len(o.uniqueRs)/o.length
    setattr(o, "HPNR", HPNR)

#Calculate the ratio of mRNA targeting vs all spacers (absolute)
for o in ORFobjects:
    if o.hits > 0:
        mRNAper = round(o.mRNAhits/o.hits,3)
        setattr(o, "mRNAhitRatio", mRNAper)

#Calculate the ratio of mRNA targeting vs all sacers (unique)
for o in ORFobjects:
    if o.hits > 0:
        totalUniqueHits = len(o.uniqueFs) + len(o.uniqueRs)
        mRNAratio = round(len(o.uniquemRNAhits)/totalUniqueHits,3)
        setattr(o, "mRNAhitRatioUnique", mRNAratio)


#Binomial tests if an ORF has bias on spacers on either strand and gives a P value. Based on unique spacers.
for o in ORFobjects:
    n = int(len(o.uniqueFs+o.uniqueRs))
    x = int(len(o.uniquemRNAhits))
    p = float(0.5)
    q = float(0.5)
    if x > 0:
        P = ((factorial(n))/((factorial(n-x)*factorial(x)))) * (p**x) * (q**(n-x))
        setattr(o, "P_value", round(P,5))

#Sort the ORF list
ORFobjects.sort(key=lambda x: x.P_value, reverse=False)
#for i in ORFobjects:
 #   print (i.id + " (" + i.annotation + ") " + ". Hits: " + str((i.hits)) + ". HPN: " + str(i.HPN) + ", mRNAHPN: " + str(i.mRNAHPN) + ", mRNA hit ratio: " + str(i.mRNAhitRatio))

#Output csv file with the following header row
outList = []
outTable = [["ORF id","annotation","Hits","Unique F hits","Unique R hits","mRNA hits", "Unique mRNA hits","Normalized uniq","Normalized F","Normalized R","Normalized (mRNA) uniq","mRNA ratio uniq", "mRNA ratio abs","Length","Start","Stop","P value uniq"]]

#Print each ORF into the csv file
for o in ORFobjects:
    values = [o.id,o.annotation,str(o.hits),str(len(o.uniqueFs)),str(len(o.uniqueRs)),str(o.mRNAhits), str(len(o.uniquemRNAhits)),str(o.HPN),str(o.HPNF),str(o.HPNR),str(o.mRNAHPN),str(o.mRNAhitRatioUnique),str(o.mRNAhitRatio), str(o.length),str(o.start),str(o.stop),str(o.P_value)]
    outTable.append(values)

with open(upOutfileO, "w") as outfile:
    writer = csv.writer(outfile)
    writer.writerows(outTable)

#Also some stats for the console:
mRNACapableCount = 0
for s in spacers:
    if s.mRNAcapable == 1:
        mRNACapableCount = mRNACapableCount + 1
mRNApercent = round(mRNACapableCount/len(spacers)*100,2)

F_spacers = len(list(s for s in spacers if s.strand == "f"))
F_spacersPer = round(F_spacers/len(spacers)*100,2)
R_spacers = len(list(s for s in spacers if s.strand == "r"))
R_spacersPer = round(R_spacers/len(spacers)*100,2)

print ("Number of targeting spacers: " + str(len(spacers)))
print ("Mapped to F strand: " + str(F_spacers) + " (" + str(F_spacersPer) + ")")
print ("Mapped to R strand: " + str(R_spacers) + " (" + str(R_spacersPer) + ")")
print ("mRNAcapable spacers: " + str(mRNACapableCount)+ " (" + str(mRNApercent) + "%)")

spacerTable = [["spacer_id","strand","ORF_id","ORF_description","mRNA_capable","intergenic","length","start"]]

for s in spacers:
    spacerStats = [s.id,s.strand,s.ORF,s.ORF_description,s.mRNAcapable,s.intergenic,s.length,s.start]
    spacerTable.append(spacerStats)

with open(upOutfileS, "w") as outfileS:
    writer = csv.writer(outfileS)
    writer.writerows(spacerTable)


print ("----------------")
