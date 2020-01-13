'''
Extracts PAMs from an input SAM and reference genome file
@author: Ville Hoikkala
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import csv

class spacer:
    start = 0
    stop = 0
    ID = ""
    length = 0
    strand = ""
    pamDown = ""
    pamUp = ""

def getSpacerPositions(file,direction):
    parsed = list(csv.reader(open(file, 'r'), delimiter='\t')) #opens file as tab-delimited
    numberOfEntries = len(parsed)-3 #Count the number of entries. Sam files contain three rows of non-entry data
    #spacerList = [spacer() for i in range(numberOfEntries)] #Create list to store spacers as objects
    spacerList = [] #Create list to store spacers as objects
    lineCount = 0 #for skipping the header lines of the SAM file
    spacerCount = 0 #keep count of spacers
    if direction == "R":
        for line in parsed:
            #print(lineCount)
            if lineCount > 2: #Skip the headers
                if int(line[1]) == 16: #0 for F hits, 16 for reverse hits
                    #print(lineCount)
                    #print("Name of spacer: " + line[0])
                    #print("Spacer start pos: " + str(line[3]))
                    ID = line[0] #name of query
                    length = int(len((line[9]))) #spacer sequence is stored in column 9
                    start = int(line[3]) #protospacer start is in column 3
                    stop = int(line[3]) + length #we can deduce the stop by knowing the start and length
                    spacerList.insert(spacerCount, spacer()) #insert a new spacer object into the list
                    setattr(spacerList[spacerCount], 'start', start)
                    setattr(spacerList[spacerCount], 'stop', stop)
                    setattr(spacerList[spacerCount], 'ID', ID)
                    setattr(spacerList[spacerCount], 'length', length)
                    setattr(spacerList[spacerCount], 'strand', 'R')
                    spacerCount = spacerCount + 1
            lineCount = lineCount + 1
    if direction == "F":
        for line in parsed:
            #print(lineCount)
            if lineCount > 2: #Skip the headers
                if int(line[1]) == 0: #0 for F hits, 16 for reverse hits
                    ID = line[0] #name of query
                    length = int(len((line[9]))) #spacer sequence is stored in column 9
                    start = int(line[3]) #protospacer start is in column 3
                    if start < 15:
                        break
                    stop = int(line[3]) + length #we can deduce the stop by knowing the start and length
                    spacerList.insert(spacerCount, spacer()) #insert a new spacer object into the list
                    setattr(spacerList[spacerCount], 'start', start)
                    setattr(spacerList[spacerCount], 'stop', stop)
                    setattr(spacerList[spacerCount], 'ID', ID)
                    setattr(spacerList[spacerCount], 'length', length)
                    setattr(spacerList[spacerCount], 'strand', 'F')
                    #print("Name of spacer: " + line[0])
                    #print("Spacer start pos: " + str(start))
                #if direction == "R" and int(line[1]) == 16: #0 for F hits, 16 for reverse hits
                #   posDic[line[0]] = int(line[3]) #"Position" in each line on column 3
                #print line[3]
                    spacerCount = spacerCount + 1
            lineCount = lineCount + 1
    return spacerList

def removeOverlaps(spacers,limit):
    dir = spacers[0].strand
    print("Removing overlapping spacers from " + dir + " hits")
    uniqueSpacers = []

    simiCounter = 0 #For counting overlapping protospacers
    for s in spacers:
        overlap = False
        for uniS in uniqueSpacers:
            distance = abs(s.start - uniS.start) #Calculate distance between protospacers
            #print ("Distance: " + str(distance))
            if distance < limit: #If limit is exceeded...
                overlap = True
                simiCounter = simiCounter + 1
                break #... break the loop
        if overlap == False: #If no overlap was found...
            uniqueSpacers.append(s)

    print ("    Number of spacers: " + str(len(spacers)) + ", of these skipped due to overlap: " + str(simiCounter) + "\r", end="")
    print("")
    return uniqueSpacers

def pamSeek(spacers, genome, pamSize,dir):
    #print (genome.id)
    if dir == "F":
        for i in spacers:
            upPamStart = i.start-pamSize-1 #search for start position of PAM
            #print("Start pos of upPAM: " + str(upPamStart))
            upPamStop = i.start-1 #search for stop pos of PAM
            upPam = genome.seq[upPamStart:upPamStop] #PAM is area between these
            setattr(i, 'pamUp', upPam) #set as attribute for the spacer

            #print(i.ID + "Upstream PAM: " + upPam)

            downPamStart = i.stop-1
            downPamStop = i.stop+pamSize-1
            downPam = genome.seq[downPamStart:downPamStop]
            setattr(i, 'pamDown', downPam)

            #print(i.ID + "Downstream PAM: " + downPam)


        print("Number of F spacers: " + str(len(spacers)))
        return spacers
    if dir == "R":
        for i in spacers:
            upPamStart = i.stop-1 #search for reversed start position of PAM
            #print("Start pos of upPAM: " + str(upPamStart))
            upPamStop = i.stop+pamSize-1 #search for stop pos of PAM
            upPam = SeqRecord(genome.seq[upPamStart:upPamStop]) #PAM is area between these
            upPam = upPam.reverse_complement().seq
            setattr(i, 'pamUp', upPam) #set as attribute for the spacer

            #print(i.ID + "Upstream PAM: " + upPam.seq)

            downPamStart = i.start-pamSize-1
            downPamStop = i.start-1
            downPam = SeqRecord(genome.seq[downPamStart:downPamStop])
            downPam = downPam.reverse_complement().seq
            setattr(i, 'pamDown', downPam)

            #print(i.ID + "Downstream PAM: " + downPam.seq)


        print("Number of R spacers: " + str(len(spacers)))
        return spacers


""" *** PROGRAM STARTS HERE *** """
print("Starting PAMseek")
upStreamPAMsOrigin = [] #List for all spacers where the PAMs come from
downStreamPAMsOrigin = [] #And downstream
totalHits = 0

inputFile = sys.argv[1] #1st argument: Bowtie2 generated SAM-file
trgt_name = sys.argv[2] #genome to be mapped on (fasta file in the same folder)
outputFilename = sys.argv[3] #output file which contains PAM sequences
overlapLimit = int(sys.argv[4]) #Maximum allowed overlap. Zero disables (retains all spacers)
pamSize = int(sys.argv[5]) #Size of PAM region

trgt = trgt_name + ".fasta" #add extension to target genome

#Define output filenames
upOutfile = outputFilename+"_PAMs_upstream_"+trgt
downOutfile = outputFilename+"_PAMs_downstream_"+trgt

#Turn the genome into a SeqRecord object
genome = SeqIO.read(trgt, "fasta")

'''Obtaining positions on the target genome'''

#First, create spacer objects from each protospacer and make list
spacersF = getSpacerPositions(inputFile,"F") #This list contains spacer objects that hit the phage in F strand
spacersR = getSpacerPositions(inputFile,"R") #Same for other strand

#Second, feed the spacer objects into pamSeek function. This returns the set of spacers updated with PAM information
spacersF = pamSeek(spacersF, genome, pamSize,"F")
spacersR = pamSeek(spacersR, genome, pamSize,"R")

#Third, remove overlapping spacers to get unique PAMs
if overlapLimit > 0:
    spacersF = removeOverlaps(spacersF,overlapLimit)
    spacersR = removeOverlaps(spacersR,overlapLimit)

#Fourth, get PAMs from the spacer objects and add to lists
upstreamPAMs = []
downstreamPAMs = []
print("Fetching PAMs from target genome...")
for i in spacersF: #Fetch PAM sequences from SeqRecords (spacers) and add them to lists
    upstreamPAMs.append(SeqRecord(Seq(str(i.pamUp)),id=i.ID, description="upstream PAM"))
    downstreamPAMs.append(SeqRecord(Seq(str(i.pamDown)), id=i.ID, description="downstream PAM"))

for i in spacersR:
    upstreamPAMs.append(SeqRecord(Seq(str(i.pamUp)),id=i.ID, description="upstream PAM"))
    downstreamPAMs.append(SeqRecord(Seq(str(i.pamDown)), id=i.ID, description="downstream PAM"))

print("Writing PAM fasta files")

#Write the files
SeqIO.write(upstreamPAMs, upOutfile, "fasta")
SeqIO.write(downstreamPAMs, downOutfile, "fasta")

print("PAMseek finished!")
print("------------------")
