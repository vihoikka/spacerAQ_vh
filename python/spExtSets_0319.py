#Ville Hoikkala 2020
#Finds spacers in ion torrent reads. Removes a given list of unwanted spacers. Tuned for F. columnare C1 (II-C) and C2 (VI-B).
#Commented lines for debugging purposes

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from subprocess import call #For command line tools
import sys
from fuzzywuzzy import fuzz
import time
from difflib import SequenceMatcher


print("-----------")
print("Starting SpacerExtractor")
userfile = sys.argv[1]
useroutfile = sys.argv[2]
output_filename = sys.argv[2]
locus = sys.argv[3]
duplicateRemoval = sys.argv[4]
discard = sys.argv[5]
print(output_filename)
output_handle1 = open(output_filename+".fasta", "w")
output_handle2 = open(output_filename+".fastq", "w")
discardable_spacers = {'C1':['TTAATTAAATCACAAAAAATGACAAGAGA','TGTTCTATCGGCGTAAAATAAGTTAGCAG'],'C2':['TATTCTAATGGATCAATTAAAAAAGCAAGA', 'TCGAGAGGTGTTAATGAAACGAACATAGTT','AATTACTGCTAAATGATACAATGTGTTTG']}
discardThreshold = 90 #When removing discardable spacers from results, use this threshold (0-100)

readsList = list(SeqIO.parse(userfile, "fastq")) #Read fastq file into a list

def removeFuzzyDuplicates(unFuzziedList):
    unFuzzied = set()
    for i in unFuzziedList:
        unFuzzied.add(i.seq)

    cycleCountDuplicate = 0
    progressCountDuplicate = 0
    tenPercentDuplicate = len(unFuzzied)*0.01

    fuzziness = 90
    print("Removing similar spacers using fuzzyness threshold of " + str(fuzziness))
    uniqueSpacers = set()
    simiCounter = 0
    for sequence in unFuzzied:
        isSimilar = False
        #if len(uniqueSpacers) > 0:
        for uniqueSequence in uniqueSpacers:
            #print ("Looking at similarity of " + str(value.seq))
            #p = Pool()
            #result = p.map(fuzz.ratio(sequence.seq, uniqueSequence.seq))
            similarity = fuzz.QRatio(sequence, uniqueSequence) #This is done to remove any remaining duplicates using a fuzzy search
            #print (similarity)
            #print("Similarity between")
            #print(sequence.seq)
            #print(uniqueSequence.seq)
            #print(str(similarity))
            if similarity > fuzziness:
                isSimilar = True
                simiCounter = simiCounter + 1
                break
        if isSimilar == False:
            #print ("Found unique spacer")
            uniqueSpacers.add(sequence)
        #else:
            #uniqueSpacers.add(sequence)
        #time.sleep(0.5)
        cycleCountDuplicate = cycleCountDuplicate + 1

        if cycleCountDuplicate > tenPercentDuplicate:
            cycleCountDuplicate = 0
            progressCountDuplicate = progressCountDuplicate + 1
            print ("   " + str(progressCountDuplicate) + " % of spacers fuzzy-checked. Number of spacers:" + str(len(unFuzzied)) + ", skipped due to duplicicity: " + str(simiCounter) + ", unique sequences: " + str(len(uniqueSpacers)) + "\r", end="")

    return uniqueSpacers

def lookForSpacers(locus, readsList):
    #Define repeats
    startOfRepeats = {"C2":"GTTG","C1":"AAT"}
    endOfRepeats = {"C2":"ACAAC","C1":"CAAC"}
    searchBound1sts = {"C2":115,"C1":108}
    searchBound2nds = {"C2":40,"C1":40}

    #Define searchbounds
    startOfRepeat = startOfRepeats[locus] #First bases of repeat sequence
    endOfRepeat = endOfRepeats[locus] #Last bases of repeat sequence
    searchBound1st = searchBound1sts[locus] #Max coordinate for searching start of repeat. For non-trimmed ~60, trimmed 40.
    searchBound2nd = searchBound2nds[locus]

    print ("Looking for spacers in " + locus + ". Found " + str(len(readsList)) + " reads.")

    #Setup lists for spacers
    spacersList = [] #Initiate spacer list
    skipCounter = 0
    finalSpacers = []

    #Setup progress counter
    cycleCount = 0
    progressCount = 0
    percentOfTotal = len(readsList)*0.01

    for read in readsList:
        #print ("Read: " + read.seq)
        #print ("Searching for end of first repeat: " + endOfRepeat + " between 0 and " + str(searchBound1st))
        endRep = read.seq.find(endOfRepeat,0,searchBound1st) #Search for the end of first repeat
        #print ("Starting coordinate of spacer: " + str(endRep+len(endOfRepeat)))
        startCut = read[endRep+len(endOfRepeat):] #Cut all before the first spacer. This also cuts the corresponding phred scores
        #print ("Starting with spacer: " + startCut.seq)
        startRep = startCut.seq.find(startOfRepeat,26,searchBound2nd) #same for second repeat
        #print("Startrep: " + str(startRep))
        #if startRep == -1: #if the next repeat was not found, assume spacer length of 30
            #print ("Did not find second repeat, assuming spacer length 30")
            #startRep = 30
        sequence = startCut[:startRep]
        #print ("FINAL SPACER: " + sequence.seq)
        #time.sleep(2)
        #print sequence.format("fastq")
        #print ("Extraction length: " + str(len(sequence.seq)))
        if len(sequence.seq) < 34 and len(sequence.seq) > 27 and startRep != -1: #check length
            spacersList.append(sequence) #Copy the entry into the new list (as a sequence object)
        else:
            skipCounter = skipCounter + 1
                        #print ("Duplicate sequence found")
                        #time.sleep(0.2)

        #for key2, value2 in dupCheckList.items():
            #print (value2)
            #time.sleep(2)
            #print ("Spacer added to list!")

        cycleCount = cycleCount + 1
        if cycleCount > percentOfTotal:
            cycleCount = 0
            progressCount = progressCount + 1
            print ("   " + str(progressCount) + " % of spacers checked. Number of spacers:" + str(len(spacersList)) + ", skipped due to length: " + str(skipCounter) + "\r", end="")


    if duplicateRemoval == "dup":
        counter = 0
        print ("Checking for duplicates...")
        readSet = set()
        unFuzziedList = []
        for read in spacersList: #convert list to a set, which removes any 100% duplicates fast
            readSet.add(read.seq)
            counter = counter + 1
        uniqueReadsList = list(readSet)
        for read in uniqueReadsList:
            readObject = SeqRecord(read, id = str(counter)) #Convert the set back to SeqRecord objects for easier writing using SeqIO.write (see below)
            unFuzziedList.append(readObject)
        finalSpacers = removeFuzzyDuplicates(unFuzziedList)
        #finalSpacers = removeDuplicatesSequenceMatcher(unFuzziedList)

        noOfUniques = len(readSet)
        noOfRemoved = len(spacersList) - noOfUniques
        print ("Original data contained " + str(len(spacersList)) + " spacers; " + str(noOfRemoved) + " were removed as duplicates; continuing with " + str(noOfUniques) + " spacers")
    else:
        print ("Duplicates not removed")
        finalSpacers = spacersList

    if locus == "C1": #II-C spacers are flipped to reflect transcribed spacers
        newList = []
        print("Reverse-complementing C1 spacers")
        for entry in finalSpacers: #Reverse complement all reads in C1
            #print (entry.seq)
            entry.seq = entry.seq.reverse_complement()
            newList.append(entry)
            #print (entry)
            #print ("---")
        finalSpacers = newList

    return finalSpacers

def removeDiscardableSpacers(spacerlist,locus,discardDict,fuzziness):
    print("Removing discardable spacers as requested " + "(" + locus + ")")
    newList = []
    simiCounter = 0
    for sequence in spacerlist:
        isSimilar = False
        for oldspacer in discardDict[locus]:
            similarity = fuzz.ratio(sequence.seq, oldspacer) #Uses fuzzywuzzy to compare spacer to a discardable spacer
            if similarity > fuzziness:
                #print(sequence.seq + " vs " + oldspacer)
                #print(str(similarity))
                #time.sleep(3)
                isSimilar = True
                simiCounter = simiCounter + 1
                break
        if isSimilar == False:
            #print ("Found unique spacer")
            newList.append(sequence)
    print ("    Discarded " + str(simiCounter) + " spacers")
    return newList

def writeToDisk(finalSpacers):
    SeqIO.write(finalSpacers, output_handle1, "fasta")
    output_handle1.close()
    SeqIO.write(finalSpacers, output_handle2, "fastq")
    output_handle2.close()

'''MAIN PROGRAM'''

finalSpacers = lookForSpacers(locus, readsList)
if discard == 'discard':
    finalSpacers = removeDiscardableSpacers(finalSpacers,locus,discardable_spacers,discardThreshold)
print("Total number of spacers: " + str(len(finalSpacers)))
print("SpacerExtractor finished")
writeToDisk(finalSpacers)

print("-----------")
