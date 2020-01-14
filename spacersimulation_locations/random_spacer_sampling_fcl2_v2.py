from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse
import random


parser = argparse.ArgumentParser(description='Generate random spacers across genome') #Initiates parser

parser.add_argument('-g', '--genome', help='input genome', type=str) #Adds genome argument
parser.add_argument('-m', '--mode', help='sample randomly or from PAMs (random or PAM)', type=str) #Adds genome argument
parser.add_argument('-c', '--count', help='number of spacers to be sampled', type=int) #Adds genome argument
parser.add_argument('-ps', '--PAMseq', help='PAM sequence', type=str) #Adds genome argument
parser.add_argument('-pp', '--PAMpolarity', help='Is the PAM downstream or upstream of the protospacer', type=str) #Adds genome argument


def createPAMspacers(g,c,p,pp):
    PAMlength = len(p)
    PAMlocationDic = {}
    strands = ["f","r"]
    strandCount = {}
    NCounter = 0
    while (p.find("N") != -1): #Prune PAM to get rid of NNNN
    #    print ("Splicing an N from the start of PAM")
        p = p[1:len(p)]
        NCounter = NCounter + 1
    #    print ("NCounter " + str(NCounter))

    for s in strands: #Go through both strands
        pamCount = 0
        ntCount = 0 #Nucleotide counter for progressing through genome
        if (s == "r"): #Flip the PAM when doing reverse complement
            p = p.reverse_complement()
        #    print ("Flipped PAM to revcomp (" + p + ")")
        for nt in g.seq: #Go through genome one nucleotide at a time
            if (g.seq[ntCount:ntCount+len(p)] == p): #If a sequence identical to PAM is found...
                PAMlocationDic[str(ntCount)] = s #Store its position and strand in the dictionary. Position as key, strand as value
                pamCount = pamCount + 1
            ntCount = ntCount + 1
        strandCount[s] = pamCount
#    print (str(strandCount["f"]) + " F PAMs and " + str(strandCount["r"]) + " R PAMs. Picking " + str(c) + " random PAMs")
    randomPAMlocations = random.sample(set(list(PAMlocationDic)), c) #random() requires a list, so we must break down the dictionary structure
    randomPAMstrands = [PAMlocationDic[k] for k in randomPAMlocations]
    randomPAMdict = {}
    newDicCycler = 0
    for i in randomPAMlocations: #Recompiling the dictionary
        randomPAMdict[i] = randomPAMstrands[newDicCycler]
        newDicCycler = newDicCycler + 1
#    print(len(randomPAMdict))

    spacerCounter = 1
    spacerList = []
    for key in randomPAMdict:
        if pp == "downstream": #implement upstream some other time if necessary
            if randomPAMdict[key] == "f":
                startPos = int(key) - 30 - NCounter #start position is 30 nucleotides upstream of PAM + the N's that were pruned in the beginning
                stopPos = int(key)-NCounter #stop is simply the start of PAM (taking N-pruning into account)
                spacer = g[startPos:stopPos]
                spacer_record = SeqRecord(spacer.seq)
                spacer_record.id = str(spacerCounter)
                spacer_record.description = "forward"
                #print (spacer)
                spacerList.append(spacer_record)
                #print (spacer)

            if randomPAMdict[key] == "r": #Slightly more complicated for reverse hits
                startPos = int(key) + NCounter + len(p)
                stopPos = int(key) + 30 + NCounter + len(p)
                spacer = g.seq[startPos:stopPos]
                #print (spacer)
                spacer = SeqRecord(spacer)
                spacerR = spacer.reverse_complement() #Important to reverse complement final sequence
                spacerR_record = SeqRecord(spacerR.seq)
                spacerR_record.id = str(spacerCounter)
                spacerR_record.description = "reverse"
                #print ("Flipping")
                #print (spacerR.seq)
                spacerList.append(spacerR_record)
                #print (spacerR_record)
                #print (spacerList)
        spacerCounter = spacerCounter + 1
        #print(str(spacerCounter))

    #print (spacerList)
    return spacerList

#    print (randomPAMdict)

def spacersToFasta(s):
    SeqIO.write(s, "randomizedSpacers.fasta", "fasta")


args = parser.parse_args() #Compiles the parser

genomes =  list(SeqIO.parse(args.genome, "fasta"))
genomeParsed = genomes[0]


if (args.mode == "PAM"):
#    print ("Sampling " + str(args.count) + " spacers from " + args.PAMpolarity + " PAM (" + PAMseq + ") sites from " + args.genome)
    PAMseq = Seq(args.PAMseq, generic_dna)
    spacerList = createPAMspacers(g = genomeParsed, c = args.count, p = PAMseq, pp = args.PAMpolarity)

if (args.mode == "random"):
    print ("Sampling " + str(args.count) + " spacers across genome from "  + args.genome)
    spacerList = []
    strands = ["f","r"]
    strandCount = {}
    spacerCounter = 0
    allLocations = []
    ntCounter = 0
    for nt in range(len(genomeParsed.seq)): #Get possible positions
        allLocations.append(ntCounter)
        ntCounter = ntCounter + 1
    #print(allLocations)

    randomLocations = random.sample(set(list(allLocations)), args.count) #random() requires a list, so we must break down the dictionary structure

    randomSpacersDic = {}

    for loc in randomLocations:
        strand = random.choice(["f","r"])
        randomSpacersDic[loc] = strand

        if randomSpacersDic[loc] == "f":
            startPos = int(loc)
            stopPos = startPos + 30
            spacer = genomeParsed[startPos:stopPos]
            spacer_record = SeqRecord(spacer.seq)
            spacer_record.id = str(spacerCounter)
            spacer_record.description = "forward"
            #print (spacer)
            spacerList.append(spacer_record)
            #print (spacer)

        if randomSpacersDic[loc] == "r":
            startPos = int(loc-30)
            stopPos = startPos + 30
            spacer = genomeParsed[startPos:stopPos]
            spacerR = spacer.reverse_complement() #Reverse complement final sequence
            spacerR_record = SeqRecord(spacerR.seq)
            spacerR_record.id = str(spacerCounter)
            spacerR_record.description = "reverse"
            #print (spacer)
            spacerList.append(spacerR_record)
            #print (spacer)

        spacerCounter = spacerCounter + 1
spacersToFasta(spacerList)

#spacersToFasta(spacers)
