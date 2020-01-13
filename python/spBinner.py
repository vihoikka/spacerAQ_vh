'''
First created Aug 29, 2017
Takes a list of spacer sequences and maps them onto target genome.
Extracts start position of each spacer on the target genome.
Can be used to plot protospacer density plots in e.g R.

@author: ville hoikkala 2020
'''
import csv
import sys

def getSpacerPositions(file,direction):
	parsed = list(csv.reader(open(file, 'r'), delimiter='\t')) #opens file as tab-delimited
	posDic = {} #Create dictionary to store IDs as keys and positions as values
	lineCount = 0 #for skipping the header lines of the SAM file
	for line in parsed:
		if lineCount > 2: #Skip the headers
			if direction == "F" and int(line[1]) == 0: #0 for F hits, 16 for reverse hits
				posDic[line[0]] = int(line[3]) #"Position" in each line on column 3.
				#print(line[3])
			if direction == "R" and int(line[1]) == 16: #0 for F hits, 16 for reverse hits
				posDic[line[0]] = int(line[3]) #"Position" in each line on column 3
			#print line[3]
		lineCount = lineCount + 1
	return posDic

def fastaToDictionary(filename): #This takes a fasta file and creates a Python dictionary from it
	spacerD = {} #create empty dictionary
	spacers = open(filename) #open file
	spacerLines = spacers.readlines() #turn each line in file into a list item
	lineForLoopCounter = 0 #for iterating loop
	for i in spacerLines: #go through all lines
		if i[0] == ">": #if line starts with >, it's a header and dictionary key
			spacerD[i] = spacerLines[lineForLoopCounter+1] #set i as dictionary key and the following line as the value
		lineForLoopCounter += 1 #for iterating
	for i, value in spacerD.items(): #Remove \n from each value
		spacerD[i] = value.rstrip('\n')
	spacerDcorr = { x.replace('\n', ''): spacerD[x] for x in list(spacerD.keys()) } #same thing for all keys (different method here)

	return spacerDcorr

def openTargetFile(filename):
	genomeFile = open(filename) #open file
	genomeLines = genomeFile.readlines()
	genome = genomeLines[1]
	return genome

def removerOverlaps(positions,limit):
	#print("Removing overlapping spacers")
	uniquePositions = {}
	#Progress counters
	cycleCountDuplicate = 0
	tenPercentDuplicate = len(positions)*0.01
	progressCountDuplicate = 0

	simiCounter = 0 #For counting overlapping protospacers
	for key, value in positions.items():
		overlap = False
		for key2, value2 in uniquePositions.items():
			distance = abs(value - value2)
			#print (distance)
			if distance < limit:
				overlap = True
				simiCounter = simiCounter + 1
				break
		if overlap == False:
			uniquePositions[key] = value
		cycleCountDuplicate = cycleCountDuplicate + 1
		if cycleCountDuplicate > tenPercentDuplicate:
			cycleCountDuplicate = 0
			progressCountDuplicate = progressCountDuplicate + 1
			print ("   " + str(progressCountDuplicate) + " % of protospacers checked. Number of spacers:" + str(len(positions)) + ", of these skipped due to overlap: " + str(simiCounter) + "\r", end="")

	return uniquePositions

''' PROGRAM STARTS HERE '''

inputFile = sys.argv[1] #1st argument: input file
windowSize = int(sys.argv[2]) #size of bins in target genome
overlapLimit = int(sys.argv[3]) #3rd argument: amount of overlap. If starting coordinates of two spacers are closer than this, discard the other one (see function removerOverlaps). 0 for disabling.
perOrAbsolute = sys.argv[4] #4th argument: whether bin numbers are reported as absolute values or percentages (absolute or percent)
filename = sys.argv[5] #output filename
trgt = sys.argv[6] #genome to be mapped on (fasta file in the same folder)
numberOfWindows = int(sys.argv[7]) #if user determines window size manually, set this to 0 (number of bins is then automatically calculated. If user determines bin size, then this = 1 (number of bins is then automatically calculated - probably results in an incomplete last bin)
binFraction = float(sys.argv[8])
print(inputFile)

'''Obtaining positions on the target genome'''
spacersPosF = getSpacerPositions(inputFile,"F") #The values in this copy of the dictionary will contain spacer positions
spacersPosR = getSpacerPositions(inputFile,"R") #The values in this copy of the dictionary will contain spacer positions

'''If user wants to remove overlaps (e.g. duplicate protospacers), go through all positions and remove overlapping protospacers'''
if overlapLimit > 0:
	spacersPosF = removerOverlaps(spacersPosF,overlapLimit)
	spacersPosR = removerOverlaps(spacersPosR,overlapLimit)
print("---------------------")
print("Spacerbinner")


directions = [spacersPosF, spacersPosR]
spacersPosF[str] = "F"
spacersPosR[str] = "R"
mapped = {}

targetSequence = openTargetFile(trgt) #Open target genome from fasta file and

'''Binner. Uses a sliding window approach to bin protospacers in a given window size'''
if numberOfWindows == 0 and binFraction == 0: #If user wants to determine window size
	numberOfWindows = len(targetSequence)/windowSize
	print("Window size determined manually (" + str(windowSize + "). Number of windows: " + str(numberOfWindows)))
if numberOfWindows != 0 and binFraction == 0: #if user wants to determine number of windows..
	windowSize = round(len(targetSequence)/numberOfWindows) #.. then calculate appropriate window size
	print("Number of windows (" + str(numberOfWindows) + ") determined manually. Window size: " + str(windowSize))

if numberOfWindows == 0 and binFraction != 0:
	windowSizeRound = round(len(targetSequence)*binFraction)
	windowSize = len(targetSequence)*binFraction
	numberOfWindows = 1/binFraction
	excess = round(windowSize-windowSizeRound,3)*numberOfWindows
	lastBinSize = windowSize + excess*numberOfWindows #add the jakojaannos to the last bin
	print("Window size is " + str(round(windowSize)) + " with an excess of " + str(excess)  +" bp, which has been allocated to the last bin")

print("Binning genome " + trgt +  " (length " + str(len(targetSequence)) + " bp" +
	"). R hits: " + str(len(spacersPosR)) +
	", F hits: " + str(len(spacersPosF)) +
	". Bin size: "+ str(windowSize) +
	"; number of bins: " + str(numberOfWindows))
loopCount = 0
totalPercent = 0
binSizes = []
locations = []
for i in directions:
	strand = i[str]
	windows = {}
	lastWindow = max(range(int(numberOfWindows)))
	#print("Last window is " + str(lastWindow))
	for j in range(int(numberOfWindows)):
		#print(j)
		searchMax = (j+1)*windowSizeRound #Define search boundaries for the window in question
		if (j != 0 and j != lastWindow):
			searchMin = searchMax-windowSizeRound-1
			binSizes.append(searchMax-searchMin)
			locations.append(searchMin)
		if (j == lastWindow): #If we are at the last bin, add the ylijaama of ther bins to this one
			searchMax = (j+1)*windowSizeRound+excess
			#print("Searching for last. Size " + str(searchMax))
			searchMin = searchMax-windowSizeRound-excess-1
			#print ("Low limit " + str(searchMin))
			binSizes.append(searchMax-searchMin)
			locations.append(searchMin)
		if (j == 0):
			searchMin = searchMax-windowSizeRound
			binSizes.append(searchMax-searchMin)
			locations.append(searchMin)
		#print(searchMax)
		#print("min: " + str(searchMin))
		#print searchMax
		windows[j] = 0
		for key, value in i.items():
			if isinstance(value, int):
				if (value > searchMin and value <= searchMax):
					#print("FOUND: " + str(value))
					windows[j] = windows[j] + 1
					#print(windows[j])
		#windows[i] = None
		mapped[loopCount] = windows
	loopCount = loopCount + 1

	mostProtos = max(windows, key=windows.get)
	print("	Bin with most absolute spacers in " + strand + ": " + str(mostProtos) + " with " + str(windows[mostProtos]) + " hits") #Get window with most protospacers

	#print windows
if perOrAbsolute == "percent": #If user wants to report percentage
	total = 0
	print("Calculating percentages...")
	#Sum up all hits
	for totkey, totvalue in mapped.items():
		print ("Checking total from " + str(totkey))
		for totkey2, totvalue2 in totvalue.items():
			#print (totvalue2)
			total = total + totvalue2
		print ("Strand total " + str(total))
	print("	Total hits " + str(total))
	for strand, bins in mapped.items(): #After obtaining the total, individually calculate percentages for each window on both strands
		#print("Calculating percentage " + str(strand))
		for individualBin, count in bins.items(): #Go through each bin in the strand
			percentage = 0
			if total > 0:
				percentage = count/total*100
			bins[individualBin] = percentage
			totalPercent = totalPercent + percentage
			#print (str(bins[individualBin]))
print ("	Total percentage (for debugging purposes): " + str(totalPercent))

loopCounter = 0
binCounter = 0
for key, value in mapped.items(): #Save both strand files separately
	#print("Loopcounter " + str(loopCounter))
	if loopCounter == 0: outputname = filename + "_mappedF.csv"
	else: outputname = filename + "_mappedR.csv"
	with open(outputname, 'w') as csv_file:
		writer = csv.writer(csv_file)
		writer.writerow(["Pos","target","binSize","location"])
		for key2, value2 in value.items():
			writer.writerow([key2, value2, str(binSizes[binCounter]), str(locations[binCounter])])
			binCounter = binCounter + 1
	loopCounter = loopCounter + 1

print("Binning finished. Files saved as " + filename + "_mappedF.csv and " + filename + "_mappedR.csv")
print("---------")
