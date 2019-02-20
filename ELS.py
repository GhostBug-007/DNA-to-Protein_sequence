#Contributed by GhostBug

data_file = open("9915_ref_Bos_indicus_1.0_chr3.fa",'r')
file_write = open("Output_proteinSequence.txt",'a')

data_file=data_file.read()
length =len(data_file)


'''
This function is for finding the next start point 
'''
def findActive(i):
	n = i
	for l in range(n,length):
		if(i > length-4):
			break
		if data_file[i] == "\n":
			i=i+1
			continue
		elif data_file[i+1] == "\n" :
			code = data_file[i]+data_file[i+2:i+4]
			if code == "ATG" :
				return i
			i=i+1
			continue
		elif data_file[i+2] == "\n" :
			code = data_file[i:i+2]+data_file[i+3]
			if code == "ATG" :
				return i
			i=i+1
			continue

		code = data_file[i:i+3]
		if code == "ATG" :
			return i

		i=i+1
	return length


'''
This function is to read the output file and find the minimum and the maximum length protein
'''	
def findMinMax():
	fileOpen = open("Output_proteinSequence.txt",'r')
	fileRead = fileOpen.read()
	minIndex , maxIndex = 0,0
	count =0
	lenFile = len(fileRead)
	minCount , maxCount = lenFile ,0
	for i in range(lenFile):
		count= count + 1
		if fileRead[i] == "_":
			if count > maxCount:
				maxCount = count 
				maxIndex = i+1 - count
			elif count < minCount and count > 3:			#count > 3 ; This condition has been used to ensure that we obtain a valid protein
				minCount = count 
				minIndex = i+1 - count 
			count = 0


	newFile = open("Minimum_Maximum_sequence.txt",'a')

	newFile.write("\nThe protein with maximum length is \n ")
	newFile.write(fileRead[maxIndex:maxIndex+maxCount])
	x = polarityCheck(fileRead[maxIndex:maxIndex+maxCount])
	newFile.write("\nThe percentage of Polar amino acids is "+str(x))

	newFile.write("\n\nThe protein with minimum length is \n")
	newFile.write(fileRead[minIndex:minIndex+minCount])
	x = polarityCheck(fileRead[minIndex:minIndex+minCount])
	newFile.write("\nThe percentage of Polar amino acids is "+str(x))

	fileOpen.close()
	return 0


'''
The function polarityCheck is to find the polar or non-polar nature of the protein with the maximum length
'''
def polarityCheck(polarList):
	polarityFile = open("polarity.txt",'r')
	polarLine = polarityFile.readline()
	polarDict = dict()
	countPolar = 0
	countNonP = 0
	j=20						# j = 20, because the total number of amino acids are 20
	while(j > 0):
		x = polarLine.split()
		polarLine = polarityFile.readline()
		j = j -1
		polarDict.update({x[0]:x[1]})

	for value in polarList:
		if value in polarDict and (polarDict[value] == "polar") :
			countPolar = countPolar+1
		elif value in polarDict and (polarDict[value] == "nonpolar"): 
			countNonP = countNonP+1
	ratio = (countPolar / (countNonP+countPolar))*100
	return ratio			

'''
This is the code for Retrieving the genetic code 
for various patterns and storing them in a dictonary
and using it as a key
'''
file = open("els_Genetic_code.txt")
line = file.readline()
line = line.split(",") 

for j in range(len(line)) :
	line[j] = line[j].strip()

geneticCode = dict()

for value in line:
	value = value.split("'")
	geneticCode.update({value[1]:value[3]})

'''
This part of the code is for matching of sequence
'''
code = 0
active = 0
i=0
for x in data_file :
	if active == 0:
		active = 1
		code = 0
		i=findActive(i)
	if code == "TAA" or code == "TAG" or code == "TGA" :
		active = 0
		continue

	if(i > length-4):
		break
	if data_file[i] == "\n":
		i=i+1
		continue
	elif data_file[i+1] == "\n" :
		code = data_file[i]+data_file[i+2:i+4]
		if code in geneticCode:
			file_write.write(geneticCode[code])
		i=i+4
		continue
	elif data_file[i+2] == "\n" :
		code = data_file[i:i+2]+data_file[i+3]
		if code in geneticCode:
			file_write.write(geneticCode[code])
		i=i+4
		continue

	code = data_file[i:i+3]
	if code in geneticCode:
		file_write.write(geneticCode[code])
	i=i+3


print("File Read Successfully \nYou will find the output protein sequence\nin the file named \"Output_proteinSequence.txt\"")
findMinMax()
print("\nYou will also find another file named\n \"Minimum_Maximum.txt\"\nwhich will contain the maximum\nand minimum length protein in the protein sequence")
print()

# *******************