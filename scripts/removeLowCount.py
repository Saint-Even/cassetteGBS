
import sys, os, fileinput, math, glob
import statistics as stats

print("Begin: removeLowCount.py")

op= snakemake.output[0]
thresh= snakemake.params.thresh

#declare lists
fileNames= []
countReads= []

fileNamesDel= []
countReadsDel= []

#use op to glob fastq files
#data/{rn}/reads/" +"{variety}_R{r}.fastq"
actDir = os.path.dirname(op)
getStr= os.path.join(actDir, "reads/*.fastq*")
fileNames= glob.glob(getStr)




# for each fastq file
for i in range(len(fileNames)):
	#get file name
	fileName= fileNames[i]
	#add count of reads to list
	lines= 0
	with open(fileName, 'r') as counter:
		#count lines
		for index, _ in enumerate(counter):
			lines= index
	if lines != 0:
		lines += 1
	reads = lines/4
	countReads.append(reads)
	#print(f"TEST fileName:{fileName} lines:{lines} reads:{reads}")

#get non zero read counts
values= []
for i in range(len(countReads)):
	value= countReads[i]
	#print(f"TEST value: {value}")
	if value != 0:
		values.append(value)
		#print(f"TEST append value: {value}")


#Compute threshold
avgReads= stats.mean(values)
stdev= stats.stdev(values, avgReads)
cut= math.ceil(avgReads - thresh*stdev)
if cut < 1:
	cut = 1
#print(f"TEST avgReads:{avgReads} stdev:{stdev} cut:{cut}")

#identify samples files below the threshold
for i in range(len(countReads)):
	if countReads[i] < cut:
		fileNamesDel.append(fileNames[i])
		countReadsDel.append(countReads[i])

#mark samples files that do not exceed the threshold
for i in range(len(fileNamesDel)):
	fileName= fileNamesDel[i]
	#print(f"TEST fileName: {fileName}")
	#remove from passing list
	j= fileNames.index(fileName)
	fileNames.pop(j)
	countReads.pop(j)

	#asses extension for previous masking
	base= os.path.basename(fileName)
	#add exclusion mask   
	if not base.endswith(".excluded"):
		fileNameMask= f"{fileName}.excluded"
		os.rename(fileName, fileNameMask)

#present stats report
with open(op, 'w') as writer:
	writer.write("#Read-Count	Sample-File\n")
	
	writer.write("########Passing#########\n")
	for i in range(len(fileNames)):
		fileName= fileNames[i]
		count= countReads[i]
		writer.write(f"{count}	{fileName}\n")
	
	writer.write("########Excluded########\n")
	for i in range(len(fileNamesDel)):
		fileName= fileNamesDel[i]
		count= countReadsDel[i]
		writer.write(f"{count}	{fileName}\n")

#print("TEST end")
	