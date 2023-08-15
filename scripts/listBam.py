
import sys, os, fileinput, glob

print("Begin: listBam.py")

op= snakemake.output[0]

#use op to get action directory
actDir = os.path.dirname(op)
getStr= os.path.join(actDir, "*.bam")

#declare list
fileNames= glob.glob(getStr)

#generate list of bam files
#open bamList file
with open(op, 'w') as writer:
    #get full path of each file in fileNames
    for i in range(len(fileNames)):
        #get file name
        fileName= fileNames[i]
        apath= os.path.abspath(fileName)
        writer.write(f"{apath}\n")

