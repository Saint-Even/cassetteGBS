import string
import sys, os, fileinput, glob

print("Begin: merge.py")

op= snakemake.output[0]
threads= snakemake.threads
log= snakemake.log

#use op to get action directory
actDir = os.path.dirname(op)
getStr= os.path.join(actDir, "chr*H/variants.bcf.sort")

#get list of files to merge
fileNames= glob.glob(getStr)
fileNames.sort()

fileStr= ""
for i in range(0, len(fileNames)):
    fileStr += f"{fileNames[i]} "

#call shell command
cmdStr= f"bcftools concat --threads {threads} -o {op} -O v {fileStr} 2> {log}"
#print (f"TEST cmdStr:\n{cmdStr}")
os.system(cmdStr)
