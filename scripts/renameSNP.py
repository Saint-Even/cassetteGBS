#!/usr/bin/env python3

# %%
import pandas
import sys
import re
import os
import glob

# %%
#standalone usage: renameSNP.py <fileName.vcf> <outfileName.vcf>
#fileName = sys.argv[1]
#outName = sys.argv[2]

#snakemake embedded usage
fileName= snakemake.input[0]
outName= snakemake.output[0]

#use outName to get action directory
actDir = os.path.dirname(outName)
#get list of chromosome Dirs
getStr= os.path.join(actDir, "chr*H")
fileNames= glob.glob(getStr)
#get count of chromosomes
numChromosomes = len(fileNames)

# %%
#take header off the top
maxHeader = 50
header = list()
i=0
with open(fileName, 'r') as reader:
    while i < maxHeader:
        line = reader.readline()
        if line.startswith('#'):
            header.append(line)
        else:
            i += 1


# %%
#keep body without header
f = open(fileName, 'r')
vcf = pandas.read_table(f, skiprows= len(header),
    delim_whitespace=True, skip_blank_lines=True)

# %%
#identify the chromsome indicator character position
head = vcf.head()
tail = vcf.tail()

first = head.iloc[0,0]
last = tail.iloc[4,0]

indChromosome = 0
for i in range(len(first)):
    if first[i] == "0" and last[i] == str(numChromosomes-1):
        indChromosome = i

# %%
#convert SNP name col to only chr num
vcf.iloc[:,0] = vcf.iloc[:,0].str[indChromosome]
vcf.iloc[:,0] = vcf.iloc[:,0].astype('int')
vcf.iloc[:,0] = vcf.iloc[:,0]+1

# %%
#make full SNP name
vcf.iloc[:,2] = "Chr" + vcf.iloc[:,0].astype(str) + "H_" + vcf.iloc[:,1].astype(str)

# %%
#write header and dataframe to new file
with open(outName, 'w') as writer:
    for line in header:
        writer.write(line)
    writer.write(vcf.to_csv(sep="\t", index=False, header=False))
