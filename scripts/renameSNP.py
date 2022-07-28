#!/usr/bin/env python3

# %%
import pandas
import sys
import re

# %%
#usage: convertImputedForTassel <fileName.vcf> <outfileName.vcf>
fileName = sys.argv[1]
outName = sys.argv[2]

maxHeader = 50
numChromosomes = 7

# %%
#take header off the top 
header = list()

with open(fileName, 'r') as reader:
    for i in range(maxHeader):
        line = reader.readline()
        if line.startswith('#'):
            header.append(line)


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
vcf.iloc[:,0] = vcf.iloc[:,0].astype(str)

# %%
#make full SNP name 
vcf.iloc[:,0] = "Chr" + vcf.iloc[:,0] + "H_" + vcf.iloc[:,1].astype(str)

# %%
#write header and dataframe to new file
with open(outName, 'a') as appender:
    for line in header:
        appender.write(line)
    appender.write(vcf.to_csv(sep="\t", index=False, header=False))
    
