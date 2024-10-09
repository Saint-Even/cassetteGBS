#!/usr/bin/env python3

# %%
import pandas
import os
import glob

#This script takes a vcf from cassetteGBS and renames the id column using the columns filled by platypus variant caller. The specifics of how this is done need to be adapted to each reference genome, due to how platypus handles each in an idiosyncratic way. I have solved renaming for synergySRG and morex_V3

#We do not switch genomes frequently enough where it would be worthwhile to automate the many possibilities. Instead the script has a pipeline embedded mode and a standalone mode. If you get an error from the renaming step, copy this script to a separate dir, along with a copy of the input vcf as reported by the error. You may then comment in  the standalone mode lines, comment out the embedded lines and debug the namingcolumn manipulation.

# %%
# standalone usage: renameSNP.py <fileName.vcf> <outfileName.vcf>
#fileName = sys.argv[1]
#outName = sys.argv[2]
#numChromosomes = 7  # disable numCromosomes below

# snakemake embedded usage
fileName = snakemake.input[0]
outName = snakemake.output[0]

# use outName to get action directory, embedded use only
actDir = os.path.dirname(outName)
# get list of chromosome Dirs
getStr = os.path.join(actDir, "chr*H")
fileNames = glob.glob(getStr)
# get count of chromosomes
numChromosomes = len(fileNames)

# %%
# take header off the top
maxHeader = 150
header = list()
i = 0
with open(fileName, 'r') as reader:
    while i < maxHeader:
        line = reader.readline()
        if line.startswith('#'):
            header.append(line)
        else:
            i += 1

# %%
# keep body without header
f = open(fileName, 'r')
vcf = pandas.read_table(f, skiprows=len(header),
                        delim_whitespace=True, skip_blank_lines=True)

# %%
# identify the chromsome indicator character position
# manually check col 0 for morex
head = vcf.head()
tail = vcf.tail()

first = head.iloc[0, 0]
last = tail.iloc[4, 0]
indChromosome = 0
for i in range(len(first)):
    if first[i] == "0" and last[i] == str(numChromosomes-1):
        indChromosome = i

# %%
# correct platypus problematic SNP name col to only chr num

# VCF columns format
# 0:Chr, 1:postion, 2:id

# use for synergySRG
#vcf.iloc[:, 0] = vcf.iloc[:, 0].str[indChromosome]
#vcf.iloc[:, 0] = vcf.iloc[:, 0].astype('int')
#vcf.iloc[:, 0] = vcf.iloc[:, 0]+1

# use for morex
# str slice must be confirmed by manual check
# subtraction correcton must confirmed by manual check
vcf.iloc[:, 0] = vcf.iloc[:, 0].str[4:8]
vcf.iloc[:, 0] = vcf.iloc[:, 0].astype('int')
vcf.iloc[:, 0] = vcf.iloc[:, 0]-95

# confirm chr col is correct
# c = vcf.iloc[:, 0]
# print(c)

# %%
#make full SNP name
vcf.iloc[:,2] = "Chr" + vcf.iloc[:,0].astype(str) + "H_" + vcf.iloc[:,1].astype(str)

# %%
#write header and dataframe to new file
with open(outName, 'w') as writer:
    for line in header:
        writer.write(line)
    writer.write(vcf.to_csv(sep="\t", index=False, header=False))
