#!/usr/bin/env python3
import os
import sys 
import fileinput 

print("Begin: barcodes.py")

ip= snakemake.input[1]
op= snakemake.output[0]
rn= snakemake.wildcards.runName

#the barcode file ip is known to exist
with open(ip, 'r') as reader, open(op, 'w') as writer:
	line = reader.readline()
	while line:
		list = line.split()
		#Barcode runName_varietyName_R1.fq runName_varietyName_R2.fq
		#writer.write(list[0]+f" {rn}_"+list[1]+"_R1.fq"+f" {rn}_"+list[1]+"_R2.fq\n")
		#Easier to read
		#Barcode varietyName_R1.fq varietyName_R2.fq
		writer.write(list[0]+" "+list[1]+"_R1.fq"+" "+list[1]+"_R2.fq\n")
		line = reader.readline()
	