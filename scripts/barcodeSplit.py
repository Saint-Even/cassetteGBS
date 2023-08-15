#!/usr/bin/env python3

import os
import sys 
import fileinput 
import math

print("Begin: barcodeSplit.py")

ip= snakemake.input[0]
check= snakemake.output[0]
op= os.path.dirname(check)
#get parallel Num
par= snakemake.params.par

####example of specific config parameters####
##open specific config
#import yaml
#configFile = snakemake.input[1]
#with open(configFile, "r") as f:
#    config = yaml.safe_load(f)
#str= config["specific"]
##print(str)

#count lines in input
with  open(ip,'r') as counter:
    for i, _ in enumerate(counter):
        pass
lines = i+1

#alert unknown behaivour
f= math.floor(lines/par)
if f==0:
    print(f"ALERT: parallel is set to {par} \n the number of lines in {ip} is {lines}")
    print(f"A single file of barcode info will be sent to Sabre and {par-1} Empty files")
    print("The behaviour of empty files has not been tested at this time.")

if not os.path.isdir(op):
    os.mkdir(f"{op}")

copiedLines =0
#copy ip lines into files, leaving the last file to get all remaining lines
with open(ip, "r") as reader:
    for i in range(par-1):
        with open(f"{op}/parallelBarcodes_{i}", "w") as writer:
            for l in range(f):
                writer.write(reader.readline())
                copiedLines += 1
    #finish all remaining lines
    with open(f"{op}/parallelBarcodes_{par-1}", "w") as writer:
        while copiedLines < lines:
            writer.write(reader.readline())
            copiedLines += 1

with open(check, 'w') as writer:
    writer.write("done")