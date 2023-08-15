import os, sys, fileinput

print("Begin: initialize.py")


ip0= snakemake.input[0]
ip1= snakemake.input[1]
op= snakemake.output[0]
rn= snakemake.wildcards.runName

#print(f"TEST ip0:{ip0} ip1:{ip1} op:{op} rn:{rn}")


os.system(f"cp {ip0} data/{rn}/{rn}_R1.fastq.gz")
os.system(f"cp {ip1} data/{rn}/{rn}_R2.fastq.gz")
with open(op, "w") as writer:
    writer.write("initialized")
       


