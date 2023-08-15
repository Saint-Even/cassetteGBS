#!/usr/bin/python

"""
Requires an input file from vcftools produced by:
    vcftools --vcf <filteredVariantsNotImputed.vcf> --extract-FORMAT-info GT
"""

import sys
from decimal import Decimal

print("Begin: summarize.py")

ip= snakemake.input[0]
op= snakemake.output[0]
#print(f"ip: {ip}")
#print(f"op: {op}")

#declare lists and dictionaries
samples={}
listSamples=[]
sites={}
listSites=[]

with open(ip, 'r') as f:
    #Get sample name info from the first line
    line = f.readline()
    listSamples = line.split()[2:]
    for y in listSamples:
        samples[y]=[]

    #Take the info for the loci (sites)
    line = f.readline()
    while line:
            linesplit=line.split()
            chr = linesplit[0]
            pos = linesplit[1]
            listGT = linesplit[2:]
            #Fill the sites dict
            listSites.append(chr+'_'+pos)
            sites[chr+'_'+pos] = listGT
            #Fill the samples dict
            a=0
            while a < len(listGT):
                    samples[listSamples[a]].append(listGT[a])
                    a+=1

            line = f.readline()

#Output the info
out=open(op,'w')

#FIRST section

#out.write("Individual info for Sites\n")
#out.write('Chr_Pos\t#Samples\t#homozygotes\t#heterozygotes\thetZ/homoZ\t#other\t#missing\n')

total_homo1= 0
total_homo2= 0
total_het1= 0
total_het2= 0
total_miss= 0
total_others= 0
total_sites= 0

for s in listSites:

        nb_homo1=sites[s].count('0/0')
        total_homo1+= nb_homo1
        nb_homo2=sites[s].count('1/1')
        total_homo2+= nb_homo2
        nb_het1=sites[s].count('1/0')
        total_het1+= nb_het1
        nb_het2=sites[s].count('0/1')
        total_het2+= nb_het2
        nb_miss=sites[s].count('./.')
        total_miss+= nb_miss
        nb_others=len(sites[s])-nb_homo1-nb_homo2-nb_het1-nb_het2-nb_miss
        total_others+= nb_others
        total_sites+= len(sites[s])

        if (nb_homo1+nb_homo2 == 0 ):
            ratio= Decimal('nan')
        else:
            ratio=(nb_het1+nb_het2+0.0)/(nb_homo1+nb_homo2+0.0)

#        Inididual info
#        out.write(s+'\t'+str(len(sites[s]))+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(ratio)+'\t'+str(nb_others)+'\t'+str(nb_miss)+'\n')


#convert ints to float for division
total_homo1+= 0.0
total_homo2+= 0.0
total_het1+= 0.0
total_het2+= 0.0
total_miss+= 0.0
total_others+= 0.0
total_sites+= 0.0


out.write('Accumulated date by site\n\n')

out.write('#Samples\tHomozygotes\tHeterozygotes\tOtherGenotypes\tMissing\n')
out.write(str(total_sites)+'\t'+str(total_homo1+total_homo2)+'\t'+str(total_het1+total_het2)+'\t'+str(total_others)+'\t'+str(total_miss)+'\n')

hetZoT= (total_het1+total_het2)/total_sites
homZoT= (total_homo1+total_homo2)/total_sites
hetZhomZ= (total_het1+total_het2)/(total_homo1+total_homo2)
OoT= total_others/total_sites
MoT= total_miss/total_sites

out.write('HetZ/Total\tHomZ/Total\tHetZ/HomZ\tOther/Total\tMissing/Total\n')
out.write(str(hetZoT)+'\t'+str(homZoT)+'\t'+str(hetZhomZ)+'\t'+str(OoT)+'\t'+str(MoT)+'\n\n')



#SECOND section
out.write("Samples Data\n")
out.write('Samples\t#Sites\t#homozygotes\t#heterozygotes\thetZ/homoZ\t#other\t#missing\n')

total_homo1= 0
total_homo2= 0
total_het1= 0
total_het2= 0
total_miss= 0
total_others= 0
total_samples= 0

for s in listSamples:
        nb_homo1=samples[s].count('0/0')
        total_homo1+= nb_homo1
        nb_homo2=samples[s].count('1/1')
        total_homo2+= nb_homo2
        nb_het1=samples[s].count('1/0')
        total_het1+= nb_het1
        nb_het2=samples[s].count('0/1')
        total_het2+= nb_het2
        nb_miss=samples[s].count('./.')
        total_miss+= nb_miss
        nb_others=len(samples[s])-nb_homo1-nb_homo2-nb_het1-nb_het2-nb_miss
        total_others+= nb_others
        total_samples+= len(samples[s])
        if (nb_homo1+nb_homo2 == 0 ):
            ratio= Decimal('nan')
        else:
            ratio=(nb_het1+nb_het2+0.0)/(nb_homo1+nb_homo2+0.0)

        out.write(s+'\t'+str(len(samples[s]))+'\t'+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(ratio)+'\t'+str(nb_others)+'\t'+str(nb_miss)+'\n')

#convert int to float
total_homo1+= 0.0
total_homo2+= 0.0
total_het1+= 0.0
total_het2+= 0.0
total_miss+= 0.0
total_others+= 0.0
total_samples+= 0.0

out.write('accumulated data by sample\n\n')

out.write('#Samples\t#homozygotes\t#heterozygotes\t#other\t#missing\n')
out.write(str(total_samples)+'\t'+str(total_homo1+total_homo2)+'\t'+str(total_het1+total_het2)+'\t'+str(total_others)+'\t'+str(total_miss)+'\n')

hetZoT= (total_het1+total_het2)/total_samples
homZoT= (total_homo1+total_homo2)/total_samples
hetZhomZ= (total_het1+total_het2)/(total_homo1+total_homo2)
OoT= total_others/total_samples
MoT= total_miss/total_samples

out.write('HetZ/Total\tHomZ/Total\tHetZ/HomZ\tOther/Total\tMissing/Total\n')
out.write(str(hetZoT)+'\t'+str(homZoT)+'\t'+str(hetZhomZ)+'\t'+str(OoT)+'\t'+str(MoT)+'\n\n')

out.close()
