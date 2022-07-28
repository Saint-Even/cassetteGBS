#!/usr/bin/python

"""
Requires input file from vcftools produced by:
	vcftools --vcf <filteredVariantsNotImputed.vcf> --extract-FORMAT-info GT

The output file wiil be: out.GT.FORMAT

Then, submit the out.GT.FORMAT to Summary4VCF.py:

Summary4VCF.py out.GT.FORMAT

"""

import sys
from decimal import Decimal


try:
	outgt=sys.argv[1]
	f = open(outgt,'r')
except:
	print(__doc__)
	sys.exit(1)

#Take the info for the samples
samples={}
listSamples=[]
line = f.readline()
listSamples = line.split()[2:]
for y in listSamples:
	samples[y]=[]

line = f.readline()

#Take the info for the loci (sites)
sites={}
listSites=[]
while line:
	chr = line.split()[0]
	pos = line.split()[1]
	listSites.append(chr+'_'+pos)
	listGT = line.split()[2:]
	#Fill the sites dict
	sites[chr+'_'+pos] = listGT
	
	#Fill the samples dict
	a=0
	while a < len(listGT):
		samples[listSamples[a]].append(listGT[a])
		a+=1
	
	line = f.readline()

#Output the info in 2 files
#one for the sites
out1=open('Summary_By_Sites_python.txt','w')
out1.write('Chr_Pos\t#Samples\t#homozygotes\t#heterozygotes\thetZ/homoZ\t#other\t#missing\n')

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


#	out1.write(s+'\t'+str(len(sites[s]))+'\t'+str(nb_homo1+nb_homo2+nb_het1+nb_het2+nb_miss)+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(nb_miss)+'\n')
	out1.write(s+'\t'+str(len(sites[s]))+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(ratio)+'\t'+str(nb_others)+'\t'+str(nb_miss)+'\n')


#convert ints to float for division
total_homo1+= 0.0
total_homo2+= 0.0
total_het1+= 0.0
total_het2+= 0.0
total_miss+= 0.0
total_others+= 0.0
total_sites+= 0.0




out3=open('Summary_Ratios_SitesAndPosition_python.txt','w')

out3.write('accumulated date by site\n\n')

out3.write('#Samples\t#homozygotes\t#heterozygotes\t#otherGenotypes\t#missing\n')
out3.write(str(total_sites)+'\t'+str(total_homo1+total_homo2)+'\t'+str(total_het1+total_het2)+'\t'+str(total_others)+'\t'+str(total_miss)+'\n')

hetZoT= (total_het1+total_het2)/total_sites
homZoT= (total_homo1+total_homo2)/total_sites
hetZhomZ= (total_het1+total_het2)/(total_homo1+total_homo2)
OoT= total_others/total_sites
MoT= total_miss/total_sites

out3.write('HetZ/Total\tHomZ/Total\tHetZ/HomZ\tOther/Total\tMissing/Total\n')
out3.write(str(hetZoT)+'\t'+str(homZoT)+'\t'+str(hetZhomZ)+'\t'+str(OoT)+'\t'+str(MoT)+'\n\n')

#the other for the samples
out2=open('Summary_By_Samples_python.txt','w')
out2.write('Samples\t#Sites\t#homozygotes\t#heterozygotes\thetZ/homoZ\t#other\t#missing\n')

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

#	out2.write(s+'\t'+str(len(samples[s]))+'\t'+str(nb_homo1+nb_homo2+nb_het1+nb_het2+nb_miss)+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(nb_miss)+'\n')
	out2.write(s+'\t'+str(len(samples[s]))+'\t'+'\t'+str(nb_homo1+nb_homo2)+'\t'+str(nb_het1+nb_het2)+'\t'+str(ratio)+'\t'+str(nb_others)+'\t'+str(nb_miss)+'\n')

#convert int to float
total_homo1+= 0.0
total_homo2+= 0.0
total_het1+= 0.0
total_het2+= 0.0
total_miss+= 0.0
total_others+= 0.0
total_samples+= 0.0

out3.write('accumulated date by sample\n\n')

out3.write('#Samples\t#homozygotes\t#heterozygotes\t#other\t#missing\n')
out3.write(str(total_samples)+'\t'+str(total_homo1+total_homo2)+'\t'+str(total_het1+total_het2)+'\t'+str(total_others)+'\t'+str(total_miss)+'\n')

hetZoT= (total_het1+total_het2)/total_samples
homZoT= (total_homo1+total_homo2)/total_samples
hetZhomZ= (total_het1+total_het2)/(total_homo1+total_homo2)
OoT= total_others/total_samples
MoT= total_miss/total_samples

out3.write('HetZ/Total\tHomZ/Total\tHetZ/HomZ\tOther/Total\tMissing/Total\n')
out3.write(str(hetZoT)+'\t'+str(homZoT)+'\t'+str(hetZhomZ)+'\t'+str(OoT)+'\t'+str(MoT)+'\n\n')

f.close()
out1.close()
out2.close()
out3.close()