
configfile: "globalConfig.yaml"

####################################
#check start contents for name matching, copy fastq files to data directory
rule initialize:
    input:
        "cassettes/cassette_{runName}/{runName}_R1.fastq.gz",
        "cassettes/cassette_{runName}/{runName}_R2.fastq.gz",
        "cassettes/cassette_{runName}/{runName}_barcodes.txt"
    output:
        "data/{runName}/done.initialize",
        "data/{runName}/{runName}_R1.fastq.gz",
        "data/{runName}/{runName}_R2.fastq.gz"
    
    script:
        "scripts/initialize.py"
        #print(f"{wildcards.runName}")
        #print(f"{output[0]}")
        
rule checkFastqFiles:
    input:
        fw= "data/{runName}/{runName}_R1.fastq.gz",
        rv= "data/{runName}/{runName}_R2.fastq.gz"
    output:
        "data/{runName}/fastqFileValidation.txt"
    
    conda:
        "envs/qcCheck.yaml"
    shell:
        '''
        #fastq_info {input.fw} {input.rv} > {output} 2>&1
        touch {output}
        echo "Skipping fastqcheck"
        if [ $? -ne 0 ]; then
            echo "The fastq file validation failed for run: {wildcards.runName}"
            exit 0
        fi
        '''

rule fastqc:
    input:
        fw= "data/{runName}/{runName}_R1.fastq.gz",
        rv= "data/{runName}/{runName}_R2.fastq.gz",
        check= "data/{runName}/fastqFileValidation.txt"
    output:
        "data/{runName}/done.fastQC"
    log:
        "logs/{runName}/QC/fastQC.log"
    conda:
        "envs/qcFast.yaml"
    threads:
        config["qcThreads"]
    shell:
        '''
        mkdir -p data/{wildcards.runName}/fastQC/
        fastqc \
            -t {threads} \
            -o data/{wildcards.runName}/fastQC \
            {input.fw} \
            {input.rv} \
            > {log} 2>&1
        
        touch {output}
        '''
rule multiQC:
    input:
        "data/{runName}/done.fastQC"
    output:
        "data/{runName}/done.multiQC"
    log:
        "logs/{runName}/QC/multiQC.log"
    conda:
        "envs/qcMulti.yaml"
    threads:
        config["qcThreads"]
    shell:
        '''
        cd data/{wildcards.runName}/fastQC
	print "TEST: before multiqc"

        multiqc ./ 2> ../../../{log}

        if [ $? -eq 0 ]; then
            touch ../../../{output}
        fi

	print "TEST: after multiqc"

        '''

#create barcode file with run data for sabre
rule barcodes:
    input:
        "data/{runName}/done.initialize",
        "cassettes/cassette_{runName}/{runName}_barcodes.txt"
    output:
        "data/{runName}/sabreBarcodes.txt"

    script:
        "scripts/barcodes.py"


#split barcode file into multiple files, set by config parallel
#this rule pulls a runConfig parameter through input. 
#and has a runConfig usage example in script
rule parallelBarcodes:
    input:
        "data/{runName}/sabreBarcodes.txt",
        "cassettes/cassette_{runName}/runConfig.yaml"
    output:
        "data/{runName}/sabreParallelBarcodes/done.parallelBarcodes"

    params:
        par= config["sabreParallel"]
    script:
        "scripts/barcodeSplit.py"

def gather_setup(wildcards):
    
    #collect wildcards from propogation
    rn= wildcards.runName
    
    #construct output requirements and propogate wildcards
    #    initialize
    #rstr= f"data/{rn}/done.initialize"
    #    check
    #rstr= f"data/{rn}/fastqFileValidation.txt"
    #    fastQC
    #rstr= f"data/{rn}/done.fastQC"
    #    multiQC
    #rstr= f"data/{rn}/done.multiQC"
    #    barcodes
    #rstr = f"data/{rn}/sabreBarcodes.txt"
    #    parallelBarcodes
    #rstr = f"data/{rn}/sabreParallelBarcodes/done.parallelBarcodes"
    #print(f"TEST rstr: {rstr}")
    #return rstr
    
    #    all endpoints
    #    ;;;switched done.fastQC for done.multiQC
    ends= ["done.fastQC", "sabreParallelBarcodes/done.parallelBarcodes"]
    rstr= f"data/{rn}/"+"{ends}"
    return expand(rstr, ends= ends)

    
checkpoint setup:
    input:
        gather_setup
    output:
        "data/{runName}/done.setup"
    shell:
        '''
        echo "marker create: {output}"
        touch {output}
        '''

################################


#produces <runName>_<varietyName>_R<1|2>.fq files
rule demultiplexReads:
    input:
        fw= "data/{runName}/{runName}_R1.fastq.gz",
        rv= "data/{runName}/{runName}_R2.fastq.gz",
        bc= "data/{runName}/sabreParallelBarcodes/parallelBarcodes_{n}",
        cp= "data/{runName}/done.setup"
    output:
        d = "data/{runName}/done.demultiplexReads_{n}",
        u = "data/{runName}/sabreParallelBarcodes/unknown_R1_{n}.barcodes",
        w = "data/{runName}/sabreParallelBarcodes/unknown_R2_{n}.barcodes"
    log:
        "logs/{runName}/sabre/demultiplexReads_{n}.log"
    wildcard_constraints:
        n="\d+"
    conda:
        "envs/demultiplexReads.yaml"
    shell:
        '''
        #navigate down to force sabre execution output directory
        [ ! -d "data/{wildcards.runName}/reads" ] && \
            mkdir "data/{wildcards.runName}/reads" && \
                echo "Sucesfully created reads directory"
        cd data/{wildcards.runName}/reads
       
        #push paths up to start designation from standard execution level 
        sabre pe -c \
            -f ../../../{input.fw} \
            -r ../../../{input.rv} \
            -b ../../../{input.bc} \
            -u ../../../{output.u} \
            -w ../../../{output.w} \
            1 > ../../../{log} \
            2 >> ../../../{log}

        touch ../../../{output.d}
        '''


def gather_demultiplexReads(wildcards):
    #print("TEST begin gather_demultiplexReads")
    
    #collect wildcards from propogation
    rn= wildcards.runName
    
    #collect wildcards from output
    gstr= f"data/{rn}/sabreParallelBarcodes/"+ "parallelBarcodes_{n}"
    wc= glob_wildcards(gstr)

    #remove duplicate wildcards
    n= list(set(wc.n))

    #construct output requirements and propogate wildcards
    dstr= f"data/{rn}/done.demultiplexReads" +"_{n}"
    #print(f"TEST dstr: {dstr}")
    #print("TEST end gather_demultiplexReads")

    return expand(dstr, n=n)

checkpoint point_demultiplexReads:
    input:
        "data/{runName}/done.setup",
        gather_demultiplexReads
    output:
        "data/{runName}/done.demultiplexReads_all"
    shell:
        '''
        #echo "TEST begin point_demultiplexReads"
        #collect unknown barcodes
        #filesR1=$(ls data/{wildcards.runName}/sabreParallelBarcodes/unknown_R1_[[:digit:]]*.barcodes)
        #filesR2=$(ls data/{wildcards.runName}/sabreParallelBarcodes/unknown_R2_[[:digit:]]*.barcodes)
        #clear target files
	#echo "" > data/{wildcards.runName}/unknown_R1.barcodes
	#echo "" > data/{wildcards.runName}/unknown_R2.barcodes
	#aggregate unknown barcodes
	#for f in $filesR1; do
        #    cat $f >> data/{wildcards.runName}/unknown_R1.barcodes
        #done
        #for f in $filesR2; do
        #    cat $f >> data/{wildcards.runName}/unknown_R2.barcodes
        #done

        echo "marker create: {output}"
        touch {output}
        '''

##################################

#removeAdaptors
#produces variety x 2 trimmed reads
rule removeAdaptors:
    input:
        r1= "data/{runName}/reads/{variety}_R1.fq",
        r2= "data/{runName}/reads/{variety}_R2.fq",
        cp= "data/{runName}/done.demultiplexReads_all"
    output:
        o= "data/{runName}/reads/{variety}_R1.fastq",
        p= "data/{runName}/reads/{variety}_R2.fastq"
    log:
        "logs/{runName}/cutadapt/{variety}.log"
    conda:
        "envs/removeAdaptors.yaml"
    params:
        adfw= config["adapfor"],
        adrv= config["adaprev"],
        readlen= config["readlen"]

    shell:
        '''
        cutadapt \
            -a {params.adfw} \
            -A {params.adrv} \
            -m {params.readlen} \
            -o {output.o} \
            -p {output.p} \
            {input.r1} \
            {input.r2} \
            > {log} 2>> {log}
        '''

def gather_removeAdaptors(wildcards):
    #print("TEST begin gather_removeAdaptors")
    
    #collect wildcards from propogation
    rn= wildcards.runName
    
    #collect wildcards from output
    gstr= f"data/{rn}/reads/"+ "{variety}_R{r}.fq"
    wc= glob_wildcards(gstr)

    #remove duplicate wildcards
    var= list(set(wc.variety))
    r= list(set(wc.r))
    #print(expand("TEST var:{var}", var= var))
    #print(expand("TEST r:{r}", r=r))
    
    #construct output requirements and propogate wildcards
    rstr= f"data/{rn}/reads/"+"{variety}_R{r}.fastq"

    #print("TEST gather_removeAdaptors return string:")
    #print(expand(rstr, variety= var, r= r))
    #print("TEST end gather_removeAdaptors")

    return expand(rstr, variety= var, r= r)

checkpoint point_removeAdaptors:
    input:
        "data/{runName}/done.demultiplexReads_all",
        gather_removeAdaptors
    output:
        "data/{runName}/done.removeAdaptors"
    shell:
        '''
        #echo "TEST begin point_removeAdaptors"
        echo "marker create: {output}"
        touch {output}
        '''

###########################################

rule cleanSamplesScript:
    input:
        "data/{runName}/done.removeAdaptors"
    output:
        "data/{runName}/varietyStats.txt"
    params:
        thresh= config["stdevThreshold"]
    script:
        "scripts/removeLowCount.py"

checkpoint cleanSamples:
    input:
        "data/{runName}/done.removeAdaptors",
        "data/{runName}/varietyStats.txt"
    output:
        "data/{runName}/done.cleanSamples"
    shell:
        '''
        echo "marker create: {output}"
        touch {output}
        '''

#############################################

rule alignReads:
    input:
        r1= "data/{runName}/reads/{variety}_R1.fastq",
        r2= "data/{runName}/reads/{variety}_R2.fastq",
        cp= "data/{runName}/done.cleanSamples"
    output:
        "data/{runName}/chr{chr}H/{variety}.sam"
    log:
        "logs/{runName}/bwa/{variety}_chr{chr}H.log"
    conda:
        "envs/align.yaml"
    threads: config["bwaThreads"]
    params:
        rga= config["refGenA"],
        rgb= config["refGenB"]
    shell:
        '''
        bwa mem \
            -t {threads} \
            cassettes/cassette_{wildcards.runName}/refgenome/{params.rga}{wildcards.chr}{params.rgb} \
            {input.r1} \
            {input.r2} \
            > {output} \
            2> {log}
        '''

rule markDuplicates: 
    input:
        "data/{runName}/chr{chr}H/{variety}.sam"
    output:
        "data/{runName}/chr{chr}H/{variety}.bam"
    log:
        "logs/{runName}/markDuplicates/{variety}_chr{chr}H.log"
    conda:
        "envs/samtools.yaml"
    threads: config["deDupThreads"]
    params:
    shell:
        '''
        samtools fixmate -@{threads} -O bam,level=0 -m {input} - | \
        samtools sort -@{threads} -O bam,level=0 -T tempSort_{wildcards.variety}_{wildcards.chr} - | \
        samtools markdup -@{threads} -O bam - {output} 2> {log}
        '''

rule indexBam: 
    input:
        "data/{runName}/chr{chr}H/{variety}.bam"
    output:
        "data/{runName}/chr{chr}H/{variety}.bam.csi"
    log:
        "logs/{runName}/indexBam/{variety}_chr{chr}H.log"
    conda:
        "envs/samtools.yaml"
    threads: config["deDupThreads"]
    params:
    shell:
        '''
        samtools index -@{threads} -c {input} 2> {log}
        '''

def gather_alignReads(wildcards):

    #collect wildcards from propogation
    rn= wildcards.runName

    #collect wildcards from output
    gstr= f"data/{rn}/reads/" +"{variety}_R{r}.fastq"
    wc= glob_wildcards(gstr)
    
    #collect wildcards from config and refgenome
    rga= config["refGenA"]
    rgb= config["refGenB"]
    gstr2= f"cassettes/cassette_{rn}/refgenome/{rga}" +"{chromosome}" +f"{rgb}"
    wc2= glob_wildcards(gstr2)
    
    #remove duplicate wildcards
    var= list(set(wc.variety))
    chr= list(set(wc2.chromosome))
    #print(expand("TEST var:{var}", var= var))
    #print(expand("TEST chr:{chr}", chr= chr))

    #construct input requirements
    #   bwa
    #rstr= f"data/{rn}/"+"chr{chr}H/{variety}.sam"
    #   markDuplicates
    #rstr= f"data/{rn}/"+"chr{chr}H/{variety}.bam"
    #   indexBam
    rstr= f"data/{rn}/"+"chr{chr}H/{variety}.bam.csi"
    
    return expand(rstr, chr= chr, variety= var)

checkpoint point_alignReads:
    input:
        "data/{runName}/done.cleanSamples",
        gather_alignReads
    output:
        "data/{runName}/done.alignReads"
    shell:
        '''
        echo "marker create: {output}"
        touch {output}
        '''

#############################################

#rule listBam: produce a bam files list for every chr 
rule listBam: 
    input:
        "data/{runName}/done.alignReads"
    output:
        "data/{runName}/{chrDir}/listBam.txt"
    log:
        "logs/{runName}/listBam/{chrDir}.log"
    script:
        "scripts/listBam.py"

rule callVariants:
    input:
        "data/{runName}/chr{chr}H/listBam.txt"
    output:
        "data/{runName}/chr{chr}H/variants.vcf"
    log:
        "logs/{runName}/platypus/chr{chr}H.log"
    params:
        minMapQual=config["minMapQual"],
        genIndels=config["genIndels"],
        minBaseQual=config["minBaseQual"],
        minReads=config["minReads"],
        rgb= config["refGenB"],
        rga= config["refGenA"]
    threads:
        config["platThreads"]

    conda:
        "envs/platypus.yaml"
    shell:
        '''
        platypus callVariants \
        			--bamFiles={input} \
        	    	--nCPU={threads} \
                    --minMapQual={params.minMapQual} \
        			--minBaseQual={params.minBaseQual} \
        	    	--minGoodQualBases=5 \
        			--badReadsThreshold=10 \
        	    	--rmsmqThreshold=20 \
        			--abThreshold=0.01 \
        			--maxReadLength=250  \
        			--hapScoreThreshold=20 \
        	    	--trimAdapter=0 \
        			--maxGOF=20 \
        			--maxReads=500000000 \
        	    	--minReads={params.minReads} \
        			--genIndels={params.genIndels} \
                    --minFlank=5 \
        	    	--sbThreshold=0.01 \
        			--scThreshold=0.95 \
        			--hapScoreThreshold=15 \
        	    	--filterDuplicates=0 \
        	    	--filterVarsByCoverage=0 \
        			--filteredReadsFrac=0.7 \
        			--minVarFreq=0.002 \
        	    	--mergeClusteredVariants=0 \
        			--filterReadsWithUnmappedMates=0 \
                    --refFile=cassettes/cassette_{wildcards.runName}/refgenome/{params.rga}{wildcards.chr}{params.rgb} \
                    --logFileName={log} \
        	    	--output={output} \
                    2> {log}.stder
        '''

def gather_callVariants(wildcards):

    #collect wildcards from propogation
    rn= wildcards.runName

    #collect wildcards from output
    gstr= f"data/{rn}/" +"{chrDir}/{variety}.bam"
    wc= glob_wildcards(gstr)
    
    #remove duplicate wildcards
    var= list(set(wc.variety))
    chr= list(set(wc.chrDir))

    #construct input requirements
    #   listBam
    #rstr= f"data/{rn}/"+"{chr}/listBam.txt"
    #   callVariants
    rstr=f"data/{rn}/"+"{chr}/variants.vcf"

    return expand(rstr, chr= chr)

checkpoint point_callVariants:
    input:
        "data/{runName}/done.alignReads",
        gather_callVariants
    output:
        "data/{runName}/done.callVariants"
    shell:
        '''
        echo "marker create: {output}"
        touch {output}
        '''

#############################################

rule compress:
    input:
        "data/{runName}/done.callVariants"
    output:
        "data/{runName}/chr{chr}H/variants.vcf.gz"
    log:
        "logs/{runName}/prepareMerge/chr{chr}H_bgzip.log"

    conda:
        "envs/bcftools.yaml"
    threads:
        config["prepThreads"]
    shell:
        '''
        #echo "TEST begin compress"
        input=$(echo {output} | sed s:.gz::)
        bgzip \
            -@ {threads} \
            -stdout \
            $input \
            > {output} \
            2> {log}
        '''

rule tbiIndex:
    input:
        "data/{runName}/chr{chr}H/variants.vcf.gz"
    output:
        "data/{runName}/chr{chr}H/variants.vcf.gz.tbi"
    log:
        "logs/{runName}/prepareMerge/chr{chr}H_indextbi.log"

    conda:
        "envs/bcftools.yaml"
    threads:
        config["prepThreads"]
    shell:
        '''
        #echo "TEST begin tbiIndex"
        bcftools index \
            --threads {threads} \
            -f -t \
            -o {output} \
            {input} \
            > {log} 2>&1 \
	    || echo "tbi Index Error caught" >> {log}
	#...fake output may interfere downstream
	touch {output}
        '''

rule csiIndex:
    input:
        "data/{runName}/chr{chr}H/variants.vcf.gz"
    output:
        "data/{runName}/chr{chr}H/variants.vcf.gz.csi"
    log:
        "logs/{runName}/prepareMerge/chr{chr}H_indexcsi.log"

    conda:
        "envs/bcftools.yaml"
    threads:
        config["prepThreads"]
    shell:
        '''
        #echo "TEST begin csiIndex"
        bcftools index \
            --threads {threads} \
            -f -c \
            -o {output} \
            {input} \
            > {log} 2>&1
        '''

rule convert:
    input:
       var= "data/{runName}/chr{chr}H/variants.vcf.gz",
       tbi= "data/{runName}/chr{chr}H/variants.vcf.gz.tbi",
       csi= "data/{runName}/chr{chr}H/variants.vcf.gz.csi"
    output:
        "data/{runName}/chr{chr}H/variants.bcf"
    log:
        "logs/{runName}/prepareMerge/chr{chr}H_convert.log"

    conda:
        "envs/bcftools.yaml"
    threads:
        config["prepThreads"]
    shell:
        '''
        #echo "TEST begin convert"
        bcftools convert \
            --threads {threads} \
            -O b \
            -o {output} \
            {input.var} \
            2> {log}
        '''

rule sort:
    input:
        "data/{runName}/chr{chr}H/variants.bcf"
    output:
        "data/{runName}/chr{chr}H/variants.bcf.sort"
    log:
        "logs/{runName}/prepareMerge/chr{chr}H_sort.log"

    conda:
        "envs/bcftools.yaml"
    shell:
        '''
        #echo "TEST begin sort"
        bcftools sort \
            -O b \
            -o {output} \
            {input} \
            2> {log}
        '''

rule tbiIndexBcf:
    input:
        "data/{runName}/chr{chr}H/variants.bcf.sort"
    output:
        "data/{runName}/chr{chr}H/variants.bcf.sort.tbi"
    log:
        "logs/{runName}/prepareMerge/chr{chr}H_indextbi_bcf.log"

    conda:
        "envs/bcftools.yaml"
    threads:
        config["prepThreads"]
    shell:
        '''
        #echo "TEST begin tbiIndex"
        bcftools index \
            --threads {threads} \
            -f -t \
            -o {output} \
            {input} \
            > {log} 2>&1 \
	    || echo "tbi Index Error caught" >> {log}
	#...fake output may interfere downstream
	touch {output}
        '''

rule csiIndexBcf:
    input:
        "data/{runName}/chr{chr}H/variants.bcf.sort"
    output:
        "data/{runName}/chr{chr}H/variants.bcf.sort.csi"
    log:
        "logs/{runName}/prepareMerge/chr{chr}H_indexcsi_bcf.log"

    conda:
        "envs/bcftools.yaml"
    threads:
        config["prepThreads"]
    shell:
        '''
        #echo "TEST begin csiIndex"
        bcftools index \
            --threads {threads} \
            -f -c \
            -o {output} \
            {input} \
            > {log} 2>&1
        '''

def gather_prepareMerge(wildcards):
    #collect wildcards from propogation
    rn= wildcards.runName

    #collect wildcards from output
    gstr= f"data/{rn}/" +"chr{chr}H/variants.vcf"
    wc= glob_wildcards(gstr)
    
    #remove duplicate wildcards
    chr= list(set(wc.chr))

    #construct input requirements
    #    compress
    #rstr= f"data/{rn}/"+"chr{chr}H/variants.vcf.gz"
    #    tbiIndex
    #rstr= f"data/{rn}/"+"chr{chr}H/variants.vcf.gz.tbi"
    #    csiIndex
    #rstr= f"data/{rn}/"+"chr{chr}H/variants.vcf.gz.csi"
    #    convert
    #rstr= f"data/{rn}/"+"chr{chr}H/variants.bcf"
    #    sort
    #rstr= f"data/{rn}/"+"chr{chr}H/variants.bcf.sort"
    #    tbiIndexBcf
    #rstr= f"data/{rn}/"+"chr{chr}H/variants.bcf.sort.tbi"
    #    csiIndexBcf
    #rstr= f"data/{rn}/"+"chr{chr}H/variants.bcf.sort.csi"

    #return expand(rstr, chr= chr)
    
    #    all requirements
    #        endpoints
    #    variants.bcf.sort
    #    variants.bcf.sort.tbi
    #    variants.bcf.sort.csi
    
    ends= [".bcf.sort", ".bcf.sort.tbi", ".bcf.sort.csi"]
    rstr= f"data/{rn}/"+"chr{chr}H/variants{ends}"
    return expand(rstr, chr= chr, ends= ends)

checkpoint prepareMerge:
    input:
        "data/{runName}/done.callVariants",
        gather_prepareMerge
    output:
        "data/{runName}/done.prepareMerge"
    shell:
        '''
        echo "marker create: {output}"
        touch {output}
        '''

##################################

rule  concatenate:
    input:
        "data/{runName}/done.prepareMerge"
    output:
        "data/{runName}/variantsMerged.vcf"
    log: 
        "logs/{runName}/merge/merge.log"
    threads:
        config["mergeThreads"]
    conda:
        "envs/bcftools.yaml"
    script:
        "scripts/merge.py"

checkpoint mergeVariants:
    input:
        "data/{runName}/done.prepareMerge",
        "data/{runName}/variantsMerged.vcf"
    output:
        "data/{runName}/done.mergeVariants"
    shell:
        '''
        echo "marker create: {output}"
        touch {output}
        '''

#############################################

rule filter:
    input:
        "data/{runName}/done.mergeVariants"
    output:
        "data/{runName}/variantsFiltered.vcf"
    log:
        "logs/{runName}/vcftools/vcftools.log"
    params:
        mm= config["max-missing"],
        maf= config["minorAlleleFreq"]
    conda:
        "envs/impute.yaml"
    shell:
        '''
        vcftools \
        	--vcf data/{wildcards.runName}/variantsMerged.vcf \
        	--remove-filtered-all \
        	--max-missing {params.mm} \
        	--maf {params.maf} \
        	--remove-indels \
        	--mac 1 \
        	--min-alleles 2 \
        	--max-alleles 2 \
        	--recode \
        	--out {output} \
            2> {log}    
        
            mv {output}.recode.vcf {output}
        
        '''

rule renameVariants:
    input:
        "data/{runName}/variantsFiltered.vcf"
    output:
        "data/{runName}/variantsFilteredRenamed.vcf"
    script:
        "scripts/renameSNP.py"

rule infoGT:
    input:
        "data/{runName}/variantsFilteredRenamed.vcf"
    output:
        "data/{runName}/variants_GTinfo.vcf"
    conda:
        "envs/impute.yaml"
    shell:
        '''
        vcftools \
            --vcf {input} \
            --extract-FORMAT-info GT \
            --out {output}
        
        mv {output}.GT.FORMAT {output}
        '''

rule summary:
    input:
        "data/{runName}/variants_GTinfo.vcf"
    output:
        "data/{runName}/variantsSummary.txt"
    script:
        "scripts/summarize.py"

rule impute:
    input:
        "data/{runName}/variantsFilteredRenamed.vcf"
    output:
        "data/{runName}/variantsImputed.vcf"
    log:
        "logs/{runName}/beagle/beagle.log"
    params:
        mem= config["beagleXmx"]
    conda:
        "envs/impute.yaml"
    shell:
        '''
        beagle -Xmx{params.mem}G \
            gt={input} \
            out={output} \
            > {log} 2>&1

        rm {output}.log
        gzip -d {output}.vcf.gz
        mv {output}.vcf {output}

        '''

rule copyResults:
    input:
        "data/{runName}/variantsImputed.vcf",
	"data/{runName}/variantsSummary.txt"

    output:
        "data/{runName}/done.copyResults"
    shell:
        '''
        echo "TEST Begin: copyResults"

        #collect large results files
        files="{input}"
        loc="data/{wildcards.runName}"
        files="${{files}} ${{loc}}/variantsMerged.vcf"

        #parallel transfer of large results files to results
        echo "${{files}}" | tr -d '\n' | xargs -d ' ' -I {{}} -P 5 -n 1 rsync -Pavh {{}} results/{wildcards.runName}/

        #move logs and extract output2.log output.log
        cp data/{wildcards.runName}/variantsSummary.txt logs/{wildcards.runName}/variantsSummary.txt
        cp data/{wildcards.runName}/varietyStats.txt logs/{wildcards.runName}/varietyStats.txt

        cat logs/output2.log | grep {wildcards.runName} > logs/{wildcards.runName}/output2.run.log
        echo $(cat logs/output.log | grep {wildcards.runName}) > logs/{wildcards.runName}/output.run.log

        #sync logs to results
        rsync -Pavh logs/{wildcards.runName}/ results/{wildcards.runName}/logs

        #transfer directories to results
        rsync -Pavh ${{loc}}/fastQC results/{wildcards.runName}

        touch {output}
        '''

def gather_finalize(wildcards):
    #collect wildcards from propogation
    rn= wildcards.runName
    
    #construct input requirements
    #    filter
    #rstr= f"data/{rn}/variantsFiltered.vcf"
    #    infoGT
    #rstr= f"data/{rn}/variantsFiltered_GTinfo.vcf"
    #    summary
    #rstr= f"data/{rn}/variantsSummary.txt"
    #    impute
    #rstr= f"data/{rn}/variantsImputed.vcf"
    #    renameVariants
    #rstr= f"data/{rn}/variantsImputedRenamed.vcf"
    #    copyResults
    #rstr= f"data/{rn}/done.copyResults"
    #return rstr
    
    #    all endpoints
        #variantsSummary.txt calls up to filter
        #variantsRenamed.vcf
        #done.copyResults
    ends= ["variantsSummary.txt", "done.copyResults"]
    rstr= f"data/{rn}/"+"{ends}"
    return expand(rstr, ends= ends)

checkpoint finalize:
    input:
        "data/{runName}/done.mergeVariants",
        gather_finalize
    output:
        "data/{runName}/done.finalize"
    shell:
        '''
        echo "marker create: {output}"
        touch {output}
        '''


#############################################

def gather_all(wildcards):
    #print("BEGIN gather_all")
    
    #Welcome to the advanced parallellization functions
    #USER set to true to use keys to control which cassettes are to be run
    useKeys = True
    #USER set to true to limit the number of runs to one
    singleRun = True


    rn = []
    if useKeys:
        wc= glob_wildcards("keys/key_{rn}")
        #remove duplicate wildcards
        rn= list(set(wc.rn))
    else:
        wc= glob_wildcards("cassettes/cassette_{rn}/{junkA}/{junkB}")
        #remove duplicate wildcards
        rn= list(set(wc.rn))

    #test for correct n of run names
    n= len(rn)
    if n == 0:
        raise Exception("There must be at least one key_<runName> file. See README.md")

    if singleRun and n > 1:
        if useKeys:
            raise Exception("There must be only one file key_<runName> See README.md")
        else:
            raise Exception("There must be only one dir cassette_<runName> See README.md")
    
    #print(f"TEST gather_all runName: {rn}")

    
    #create return string
    #   end at step
    #   setup
    #rstr = "data/{runName}/done.setup"
    #   point_demultiplexReads
    #rstr = "data/{runName}/done.demultiplexReads_all"
    #   point_removeAdaptors
    #rstr = "data/{runName}/done.rmAdaptors"
    #   cleanSamples
    #rstr = "data/{runName}/done.cleanSamples"
    #   alignReads
    #rstr = "data/{runName}/done.alignReads"
    #    callVariants
    #rstr = "data/{runName}/done.callVariants"
    #    prepareMerge
    #rstr = "data/{runName}/done.prepareMerge"
    #    merge
    #rstr = "data/{runName}/done.mergeVariants"
    #    finalize
    rstr = "data/{runName}/done.finalize"
    
    frstr= expand(rstr, runName = rn)
    print(f"NOTE gather_all returns:\n{frstr}")
    return frstr

rule all:
    input: gather_all
    shell: "echo BEGIN all"

