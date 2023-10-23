
    .------------------------.
    |\\////////       90 min |
    | \/  __  ______  __     |
    |    /  \|\.....|/  \    |
    |    \__/|/_____|\__/    |
    | A    Cassette GBS      |
    |    ________________    |
    |___/_._o________o_._\___|




# Cassette GBS

The design of this pipeline was inspired by watching my daughter drop a tape into a cassette player and press play. Only two steps to hear music, a truly user friendly interface! 

### Features
- GBS processing at scale
- Run setup designed for greatest possible simplicity (1 2 3 Music!)
- Massive process parallelization to maximize hardware usage and speed
- Complete reproducibility
- Detailed logging
- Multiple run modes for advanced usage on high performance hardware
##### Future Directions
- Configuration for individual runs
- Automatically run when input data is added

### Key concepts
- Cassette -> Key -> Run
- Cassette: a directory serving as a self-contained unit containing ALL initial inputs
- Key: a single filename which directs the pipeline to select a specific cassette
- Run: a script which invokes the full command that begins the pipeline
- The major calls and processing flow are built from the fastGBS_V2 snp calling pipeline
	- That documentation and citation is available [here](https://bitbucket.org/jerlar73/fast-gbs_v2/wiki/home)  
	- This pipeline is an elabouration on theirs in the snakemake dialect, to reduce system demands and accellerate processing by highly parallell processing.

### Cassette
Within the cassettes directory There may be multiple, uniquely named cassette directories

### Cassette
Within the cassettes directory  
There may be multiple, uniquely named cassette directories
These directories must be named  
`cassette_<runName>`

Each must contain one refgenome directory and the following three name-matched files.  

- `<runName>_barcodes.txt`
- `<runName>_R1.fastq.gz`
- `<runName>_R2.fastq.gz`
- `refgenome`

There may also be an optional key file, the name of this file is important, but its content is not. 

- `key_<runName>`

### Key
Within the keys directory
For simple use there may only be one* key. This is a file which does not need to contain anything but must be named to match a cassette  
  
`key_<runName>`

### Run
`./run.sh`  
It is that easy.  
Open and edit the the run script to set the number of cores which snakemake can use. There are multiple run modes available for testing and diagnostics. 

### Notes
Individualized run configuration by runConfig.yaml is not in use at this time, although a single example of its use is available for any users who wish to adapt this. Please use the globalConfig.yaml for adjusting parameters across all pipelines. 

#### barcodes file specifications
per fastGBS V2 barcode specifications. Where each line is  
<barcode><tab><varietyName><newline>  
<varietyName> should contain only simple charcters:letters, numbers, _ and -  
	do not use ;/()<space> or other special characters.  
in a plaintext file named <runName>_barcodes.txt

#### conda environments

To create the snakemake environment required for run.sh  
In the envs directory run  
`conda env create --file snakemake.yaml`  
This can take quite a long time to complete. 

Before calling run.sh activate the snakemake environment with:  
`conda activate snakemake`  
then execute  
`./run.sh`  
All processing conda environment requirements are built into the snakemake rules

### Preparing the reference genome
This GBS pipeline is designed for use on large grain genomes. To reduce memory usage the chromosomes are separated and alignment proceeds individually.  
Decompress the reference fa files with:  
`bgzip -d refgenome.fa.gz`  
Separate the refgenome into individual chromosomes with a bash script  
Index all of the separate chromosomes with:  
`bwa index -a bwtsw refgenome.fa`  
`samtools faidx refgenome.fa`  
place all the files in a single directory  
  
Within each cassette directory the refgenome directory could contain a complete copy of these prepared reference files, however, symlinks also work and are preferred to save disk space. Ensure that the symlink paths are fully specified to allow reuse between many cassettes.  
  
In the config file set the refgenome name settings to exclude the integer, this allows the pipeline processes to sense the number of chromosomes to process.  
for example, if of chromosomes 1-7 you have a chromosome 1 file named  
`SRG_chr1H.fasta`  
then after indexing you will have  

- `SRG_chr1H.fasta`
- `SRG_chr1H.fasta.amb`
- `SRG_chr1H.fasta.ann`
- `SRG_chr1H.fasta.bwt`
- `SRG_chr1H.fasta.fai`
- `SRG_chr1H.fasta.pac`
- `SRG_chr1H.fasta.s`

To exclude the integer set the config settings as:  
`refGenA: SRG_chr`  
`refGenB: H.fasta`

### Advanced use:
*The data handling of cassetteGBS was designed to keep all data and logs fully separated and thus independent, there are variables in the snakefile that activate modes where the single key restriction can be removed, allowing multiple simultaneous runs to be selected by the key files in the keys dir. Another mode allows keyless activation which will simultaneously run every cassette in the cassettes dir. Once a user has advanced familiarity with this pipeline in standard operation they will also have found where these variables are located.  
  
Snakemake will automatically saturate only the cores provided by the run script parameter. With many pipelines, all the cores provided will remain in use for the majority of the run, but various runs will be processing at any stage. This differs from the individual run, which has parallel processing and core saturation but only within a single stage of the run. Be very aware of the memory and disk space and runtime demands when operating in these expanded parallel modes. Extreme demands against worldclass mainframe hardware are possible. The purpose of this software is to redefine GBS at scale. Have fun!


### Order of checkpoints as filed, ...add numeric order of arrival
done.alignReads
done.callVariants
done.cleanSamples
done.copyResults
done.demultiplexReads_0
done.demultiplexReads_1
done.demultiplexReads_2
done.demultiplexReads_3
done.demultiplexReads_4
done.demultiplexReads_5
done.demultiplexReads_6
done.demultiplexReads_7
done.demultiplexReads_all
done.fastQC
done.finalize
done.initialize
done.mergeVariants
done.multiQC
done.prepareMerge
done.removeAdaptors
done.setup
