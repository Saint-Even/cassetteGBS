#!/bin/bash
echo "Begin: run.sh"

#Run options
#set type to: run, clusterRun, testRun, produceDAG, showSteps, reRunInc
type=clusterRun
cores=24

#cluster qsub settings
jobs=32

echo "Run Mode is set to: ${type}"
echo "To continue enter y: "
read cont
if [ ${cont} != "y" ]; then
    echo "Exiting the run script. Edit the script to modify modes"
    exit
fi

if [ ${type} == "run" ]; then
    echo "Running..."
    snakemake \
        --cores ${cores} \
        --printshellcmds \
	--keep-going \
        --use-conda \
        --conda-frontend conda \
        all \
	> logs/output.log \
	2> logs/output2.log
fi

if [ ${type} == "clusterRun" ]; then
    echo "Running..."
    snakemake \
	--cluster "qsub" \
	--jobs ${jobs} \
	--cores \
	--latency-wait 300 \
        --printshellcmds \
	--keep-going \
        --use-conda \
        --conda-frontend conda \
        all \
	> logs/output.log \
	2> logs/output2.log
fi

if [ ${type} == "reRunInc" ]; then
    snakemake \
        --cores ${cores} \
        --printshellcmds \
	--keep-going \
	--rerun-incomplete \
        --use-conda \
        --conda-frontend conda \
        all \
	> logs/output.log \
	2> logs/output2.log
fi

if [ ${type} == "testRun" ]; then
	snakemake \
		--forceall \
		--cores ${cores} \
		--printshellcmds \
		--use-conda \
		--conda-frontend conda \
		all
fi

if [ ${type} == "produceDAG" ]; then
	 snakemake \
		--forceall \
		--dag \
		all | dot -Tpdf > pipelineDAG.pdf
fi

if [ ${type} == "showSteps" ]; then
	snakemake \
		--forceall \
		--dry-run \
		--reason \
		--cores ${cores} \
		--use-conda \
		--conda-frontend conda \
		all
fi

