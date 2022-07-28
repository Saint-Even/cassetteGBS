#!/bin/bash
echo "Begin: run.sh"

#Run options
#set type to: run, testRun, produceDAG, showSteps
type=run
cores=6

echo "Run Mode is set to: ${type}"
echo "Snakemake will use up to: ${cores} cores"


if [ ${type} == "run" ]; then
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

