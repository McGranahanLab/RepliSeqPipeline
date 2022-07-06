#!/bin/bash 

#SBATCH --mem=6G
#SBATCH --partition=cpu
#SBATCH --job-name=RepliSeq
#SBATCH --time=1-00:00:00
#SBATCH -e log_dir/_RepliSeq.wrapper.%j.stderr
#SBATCH -o log_dir/RepliSeq.wrapper.%j.stdout

# load modules
ml purge
ml cutadapt/1.18-foss-2018b-Python-3.6.6
ml FastQC/0.11.8-Java-1.8
ml R/3.6.3-foss-2020a
ml SAMtools/1.8-foss-2016b
ml BWA/0.7.17-foss-2018b
ml BEDTools/2.26.0-foss-2016b
ml Anaconda2/4.2.0

# parameters
designFile="./examples/designFile.txt"       #input path to design file 
releaseVersion="release_version"             #set name for release version
workDir="output_dir"                         #set output directory
scratchDir="scratich_dir"                    #set scratch directory
scriptDir="./pipelineScripts/"
refGenome="refGenome_dir/ucsc.hg19.fasta"    #needs to be downloaded before running this script and the total path to the reference needs to be added here
chromSize="refGenome_dir/hg19.chr.sizes.txt" #needs to be downloaded before running this script and the total path to the chromsome size file needs to be added here
adapter=default                              #set a adapter sequence for trimming here otherwise default sequence from timGalore will be used
windowSize=1000
minRPKM=0.1
maxRPKM='NA'
spanSize=300000
refCellLine="H1650"                          #choose cell line that should be used as reference for qunatile normalisation
qnormRef="$workDir/$releaseVersion/$refCellLine/Log2Ratio/$refCellLine.l2r.bedGraph"
thresholdL2R=0.03
logRepliscan=FALSE
nProcs=4
memRequest=16

configList=$workDir/$releaseVersion/configFile.list.txt


# create config files if they don't exist yet
if [[ ! -f $configList ]];
	then
		Rscript $scriptDir/generate_configFile.R $designFile $releaseVersion $workDir $scratchDir $scriptDir $refGenome $chromSize $adapter $windowSize $minRPKM $maxRPKM $spanSize $qnormRef $thresholdL2R $logRepliscan $nProcs $memRequest
		while [ ! -f $configList ];
		do
			sleep 60;
		done
	fi


# alignment
source $scriptDir/run_alignment.sh

BAMs=()
while read i; 
do
  echo "$i"
  bam=$(run_alignment $i) 
  BAMs+=( $bam )
done <$configList


# calculate normalised RPKM scores and Log2Ratio
source $scriptDir/calculate_norm.Log2Ratio.sh

PATIENTS=( $(cat $designFile | cut -f 1) )
PATIENTS=${PATIENTS[@]:1}
PATIENTS=( $(echo "${PATIENTS[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )

for p in ${PATIENTS[@]}
do
	configFile=( $(grep "$p" $configList ) )
	calculate_norm.Log2Ratio ${configFile[1]}
done


# RepliScan
source $scriptDir/run_repliscan.sh
PATIENTS=( $(cat $designFile | cut -f 1) )
PATIENTS=${PATIENTS[@]:1}
PATIENTS=( $(echo "${PATIENTS[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )

for p in ${PATIENTS[@]}
do
	log=$workDir/$releaseVersion/$p/commandHistory/${p}.repliscan.log
	configFile=( $(grep "$p" $configList ) )
	jobid=$(sbatch --parsable --time=3:00:00 -o $log $scriptDir/run_repliscan.sh ${configFile[0]})
	run_repliscan ${configFile[1]} &> $log &
done


# compare results from pipelines
PATIENTS=( $(cat $designFile | cut -f 1) )
PATIENTS=${PATIENTS[@]:1}
PATIENTS=( $(echo "${PATIENTS[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )

for p in ${PATIENTS[@]}
do
	log=$workDir/$releaseVersion/$p/commandHistory/${p}.comparePipelines.log
	configFile=( $(grep "$p" $configList ) )
	jobid=$(sbatch --parsable --time=1:00:00 -o $log $scriptDir/comparePipelineResults.R $configFile)
done


# combine log2ratio and repliscan timing classification
PATIENTS=( $(cat $designFile | cut -f 1) )
PATIENTS=${PATIENTS[@]:1}
PATIENTS=( $(echo "${PATIENTS[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )

for p in ${PATIENTS[@]}
do
	log=$workDir/$releaseVersion/$p/commandHistory/${p}.combineTiming.log
	configFile=( $(grep "$p" $configList ) )
	jobid=$(sbatch --parsable --time=1:00:00 -o $log $scriptDir/combineTiming.R $configFile)
done




 


