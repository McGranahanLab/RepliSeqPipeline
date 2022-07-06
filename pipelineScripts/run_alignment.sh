#!/bin/bash

function run_alignment {
	#this function runs the alignment steps of the Repli-Seq pipeline
	#for paired-end sequencing and needs as input a config file generated
	#by generate_config.File.R

    #read in parameters from config file
	source $1

	#number of sequenced lanes
	len=${#FQ1[@]}


	#clip adapters from reads with trim galore (wrapper tool around cutadapt and fastqc)
	jobidsClipping=()
	for (( i=0; i<$len; i++ ))
	do
		name=$SAMPLE.L$((i+1))
		log=$WORKDIR/commandHistory/$name.clip_TrimGalore.log
		jobid=$(sbatch --parsable --time=1-00:00:00 -o $log $SCRIPTDIR/clip_TrimGalore.sh ${FQ1[$i]} ${FQ2[$i]} $ADAPTER $name $CLIPFASTQDIR $FASTQCDIR $NTHREATS) 
		jobidsClipping+="afterok:$jobid,"
	done
	jobidsClipping="${jobidsClipping::-1}"


	#align reads to ref genome
	BAM=()
	jobidsAlignment=()
	for (( i=0; i<$len; i++ ))
	do
		name=$SAMPLE.L$((i+1))
		log=$WORKDIR/commandHistory/$name.alignment.log
		
		if [[ $FQ2 == "NA" ]]
		then 
			jobid=$(sbatch --parsable --time=1-00:00:00 --dependency=$jobidsClipping -o $log -c $NTHREATS --mem ${MEM}G $SCRIPTDIR/align  $CLIPFASTQDIR/${name}_trimmed.fq.gz "NA" $REFGENOME $SCRATCHDIR/BAM/$name $NTHREATS ${MEM}G)
		else
			jobid=$(sbatch --parsable --time=1-00:00:00 --dependency=$jobidsClipping -o $log -c $NTHREATS --mem ${MEM}G $SCRIPTDIR/align  $CLIPFASTQDIR/${name}_val_1.fq.gz $CLIPFASTQDIR/${name}_val_2.fq.gz $REFGENOME $SCRATCHDIR/BAM/$name $NTHREATS ${MEM}G)
		fi
		
		jobidsAlignment+="afterok:$jobid,"
		BAM+=($name.bam)
	done
	jobidsAlignment="${jobidsAlignment::-1}"


	#check stats
	jobidsStats1=()
	for (( i=0; i<$len; i++ ))
		do
			name=$SAMPLE.L$((i+1))
			log=$WORKDIR/commandHistory/$name.samstats.log
			jobid=$(sbatch --parsable --time=3:00:00 --dependency=$jobidsAlignment -o $log $SCRIPTDIR/samstats $SCRATCHDIR/BAM/${BAM[$i]} $SCRATCHDIR/BAM/$name.samstats)
			jobidsStats1+="afterok:$jobid,"
		done
	jobidsStats1="${jobidsStats1::-1}"



	#merge multiple bams if nessesary
	#-> otherwise, move them from scratch into working directory 
	log=$WORKDIR/commandHistory/$SAMPLE.mergeBAMs.log
	jobid=$(sbatch --parsable --time=3:00:00 --dependency=$jobidsStats1 -o $log $SCRIPTDIR/mergeBAMs.sh $SAMPLE $SCRATCHDIR/BAM/ $WORKDIR/BAM/full/ ${BAM[@]})
	jobidMerge="afterok:$jobid"

	
	#filter bams by alignment quality and sort by position
	log=$WORKDIR/commandHistory/$SAMPLE.filtersort.log
	jobid=$(sbatch --parsable --time=3:00:00 --dependency=$jobidMerge -o $log --mem ${MEM}G -c $NTHREATS $SCRIPTDIR/filtersort $WORKDIR/BAM/full/$SAMPLE.bam $WORKDIR/BAM/full/$SAMPLE $NTHREATS ${MEM}G)
	# jobid=$(sbatch --parsable --time=3:00:00 -o $log --mem ${MEM}G -c $NTHREATS $SCRIPTDIR/filtersort $WORKDIR/BAM/full/$SAMPLE.bam $WORKDIR/BAM/full/$SAMPLE $NTHREATS ${MEM}G)
	jobidFilter="afterok:$jobid"
		
	
	#check stats
	log=$WORKDIR/commandHistory/$SAMPLE.samstats.log
	jobid=$(sbatch --parsable --time=3:00:00 --dependency=$jobidFilter -o $log $SCRIPTDIR/samstats $WORKDIR/BAM/full/$SAMPLE.q20_sort.bam $WORKDIR/BAM/full/$SAMPLE.samstats)
	jobidsStats2+="afterok:$jobid"
	
	
	#remove duplicate reads
	log=$WORKDIR/commandHistory/$SAMPLE.dedup.log
	jobid=$(sbatch --parsable --time=3:00:00 --dependency=$jobidsStats2 -o $log $SCRIPTDIR/dedup $WORKDIR/BAM/full/$SAMPLE.q20_sort.bam $WORKDIR/BAM/processed/$SAMPLE.bam)		
	#wait $jobid
	echo $WORKDIR/BAM/processed/$SAMPLE.bam

}







