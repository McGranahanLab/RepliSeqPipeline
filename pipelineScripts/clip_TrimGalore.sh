#!/bin/bash
fastq1=$1
fastq2=$2      
adapter=$3
basename=$4  
clip_output_dir=$5
fastqc_output_dir=$6
nProcs=$7

#check that input is correct
echo $fastq1
echo $fastq2
echo $adapter
echo $basename
echo $clip_output_dir
echo $fastqc_output_dir
echo $nProcs


#single-end reads
if [[ $fastq2 == 'NA' ]]
then
	if [[ $adapter == 'default' ]]
	then
		/camp/lab/swantonc/working/dietzem/apps/TrimGalore/TrimGalore-0.6.5/trim_galore \
		--fastqc --fastqc_args "--outdir ${fastqc_output_dir}" \
		--cores $nProcs \
		--basename $basename \
		--output_dir $clip_output_dir \
		$fastq1
	else
		/camp/lab/swantonc/working/dietzem/apps/TrimGalore/TrimGalore-0.6.5/trim_galore \
		--fastqc --fastqc_args "--outdir ${fastqc_output_dir}" \
		--cores $nProcs \
		--adapter $adapter \
		--basename $basename \
		--output_dir $clip_output_dir \
		$fastq1
	fi
fi


#paired-end reads
if [[ $fastq2 != 'NA' ]]
then
	if [[ $adapter == 'default' ]]
	then
		/camp/lab/swantonc/working/dietzem/apps/TrimGalore/TrimGalore-0.6.5/trim_galore \
		--fastqc --fastqc_args "--outdir ${fastqc_output_dir}" \
		--cores $nProcs \
		--paired \
		--basename $basename \
		--output_dir $clip_output_dir \
		$fastq1 $fastq2
	else
		/camp/lab/swantonc/working/dietzem/apps/TrimGalore/TrimGalore-0.6.5/trim_galore \
		--fastqc --fastqc_args "--outdir ${fastqc_output_dir}" \
		--cores $nProcs \
		--adapter $adapter \
		--paired \
		--basename $basename \
		--output_dir $clip_output_dir \
		$fastq1 $fastq2
	fi
fi








