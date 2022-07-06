#!/bin/bash
set -o pipefail

SAMPLE=$1
shift
BAMDIR=$1
shift
OUTPUTDIR=$1
shift
BAMS=("$@")

# echo ${BAMS[@]}
# echo $BAMDIR
# echo $OUTPUT

# module load
# ml purge
# ml SAMtools/1.8-foss-2016b

# if there is only one BAM file, then just move file
nBAM=${#BAMS[@]}
if [ $nBAM = 1 ]
	then 
		cp $BAMDIR/$BAMS $OUTPUTDIR/$SAMPLE.bam
		exit
	fi

# create list of input BAM files with one BAM per line
bamList=$BAMDIR/$SAMPLE.bamList.txt
if [[ -f $bamList ]] 
	then 
		rm $bamList 
	fi

touch $bamList
for (( i=0; i<$nBAM; i++ )) 
	do
		echo "$BAMDIR/${BAMS[$i]}" >> $bamList
	done


# merging
samtools merge -n -f -b $bamList $OUTPUTDIR/$SAMPLE.bam



