#!/bin/bash

#function run_repliscan {
	#this function runs the python scripts for the RepliScan analysis 
	#(https://github.com/zyndagj/repliscan).
	#As input it needs the reference fasta file and a text file with the 
	#bam files for the different time points (G1, ES, ...).

    #parameters
	configFile=$1
	source $configFile

	create input folder
	inputDir=$REPLISCANDIR/input/
	if [[ -d $inputDir ]]
	then
		rm -r $inputDir
	fi
	mkdir -m 770 $inputDir

	#create input txt file and symlinks
	bamFolder=$BAMDIR/processed/
	ln -s $bamFolder/${PATIENT}_${TYPES[-2]}.bam $inputDir/${PATIENT}_${TYPES[-2]}.bam
	ln -s $bamFolder/${PATIENT}_${TYPES[-1]}.bam $inputDir/${PATIENT}_${TYPES[-1]}.bam

    # #uncomment following code if G1 should be included
	# if [[ ${TYPES[@]} =~ "G1" ]]
	# then
	# 	ln -s $bamFolder/${PATIENT}_G1.bam $inputDir/${PATIENT}_G1.bam
	# 	echo -e "G1\t$inputDir/${PATIENT}_G1.bam" > $inputDir/${PATIENT}_input.txt
	# else
	# 	echo -e "G1\t$inputDir/${PATIENT}_Early.bam\t$inputDir/${PATIENT}_Late.bam" > $inputDir/${PATIENT}_input.txt
	# fi

	echo -e "G1\t$inputDir/${PATIENT}_${TYPES[-2]}.bam\t$inputDir/${PATIENT}_${TYPES[-1]}.bam" > $inputDir/${PATIENT}_input.txt
	echo -e "ES\t$inputDir/${PATIENT}_${TYPES[-2]}.bam" >> $inputDir/${PATIENT}_input.txt
	echo -e "LS\t$inputDir/${PATIENT}_${TYPES[-1]}.bam" >> $inputDir/${PATIENT}_input.txt
	

	#run repliscan
	cd $REPLISCANDIR
	if [[ $LOGREPLISCAN == TRUE ]]
	then 
		repliscan.py -r $REFGENOME -w $WINDOWSIZE --log $inputDir/${PATIENT}_input.txt
	else
		repliscan.py -r $REFGENOME -w $WINDOWSIZE $inputDir/${PATIENT}_input.txt
	fi

	#plot results
	log=$CMDDIR/${PATIENT}.plotResults_repliscan.log
	jobid=$(sbatch --parsable --time=5:00:00 -o $log $SCRIPTDIR/plotResults_repliscan.R $configFile)
	#jobidSmooth="afterok:$jobid"
	
#}
