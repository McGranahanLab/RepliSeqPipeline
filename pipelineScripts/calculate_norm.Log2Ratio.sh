 #!/bin/bash

function calculate_norm.Log2Ratio {
	#this function calculates normalised RPKM scores, using scripts from 
	#https://github.com/4dn-dcic/repli-seq-pipeline. The run_alignment function 
	#needs to be completed successfully. This functions needs to be run on a patient 
	#level and it needs as input a config file for one of the patients samples,
	#minRPKM for filtering and all possible types (eg G1, Early, Late).

	# parameters
	configFile=$1
	source $configFile
	
	# calculate RPKM bedGraphs
	jobidsRPKM=()
	for t in "${TYPES[@]}"
	do
		sample=${PATIENT}_${t}
		bam=$BAMDIR/processed/$sample.bam
		log=$CMDDIR/$sample.count.log
		output=$COUNTDIR/$sample

		jobid=$(sbatch --parsable --time=3:00:00 -o $log $SCRIPTDIR/count $bam $CHROMSIZEFILE $output $WINDOWSIZE)
		jobidsRPKM+="afterok:$jobid,"
	done
	jobidsRPKM="${jobidsRPKM::-1}"

	
	# make a bed file for filtering based on sum of scores across all count bg files and a max rpkm value per type
	log=$CMDDIR/$PATIENT.make_filteredbed.log
	jobid=$(sbatch --parsable --dependency=$jobidsRPKM --time=3:00:00 -o $log $SCRIPTDIR/makeFilterBED.R $configFile)
	# jobid=$(sbatch --parsable --time=3:00:00 -o $log $SCRIPTDIR/makeFilterBED.R $configFile)
	jobidFilterBed="afterok:$jobid"


	# filter windows with a low average RPKM , outliers (if MAXRPKM provided) and windows with early or late RPKM = 0
	jobidsFilterBg=()
	for t in "${TYPES[@]}"
	do
		sample=${PATIENT}_$t
		bg=$COUNTDIR/$sample.w${WINDOWSIZE}.RPKM.bedGraph
		filterBed=$COUNTDIR/${PATIENT}.filtered.windows.bed
		output=$COUNTDIR/$sample
		log=$CMDDIR/${sample}.filterRPKM.log

		jobid=$(sbatch --parsable --dependency=$jobidFilterBed --time=3:00:00 -o $log $SCRIPTDIR/filter $bg $filterBed $output)
		#jobid=$(sbatch --parsable  --time=3:00:00 -o $log $SCRIPTDIR/filter $bg $filterBed $output)
		jobidsFilterBg+="afterok:$jobid,"

	done
	jobidsFilterBg="${jobidsFilterBg::-1}"


	# calculate log2 ratio between Early and Late
	output=$LOG2RATIODIR/$PATIENT
	earlyBg=$COUNTDIR/${PATIENT}_${TYPES[-2]}.filteredRPKM.bedGraph
	lateBg=$COUNTDIR/${PATIENT}_${TYPES[-1]}.filteredRPKM.bedGraph
	log=$CMDDIR/${PATIENT}.log2ratio.log

	jobid=$(sbatch --parsable --dependency=$jobidsFilterBg --time=3:00:00 -o $log $SCRIPTDIR/log2ratio $earlyBg $lateBg $output)
	# jobid=$(sbatch --parsable --time=3:00:00 -o $log $SCRIPTDIR/log2ratio $earlyBg $lateBg $output)
	jobidLog2Ratio="afterok:$jobid"
	
	# normalise log2ratio for coverage, apply quantile normalisation of NORMQUANTILREF != NA and loess-smooth profiles using a 300kb span size
	# --> normalise for coverage to make sure that regions with low coverage have less impact
	log=$WORKDIR/commandHistory/${PATIENT}.smooth.log
	jobid=$(sbatch --parsable --dependency=$jobidLog2Ratio --time=3:00:00 -o $log $SCRIPTDIR/smooth.R $configFile)
	# jobid=$(sbatch --parsable --time=3:00:00 -o $log $SCRIPTDIR/smooth.R $configFile)
	jobidSmooth="afterok:$jobid"

	#identify significant early and late replicated regions with threshold per chromosome
	log=$WORKDIR/commandHistory/${PATIENT}.classifyTiming.log
	jobid=$(sbatch --parsable --dependency=$jobidSmooth --time=3:00:00 -o $log $SCRIPTDIR/classifyTiming.R $configFile)
	# jobid=$(sbatch --parsable --time=3:00:00 -o $log $SCRIPTDIR/classifyTiming.R $configFile)
	jobidClassification="afterok:$jobid"

	# plot results
	log=$CMDDIR/${PATIENT}.plotResults_normLog2Ratio.log
	jobid=$(sbatch --parsable --dependency=$jobidClassification --time=5:00:00 -o $log $SCRIPTDIR/plotResults_normLog2Ratio.R $configFile)
	# jobid=$(sbatch --parsable --time=5:00:00 -o $log $SCRIPTDIR/plotResults_normLog2Ratio.R $configFile)
	jobidPlot="afterok:$jobid"

}




















