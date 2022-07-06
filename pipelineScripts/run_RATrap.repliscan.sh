#!/bin/bash

#Comparing replication between experiments

#Input: 1) First Segmentation Profile (Repliscan GFF3)
#       2) Second Segmentation Profile (Repliscan GFF3)
#       3) Output
#       4) Reference FASTA
#       5) Times (e.g. ES,MS,LS)
#       6) Minimum distance to be RAT (e.g. 0.5)
#       7) Tile Size (e.g. 1000)

#load moduled
ml Anaconda2/4.2.0

#Input
SegA=$1
SegB=$2
Output=$3
OutputDir=$4
Ref=$5
Times=$6
minDist=$7
tileSize=$8

#run comparison
cd $OutputDir
RATrap.py -d $minDist -S $tileSize -A $SegA -B $SegB -T $Times -F $Ref -O $Output --stats

#merges output from RATrap.py to simplify interpretation
mergeRATs.py --same  $Output > mergeRAT.results.txt

