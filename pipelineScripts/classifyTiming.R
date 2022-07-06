#!/usr/bin/env Rscript

## load docopt package from CRAN
options(scipen=9999)
suppressMessages(library(vroom))
suppressMessages(library(parallel))
suppressMessages(library(GenomicRanges))


### command line arguments ###
cmdArgs        <- commandArgs(trailingOnly = TRUE)
config.file    <- cmdArgs[1]

# cell.line   <- 'HMEC'
# config.file <- paste0('/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/Repli-Seq/Breast/Release_w1000_20200203/', cell.line, '/commandHistory/', cell.line, '_Early.configFile.txt')

#read in paramters
config           <- read.table(config.file, stringsAsFactors = F, sep = '=')
colnames(config) <- c('variable', 'value')
for(i in 1:nrow(config)) {
  assign(config$variable[i], config$value[i])
}


#function to calculate p-values for peaks to classify timing
fun_classifyTiming <- function(log2ratio.file, alpha = 0.05, output.prefix){
  
  #divide data into positive values = early and negative values = late 
  log2ratio       <- as.data.frame(vroom(log2ratio.file, col_names = c('chr', 'start', 'stop', 'score')))
  early.log2ratio <- log2ratio[log2ratio$score > 0,]
  late.log2ratio  <- log2ratio[log2ratio$score < 0,]
  
  #calculate threshold per chromosome
  all_early.log2ratio <- c()
  all_late.log2ratio  <- c()
  thresholds          <- c()
  for(chr in unique(log2ratio$chr)){
    chr_early.log2ratio <- early.log2ratio[early.log2ratio$chr == chr,]
    chr_late.log2ratio  <- late.log2ratio[late.log2ratio$chr == chr,]
    
    threshold.early <- quantile(chr_early.log2ratio$score, alpha, na.rm = T)
    threshold.late  <- quantile(chr_late.log2ratio$score, 1-alpha, na.rm = T)
    thresholds      <- rbind(thresholds, c(threshold.early, threshold.late))
    
    chr_early.log2ratio$timing <- 'early'
    chr_early.log2ratio$timing[chr_early.log2ratio$score < threshold.early] <- 'unknown'
    chr_late.log2ratio$timing <- 'late'
    chr_late.log2ratio$timing[chr_late.log2ratio$score > threshold.late] <- 'unknown'
    
    all_early.log2ratio <- rbind(all_early.log2ratio, chr_early.log2ratio)
    all_late.log2ratio  <- rbind(all_late.log2ratio, chr_late.log2ratio)
  }
  thresholds <- data.frame(chr = unique(log2ratio$chr), early = thresholds[,1], late = thresholds[,2], stringsAsFactors = F)
  
  #combine early and late to one table again
  combined.table       <- rbind(all_early.log2ratio, all_late.log2ratio)
  log2ratio$timing     <- 'unknown'
  ww                   <- match(paste(combined.table[,1], combined.table[,2], sep = ':'), paste(log2ratio[,1], log2ratio[,2], sep = ':'))
  log2ratio$timing[ww] <- combined.table$timing
  log2ratio            <- log2ratio[order(log2ratio$start),]
  order.chr            <- sub('chr', '', log2ratio$chr)
  order.chr            <- sub('X', 23, order.chr)
  order.chr            <- sub('Y', 24, order.chr)
  order.chr            <- order(as.numeric(order.chr))
  log2ratio            <- log2ratio[order.chr,]
  
   
  #identify timing transition regions
  transition.early.to.late <- c()
  transition.late.to.early <- c()
  for(chr in unique(combined.table$chr)){
    sub <- log2ratio[log2ratio$chr == chr,]
    #sub <- sub[order(sub$start),]
    
    early.to.late <- which(sub$score[1:(nrow(sub)-1)] > 0 & sub$score[2:nrow(sub)] < 0)
    first.early   <- min(which(sub$timing == 'early'))
    last.late     <- max(which(sub$timing == 'late'))
    early.to.late <- early.to.late[early.to.late >= first.early & early.to.late <= last.late]
    for(index in early.to.late){
      start <- ifelse(max(which((sub$score[1:(index-1)] - sub$score[2:index] < 0))) == '-Inf',
                      1,
                      max(which((sub$score[1:(index-1)] - sub$score[2:index] < 0))) + 2)
      stop  <- ifelse(min(which((sub$score[(index+1):(nrow(sub)-1)] - sub$score[(index+2):nrow(sub)] < 0))) == 'Inf',
                      nrow(sub),
                      index + min(which((sub$score[(index+1):(nrow(sub)-1)] - sub$score[(index+2):nrow(sub)] < 0))) - 1)
      
      if(stop > start & sub$timing[start] == 'early' & sub$timing[stop] == 'late'){
        transition.early.to.late <- rbind(transition.early.to.late, data.frame(chr, start = sub$start[start], stop = sub$stop[stop], stringsAsFactors = F) )
      }
    }
    
    late.to.early <- which(sub$score[1:(nrow(sub)-1)] < 0 & sub$score[2:nrow(sub)] > 0)
    first.late    <- min(which(sub$timing == 'late'))
    last.early    <- max(which(sub$timing == 'early'))
    late.to.early <- late.to.early[late.to.early >= first.late & late.to.early <= last.early]
    for(index in late.to.early){
      start <- ifelse(max(which((sub$score[1:(index-1)] - sub$score[2:index] > 0))) == '-Inf',
                      1,
                      max(which((sub$score[1:(index-1)] - sub$score[2:index] > 0))) + 2)
      stop  <- ifelse(min(which((sub$score[(index+1):(nrow(sub)-1)] - sub$score[(index+2):nrow(sub)] > 0))) == 'Inf',
                      nrow(sub),
                      index + min(which((sub$score[(index+1):(nrow(sub)-1)] - sub$score[(index+2):nrow(sub)] > 0))) - 1)

      if(stop > start & sub$timing[start] == 'late' & sub$timing[stop] == 'early'){
        transition.late.to.early <- rbind(transition.late.to.early, data.frame(chr, start = sub$start[start], stop = sub$stop[stop], stringsAsFactors = F))
        
      }
    }
  }
  
  
  #update table with transition regions
  log2ratio$timing_withTransition <- log2ratio$timing
  
  log2ratio_gr                <- makeGRangesFromDataFrame(log2ratio)
  transition.early.to.late_gr <- makeGRangesFromDataFrame(transition.early.to.late)
  transition.late.to.early_gr <- makeGRangesFromDataFrame(transition.late.to.early)
  
  overlap.early.to.late <- findOverlaps(transition.early.to.late_gr, log2ratio_gr, minoverlap = 5)
  log2ratio$timing_withTransition[subjectHits(overlap.early.to.late)] <- 'transition_early.to.late'
  overlap.late.to.early <- findOverlaps(transition.late.to.early_gr, log2ratio_gr, minoverlap = 5)
  log2ratio$timing_withTransition[subjectHits(overlap.late.to.early)] <- 'transition_late.to.early'
  
  #save table and threshold
  colnames(log2ratio)  <- sub('score', 'log2ratio', colnames(log2ratio))
  write.table(thresholds, file = paste0(output.prefix, '_thresholds.txt'),row.names = F, quote = F, sep = '\t')
  write.table(log2ratio, file = paste0(output.prefix, '.txt'), row.names = F, quote = F, sep = '\t')
  
}


########    Main   ########
log2ratio.file <- paste0(LOG2RATIODIR, PATIENT, '.loess', SPANSIZE, '.bedGraph')
output.prefix  <- paste0(LOG2RATIODIR, PATIENT, '.timingClassification')
fun_classifyTiming(log2ratio.file, alpha = as.numeric(THRESHOLDL2R), output.prefix)






