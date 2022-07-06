#!/usr/bin/env Rscript

## load docopt package from CRAN
suppressMessages(library(vroom)) 
suppressMessages(library(ggplot2)) 

### command line arguments ###
cmdArgs        <- commandArgs(trailingOnly = TRUE)
config.file    <- cmdArgs[1]

# config.file <- '/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/Repli-Seq/Lung/Release_20210727/50kb/T2Prep//commandHistory/T2Prep_ES.configFile.txt'


#read in paramters
config           <- read.table(config.file, stringsAsFactors = F, sep = '=')
colnames(config) <- c('variable', 'value')
for(i in 1:nrow(config)) {
  assign(config$variable[i], config$value[i])
}
TYPES <- unlist(strsplit(TYPES, ' '))
TYPES <- grep('[(]', TYPES, invert = T, value = T)
TYPES <- grep('[)]', TYPES, invert = T, value = T)

# #####   helper functions  ####
# #function to get average score in bins
# average.bin <- function(data.table, binSize = as.numeric(SPANSIZE)) {
#   out.table <- c()
#   for (chr in paste0('chr', c(1:22, 'X'))) {
#     sub.data <- data.table[data.table$chr == chr,]
#     binStart <- seq(1, max(sub.data$stop), by = binSize)
#     binEnd   <- binStart + binSize - 1
#     
#     bin.data <- data.frame(chr = chr, start = binStart, end = binEnd, stringsAsFactors = F)
#     bin.data$score <- 0
#     
#     for(i in 1:nrow(bin.data)) {
#       ww <- which(sub.data$start < bin.data$end[i] & sub.data$stop >= bin.data$start[i])
#       if(length(ww) == 0) {
#         next
#       }
#       average           <- mean(sub.data$score[ww])
#       bin.data$score[i] <- average
#     }
#     out.table <- rbind(out.table, bin.data)
#   }
#   return(out.table)
# }
#
# # function to adjust data so points with big gabs are not just connected
# adjust_gabs.perChr <- function(plot.data.chr, gab.size = 5000){
#   
#   diff <- plot.data.chr[2:nrow(plot.data.chr),2] - plot.data.chr[1:(nrow(plot.data.chr)-1),3]
#   ww   <- which(diff > gab.size)
#   while(length(ww) != 0) {
#     index <- ww[1]
#     sub   <- plot.data.chr[index:(index+1),]
#     add   <- c(sub[1,1], sub[1,3], sub[2,2], rep(NA, ncol(plot.data.chr)-3))
#     plot.data.chr <- rbind(plot.data.chr[1:index,], add, plot.data.chr[(index+1):nrow(plot.data.chr),])
#     diff          <- as.numeric(plot.data.chr[2:nrow(plot.data.chr),2]) - as.numeric(plot.data.chr[1:(nrow(plot.data.chr)-1),3])
#     ww            <- which(diff > gab.size)
#   }
#   
#   return(plot.data.chr)
#   
# }



#################### Main ######################
#### log2ratio ###
log2ratio.file      <- paste0(LOG2RATIODIR, '/', PATIENT, '.loess', SPANSIZE, '.bedGraph')
log2ratio           <- as.data.frame(vroom(log2ratio.file, col_names = c('chr', 'start', 'stop', 'score')))

#plot
pdf(paste0(PLOTDIR, PATIENT, '.log2ratio.pdf'), width = 10, height = 3)
for(chr in paste0('chr', c(1:22, 'X'))) {
  
  sub       <- log2ratio[log2ratio$chr == chr,]
  sub$start <- sub$start / 1000000
  #sub$log2ratio <- sub$log2ratio * 10
  
  p <- ggplot(sub, aes(x = start, y = score)) +
    geom_line(size = 0.3, col = '#d6604d') + 
    geom_hline(yintercept = 0, col = '#878787') +
    scale_x_continuous(expand = c(0,0)) +
    xlab(paste0('Chromosome ', sub('chr', '', chr), ' (mb)')) + ylab('smoothed log2(early / late)') +
    ggtitle('Log2Ratio between Early and Late replication Timing ') +
    theme_bw() 
  plot(p)
}
dev.off()



### coverage ###
if(length(grep('G1', TYPES))==0){
  early_cov <- as.data.frame(vroom(paste0(COUNTDIR, PATIENT, '_', TYPES[1],'.w', WINDOWSIZE, '.cov'),  col_names = c('chr', 'start', 'stop', 'score')))
  late_cov  <- as.data.frame(vroom(paste0(COUNTDIR, PATIENT, '_',TYPES[2],'.w', WINDOWSIZE, '.cov'), col_names = c('chr', 'start', 'stop', 'score')))

  #plot
  pdf(paste0(PLOTDIR, PATIENT, '.coverage.pdf'))
  for(chr in paste0('chr', c(1:22, 'X'))) {
    data <- rbind(data.frame(early_cov[early_cov[,1]==chr,], Timing = 'Early'),
                  data.frame(late_cov[late_cov[,1]==chr,], Timing = 'Late'))
    data$start <- data$start / 1000000
    data       <- data[data[,4] < quantile(data[,4], 0.95),] #exclude outliers in plot
    
    
    p <- ggplot(data, aes(x = start, y = score, fill = Timing)) + 
      geom_area() +
      facet_grid(Timing ~ .) + 
      theme_bw() + 
      theme(legend.position = 'none') + 
      ylab('coverage') + xlab('Chromosome Position (mb)') + 
      ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
    plot(p)
  }
  dev.off()
} else {
  g1_cov    <- as.data.frame(vroom(paste0(COUNTDIR, PATIENT, '_', TYPES[1],'.w', WINDOWSIZE, '.cov'),  col_names = c('chr', 'start', 'stop', 'score')))
  early_cov <- as.data.frame(vroom(paste0(COUNTDIR, PATIENT, '_', TYPES[2],'.w', WINDOWSIZE, '.cov'),  col_names = c('chr', 'start', 'stop', 'score')))
  late_cov  <- as.data.frame(vroom(paste0(COUNTDIR, PATIENT, '_',TYPES[3],'.w', WINDOWSIZE, '.cov'), col_names = c('chr', 'start', 'stop', 'score')))

  #plot
  pdf(paste0(PLOTDIR, PATIENT, '.coverage.pdf'))
  for(chr in paste0('chr', c(1:22, 'X'))) {
    data <- rbind(data.frame(g1_cov[g1_cov[,1]==chr,], Timing = 'G1'),
                  data.frame(early_cov[early_cov[,1]==chr,], Timing = 'Early'),
                  data.frame(late_cov[late_cov[,1]==chr,], Timing = 'Late'))
    data$start <- data$start / 1000000
    data       <- data[data[,4] < quantile(data[,4], 0.95),] #exclude outliers in plot
    
    p <- ggplot(data, aes(x = start, y = score, fill = Timing)) + 
      geom_area() +
      facet_grid(Timing ~ .) + 
      theme_bw() + 
      theme(legend.position = 'none') + 
      ylab('coverage') + xlab('Chromosome Position (mb)') + 
      ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
    plot(p)
  }
  dev.off()
}


### RPKM counts and log2ratio ###
early_RPKM <- as.data.frame(vroom(paste0(COUNTDIR, PATIENT, '_', TYPES[length(TYPES)-1] ,'.filteredRPKM.bedGraph'),  col_names =  c('chr', 'start', 'stop', 'score')))
late_RPKM  <- as.data.frame(vroom(paste0(COUNTDIR, PATIENT, '_', TYPES[length(TYPES)] ,'.filteredRPKM.bedGraph'),  col_names =  c('chr', 'start', 'stop', 'score')))

#plot per chromosome
pdf(paste0(PLOTDIR, PATIENT, '.RPKMcounts.log2ratio.pdf'), width = 10, height = 3)
for(chr in paste0('chr', c(1:22, 'X'))) {
  sub        <- early_RPKM[early_RPKM[,1]==chr,]
  x.range    <- seq(min(sub$start), max(sub$start), by = as.numeric(WINDOWSIZE))
  full.early <- data.frame(chr = chr, start = x.range, score = 0, Timing = 'Early')
  full.early$score[full.early$start %in% sub$start] <- sub$score
  full.early$start <- full.early$start / 1000000

  sub        <- late_RPKM[late_RPKM[,1]==chr,]
  full.late <- data.frame(chr = chr, start = x.range, score = 0, Timing = 'Late')
  full.late$score[full.late$start %in% sub$start] <- -1*sub$score
  full.late$start <- full.late$start / 1000000

  count.data        <- rbind(full.early, full.late)
  count.data$Timing <- factor(count.data$Timing, levels = c('Early', 'Late'))

  log2ratio.chr       <- log2ratio[log2ratio$chr == chr,]
  log2ratio.chr$start <- log2ratio.chr$start / 1000000
  # log2ratio.chr$score <- log2ratio.chr$score * 10

  
  min.y <- max(min(c(min(log2ratio.chr$score), min(count.data$score))), -10)
  min.y <- ifelse((floor(min.y) %% 2) == 0, floor(min.y), floor(min.y) - 1)
  max.y <- min(max(c(max(log2ratio.chr$score), max(count.data$score))), 10)
  max.y <- ifelse((ceiling(max.y) %% 2) == 0, ceiling(max.y), ceiling(max.y) + 1)
  y.axis.breaks <- c(rev(seq.int(0, min.y, by = -2)), seq.int(0, max.y, by = 2))
  y.axis.breaks <- unique(y.axis.breaks)
  y.axis.labels <- sub('-', '', y.axis.breaks)
  
  p <- ggplot() + 
    geom_area(data = count.data, mapping = aes(x = start, y = score, fill = Timing, group = Timing), position = 'identity') +
    guides(fill=guide_legend(title="RPKM counts")) +
    geom_line(data = log2ratio.chr, mapping = aes(x=log2ratio.chr$start, y = log2ratio.chr$score, colour='black'), size = 0.2) +
    scale_colour_manual(name = '', values = c('black'='black'), labels = c('log2-ratio')) +
    scale_y_continuous(breaks = y.axis.breaks, labels = y.axis.labels, limits = c(min.y, max.y)) +
    theme_bw() + 
    ylab('') + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  plot(p)
}
dev.off()



### Timing Classification and Log2ratio ###
output.dir <- paste0(PLOTDIR, '/classifyTiming_log2ratio/')
if(!dir.exists(output.dir)){
  dir.create(output.dir)
}

classTiming <- as.data.frame(vroom(paste0(LOG2RATIODIR, '/', PATIENT, '.timingClassification.txt')))
thresholds  <- as.data.frame(vroom(paste0(LOG2RATIODIR, '/', PATIENT, '.timingClassification_thresholds.txt')))
pdf(paste0(PLOTDIR, PATIENT, '.classifyTiming_log2ratio.pdf'), width = 13, height = 5)
#png(paste0(PLOTDIR, PATIENT, '.classifyTiming_log2ratio.png'), width = 22, height = 7, units = 'cm', res = 300)
for(chr in paste0('chr', c(1:22, 'X'))){
  
  sub             <- classTiming[classTiming$chr == chr,]
  # sub$log2ratio   <- sub$log2ratio * 10
  sub$xmin        <- sub$start / 1000000
  sub$xmax        <- sub$stop / 1000000
  sub$ymin        <- floor(min(sub$log2ratio))
  sub$ymax        <- ceiling(max(sub$log2ratio))
  sub$classTiming <- sub('_early.to.late', '', sub$timing)
  sub$classTiming <- sub('_late.to.early', '', sub$classTiming)
  sub$classTiming <- factor(sub$classTiming, levels = c('early', 'late', 'transition', 'unknown'))
  plot.thresholds <- as.numeric(thresholds[thresholds$chr == chr, c('early', 'late')]) * 10
 
  p <- ggplot() + 
    geom_rect(data = sub, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = sub$classTiming)) +
    scale_fill_manual(name = 'Timing', values = c('early'= '#de77ae', 'late' = '#7fbc41', 'transition' = '#fee0b6', 'unknown' = '#f7f7f7'), labels = c('early', 'late', 'transition', 'unknown')) + 
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = '#4d4d4d', size = 0.2) + 
    geom_hline(yintercept = as.numeric(plot.thresholds[1]), linetype = 'dashed', color = '#1a1a1a', size = 0.3) +
    geom_hline(yintercept = as.numeric(plot.thresholds[2]), linetype = 'dashed', color = '#1a1a1a', size = 0.3) +
    geom_path(data = sub, aes(x = xmin, y = log2ratio), size = 0.3) +
    theme_bw() + theme(legend.key = element_rect(colour = "black")) +
    ylab('smoothed Log2ratio') + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  # png(paste0(output.dir, PATIENT,'.', chr, '.classifyTiming_log2ratio.png'), width = 22, height = 7, units = 'cm', res = 300)
  plot(p)
  # dev.off()
}
dev.off()










