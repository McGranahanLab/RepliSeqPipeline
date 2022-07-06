#!/usr/bin/env Rscript

## load docopt package from CRAN
suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN
suppressMessages(library(ggplot2)) 
suppressMessages(library(grid)) 
suppressMessages(library(gridExtra)) 

## configuration for docopt
doc <- "Usage: plotResults_repliscan [-h] CONFIGFILE

-h --help                   show this help text"

## docopt parsing
opt           <- docopt(doc)
config.file   <- opt$CONFIGFILE

#config.file <- '/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/Repli-Seq/ENCODE/Lung_release_w5000/A549/commandHistory/A549_Early.configFile.txt'


#read in paramters
config           <- read.table(config.file, stringsAsFactors = F, sep = '=')
colnames(config) <- c('variable', 'value')
for(i in 1:nrow(config)) {
  assign(config$variable[i], config$value[i])
}


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





#################### Main ######################
#### smooth wavelet signals ###
if(LOGREPLISCAN == TRUE){
  early_smoothSignal.file <- paste0(REPLISCANDIR, '/ES_logFC_3.smooth.bedgraph')
  late_smoothSignal.file  <- paste0(REPLISCANDIR, '/LS_logFC_3.smooth.bedgraph')
} else {
  early_smoothSignal.file <- paste0(REPLISCANDIR, '/ES_ratio_3.smooth.bedgraph')
  late_smoothSignal.file  <- paste0(REPLISCANDIR, '/LS_ratio_3.smooth.bedgraph')
}
early_smoothSignal           <- read.table(early_smoothSignal.file, header = F, stringsAsFactors = F)
colnames(early_smoothSignal) <- c('chr', 'start', 'stop', 'score')
early_smoothSignal           <- early_smoothSignal[early_smoothSignal$chr %in% paste0('chr', c(1:22, 'X', 'Y')),]

late_smoothSignal           <- read.table(late_smoothSignal.file, header = F, stringsAsFactors = F)
colnames(late_smoothSignal) <- c('chr', 'start', 'stop', 'score')
late_smoothSignal           <- late_smoothSignal[late_smoothSignal$chr %in% paste0('chr', c(1:22, 'X', 'Y')),]

#create data.frame for plotting
early.data        <- data.frame(chr = early_smoothSignal[,1], y = early_smoothSignal[,4], stringsAsFactors = F)
early.data$x      <- ((early_smoothSignal[,3] - early_smoothSignal[,2]) / 2) + early_smoothSignal[,2]
early.data$Timing <- 'Early'

late.data        <- data.frame(chr = late_smoothSignal[,1], y = late_smoothSignal[,4], stringsAsFactors = F)
late.data$x      <- ((late_smoothSignal[,3] - late_smoothSignal[,2]) / 2) + late_smoothSignal[,2]
late.data$Timing <- 'Late'

plot.data <- rbind(early.data, late.data)

if(LOGREPLISCAN == TRUE){
  ww  <- match(paste(early.data$chr, early.data$x, sep = ':'), paste(late.data$chr, late.data$x, sep = ':'))
  ratio.plot.data   <- early.data[,1:3]
  ratio.plot.data$y <- early.data$y - late.data$y[ww]
} else {
  ww  <- match(paste(early.data$chr, early.data$x, sep = ':'), paste(late.data$chr, late.data$x, sep = ':'))
  ratio.plot.data   <- early.data[,1:3]
  ratio.plot.data$y <- log(early.data$y / late.data$y[ww])
}

#plot
pdf(paste0(PLOTDIR, PATIENT, '.smoothSignals_repliscan.pdf'), width = 15, height = 10)
for(chr in paste0('chr', c(1:22, 'X'))){
  p.single <- ggplot(plot.data[plot.data$chr == chr,], aes(x = x, y = y, color = Timing)) +
    geom_line(size = 0.5) + 
    facet_grid(Timing ~ .) + 
    xlab('') + ylab(ifelse(LOGREPLISCAN == TRUE, 'log(smooth wavelet signal)', 'smooth wavelet signal')) +
    ggtitle('Wavelet Smoothed Signals (Repliscan)') +
    theme_bw() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(strip.text.y = element_blank())
  
  p.combined <- ggplot(plot.data[plot.data$chr == chr,], aes(x = x, y = y, color = Timing, group = Timing)) +
    geom_line(size = 0.5) + 
    xlab('') + ylab(ifelse(LOGREPLISCAN == TRUE, 'log(smooth wavelet signal)', 'smooth wavelet signal')) +
    theme_bw() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(legend.position = 'none')
  
  p.ratio <- ggplot(ratio.plot.data[ratio.plot.data$chr == chr,], aes(x = x, y = y)) +
    geom_line(size = 0.5, color = '#de77ae') + 
    xlab(paste0('Chromosome ', sub('chr', '', chr))) + ylab(ifelse(LOGREPLISCAN == TRUE, 'log(early) - log(late)', 'log(early / late)')) +
    theme_bw() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
  
  #combine plots
  p.single    <- ggplotGrob(p.single + theme(plot.margin = unit(c(1, 0.5, 0.1, 0.5), "cm")))
  p.combined  <- ggplotGrob(p.combined + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
  p.ratio     <- ggplotGrob(p.ratio + theme(plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm")))
  
  p.single$widths   <- unit.pmax(p.single$widths, p.combined$widths, p.ratio$widths)
  p.combined$widths <- unit.pmax(p.single$widths, p.combined$widths, p.ratio$widths)
  p.ratio$widths    <- unit.pmax(p.single$widths, p.combined$widths, p.ratio$widths)

  combined.list        <- list(p.single, p.combined, p.ratio)
  combined.plot.layout <- matrix(c(rep(1, 2), rep(2,1), rep(3,1)), ncol = 1, byrow = T)
  combined.plot        <- arrangeGrob(grobs = combined.list, nrow = dim(combined.plot.layout)[1], ncol = 1, layout_matrix = combined.plot.layout)
  
  grid.draw(combined.plot)
  grid.newpage()
}

dev.off()



### Coverage ###
g1_cov    <- read.table(paste0(REPLISCANDIR, 'G1.bedgraph'),  header = F, stringsAsFactors = F)
early_cov <- read.table(paste0(REPLISCANDIR, 'ES.bedgraph'),  header = F, stringsAsFactors = F)
late_cov  <- read.table(paste0(REPLISCANDIR, 'LS.bedgraph'),  header = F, stringsAsFactors = F)
colnames(g1_cov) <- colnames(early_cov) <- colnames(late_cov) <- c('chr', 'start', 'stop', 'score')

#plot 
pdf(paste0(PLOTDIR, PATIENT, '.coverage_repliscan.pdf'), width = 10, height = 5)
for(chr in paste0('chr', c(1:22, 'X'))) {
  
  data <- rbind(data.frame(g1_cov[g1_cov[,1]==chr,], timing = 'G1'),
                data.frame(early_cov[early_cov[,1]==chr,], timing = 'Early'),
                data.frame(late_cov[late_cov[,1]==chr,], timing = 'Late'))
  data$start <- data$start / 1000000
  
  p <- ggplot(data, aes(x = start, y = score, fill = timing)) + 
    geom_area() +
    facet_grid(timing ~ .) + 
    theme_bw() + 
    ylab('coverage') + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  plot(p)
}
dev.off()




### Normalized aggregated signal from replicates per Chromosome ###
g1_counts    <- read.table(paste0(REPLISCANDIR, 'G1_norm.bedgraph'),  header = F, stringsAsFactors = F)
early_counts <- read.table(paste0(REPLISCANDIR, 'ES_norm.bedgraph'),  header = F, stringsAsFactors = F)
late_counts  <- read.table(paste0(REPLISCANDIR, 'LS_norm.bedgraph'),  header = F, stringsAsFactors = F)
colnames(g1_counts) <- colnames(early_counts) <- colnames(late_counts) <- c('chr', 'start', 'stop', 'score')

# g1_counts.bin    <- average.bin(g1_counts)
# early_counts.bin <- average.bin(early_counts)
# late_counts.bin  <- average.bin(late_counts)
# ratio.bin        <- average.bin(repTimingProfile)

#plot 
pdf(paste0(PLOTDIR, PATIENT, '.normCounts_repliscan.pdf'), width = 10, height = 5)
for(chr in paste0('chr', c(1:22, 'X'))) {
  
  data <- rbind(data.frame(g1_counts[g1_counts[,1]==chr,], timing = 'G1'),
                data.frame(early_counts[early_counts[,1]==chr,], timing = 'Early'),
                data.frame(late_counts[late_counts[,1]==chr,], timing = 'Late'))
  data$start <- data$start / 1000000
  
  # if(LOGREPLISCAN == TRUE) {
  #   data$score <- log(data$score)
  # }
  # 
  # ratio.chr <- rbind(data.frame(early_repProfile[early_repProfile[,1] == chr,1:3], score = NA, timing = 'G1'),
  #                    data.frame(early_repProfile[early_repProfile[,1] == chr,], timing = 'Early'),
  #                    data.frame(late_repProfile[late_repProfile[,1] == chr,], timing = 'Late'))
  # ratio.chr$start <- ratio.chr$start / 1000000
  
  p <- ggplot(data, aes(x = start, y = score, fill = timing)) + 
    geom_area() +
    facet_grid(timing ~ .) + 
    #geom_path(data = ratio.chr, mapping = aes(x=ratio.chr$start, y = ratio.chr$score, colour='black')) +
    #scale_colour_manual(name = '', values = c('black'='black'), labels = c('smoothed log ratio signal')) +
    theme_bw() + 
    ylab('counts') + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  plot(p)
}
dev.off()




### Smoothed ratio signal with replication timing segmentation ###
if(LOGREPLISCAN == TRUE) {
  segmentation <- read.table(paste0(REPLISCANDIR, 'logFC_segmentation.gff3'), header = F, sep = '\t', stringsAsFactors = F)
} else {
  segmentation <- read.table(paste0(REPLISCANDIR, 'ratio_segmentation.gff3'), header = F, sep = '\t', stringsAsFactors = F)
}

colour.timing           <- segmentation[,c(1,4,5)]
colnames(colour.timing) <- c('chr', 'start', 'stop')
colour.timing$timing    <- matrix(unlist(strsplit(segmentation[,9], '=')), ncol = 3, byrow = T)[,3]
colour.timing$timing    <- sub(';color', '', colour.timing$timing)
colour.timing$timing    <- factor(colour.timing$timing, levels = c('ES', 'ESLS', 'LS'))

plot.data <- rbind(early.data, late.data)

#plot 
pdf(paste0(PLOTDIR, PATIENT, '.smoothSignals.segmentation_repliscan.pdf'), width = 13, height = 5)
for(chr in paste0('chr', c(1:22, 'X'))) {
  
  data       <- data.frame(plot.data[plot.data[,1]==chr,])
  data$x     <- data$x / 1000000
  
  sub.colour       <- colour.timing[colour.timing$chr==chr,]
  sub.colour$xmin  <- sub.colour$start / 1000000
  sub.colour$xmax  <- sub.colour$stop / 1000000
  sub.colour$ymin  <- floor(min(data$y))
  sub.colour$ymax  <- ceiling(max(data$y))
  sub.colour       <- rbind(data.frame(sub.colour, Timing = 'Early'),
                            data.frame(sub.colour, Timing = 'Late'))
                            
  p <- ggplot() + 
    geom_rect(data = sub.colour, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = sub.colour$timing)) +
    scale_fill_manual(name = 'Timing', values = c('ES'='#de77ae', 'ESLS' = '#f7f7f7', 'LS' = '#7fbc41'), labels = c('ES', 'ESLS', 'LS')) + 
    geom_path(data = data, aes(x = x, y = y), size = 0.3) +
    facet_grid(Timing ~ .) +
    theme_bw() + theme(legend.key = element_rect(colour = "black")) +
    ylab(ifelse(LOGREPLISCAN == TRUE, 'log(smooth wavlete signal)', 'smooth wavelet signal')) + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  plot(p)
}
dev.off()








