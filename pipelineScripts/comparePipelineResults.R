#!/usr/bin/env Rscript

###############
### Script to compare log2ratio from 4DN pipeline with Repliscan results
###############

## load docopt package from CRAN
suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN
suppressMessages(library(ggplot2)) 
suppressMessages(library(data.table)) 
suppressMessages(library(BSgenome)) 
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer)) 



## configuration for docopt
doc <- "Usage: comparePipelineResults [-h] CONFIGFILE

-h --help                   show this help text"

## docopt parsing
opt           <- docopt(doc)
config.file   <- opt$CONFIGFILE

#config.file <- '/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/Repli-Seq/ENCODE/Release_w1000_20200106/A549/commandHistory/A549_Early.configFile.txt'


#read in paramters
config           <- read.table(config.file, stringsAsFactors = F, sep = '=')
colnames(config) <- c('variable', 'value')
for(i in 1:nrow(config)) {
  assign(config$variable[i], config$value[i])
}

### Snoothed Values ###
l2r             <- as.data.frame(fread(paste0(LOG2RATIODIR, PATIENT, '.loess300000.bedGraph'), col.names = c('chr', 'start', 'stop', 'value'), stringsAsFactors = F))
early.repliscan <- as.data.frame(fread(paste0(REPLISCANDIR, 'ES_ratio.bedgraph'), col.names = c('chr', 'start', 'stop', 'value'), stringsAsFactors = F))
late.repliscan  <- as.data.frame(fread(paste0(REPLISCANDIR, 'LS_ratio.bedgraph'), col.names = c('chr', 'start', 'stop', 'value'), stringsAsFactors = F))

early.repliscan <- early.repliscan[early.repliscan[,1]%in%paste0('chr', c(1:22, 'X', 'Y')),]
late.repliscan  <- late.repliscan[late.repliscan[,1]%in%paste0('chr', c(1:22, 'X', 'Y')),]
l2r$value       <- l2r$value 

#plot
if(LOGREPLISCAN == FALSE) {
  plot.data <- rbind(data.frame(early.repliscan[,c(1,2,3)], value = early.repliscan[,4], Type = 'repliscan_Early'),
                     data.frame(late.repliscan[,c(1,2,3)],  value = -1*late.repliscan[,4], Type = 'repliscan_Late'),
                     data.frame(l2r, Type= 'Log2Ratio'))
  plot.data <- plot.data[plot.data$value != 0,]
} else {
  plot.data <- rbind(data.frame(early.repliscan[,c(1,2,3)], value = early.repliscan[,4], Type = 'log(repliscan_Early)'),
                     data.frame(late.repliscan[,c(1,2,3)],  value = -1*late.repliscan[,4], Type = 'log(repliscan_Late)'),
                     data.frame(l2r, Type= 'Log2Ratio'))
}

pdf(paste0(PLOTDIR, PATIENT, '.comparison.smoothedValues.pdf'), width = 10, height = 3)
for (chr in paste0('chr', c(1:22, 'X', 'Y'))){
  p <- ggplot(plot.data[plot.data[,1]==chr,], aes(x = start / 1000000, y = value, group = Type, colour = Type)) + 
    geom_line(size = 0.3) +
    scale_color_manual(values = c('#f1b6da', '#b8e186', '#4d4d4d')) +
    xlab(paste0('Chromosome ', sub('chr', '', chr), ' (mb)')) + ylab('value') +
    theme_bw()
  plot(p)
}
dev.off()




### Segmentation ###
  #read data in
  l2r.classTiming                 <- as.data.frame(fread(paste0(LOG2RATIODIR, PATIENT, '.timingClassification.txt'), stringsAsFactors = F))
  l2r.classTiming$timing          <- sub('_early.to.late', '', l2r.classTiming$timing)
  l2r.classTiming$timing          <- sub('_late.to.early', '', l2r.classTiming$timing)
  repliscan.segmentation          <- as.data.frame(fread(paste0(REPLISCANDIR, grep('segmentation', list.files(REPLISCANDIR), value = T)), stringsAsFactors = F))
  repliscan.classTiming           <- repliscan.segmentation[,c(1,4,5)]
  colnames(repliscan.classTiming) <- c('chr', 'start', 'stop')
  repliscan.classTiming$timing    <- matrix(unlist(strsplit(repliscan.segmentation[,9], '=')), ncol = 4, byrow = T)[,3]
  repliscan.classTiming$timing    <- sub(';color', '', repliscan.classTiming$timing)
  repliscan.classTiming$timing    <- sub('ES$', 'early', repliscan.classTiming$timing)
  repliscan.classTiming$timing    <- sub('^LS', 'late', repliscan.classTiming$timing)
  repliscan.classTiming$timing    <- sub('ESLS', 'transition', repliscan.classTiming$timing)
  repliscan.classTiming           <- repliscan.classTiming[repliscan.classTiming$chr %in% paste0('chr', c(1:22, 'X', 'Y')),]

  #create combined table
  gr.l2r.classTiming       <- GRanges(seqnames = l2r.classTiming$chr, IRanges(start = l2r.classTiming$start, end = l2r.classTiming$stop))
  gr.repliscan.classTiming <- GRanges(seqnames = repliscan.classTiming$chr, IRanges(start = repliscan.classTiming$start, end = repliscan.classTiming$stop))
  overlap                  <- findOverlaps(gr.repliscan.classTiming, gr.l2r.classTiming)

  combined.classTiming                  <- l2r.classTiming[,c(1,2,3,5)]
  colnames(combined.classTiming)[4]     <- 'l2r.timing'
  combined.classTiming$repliscan.timing <- 'unknown'
  combined.classTiming$repliscan.timing[subjectHits(overlap)] <- repliscan.classTiming$timing[queryHits(overlap)]

  #create pie chart with fraction concordant and discordant
  ww                   <- which(combined.classTiming$l2r.timing != combined.classTiming$repliscan.timing)
  discordant.table     <- combined.classTiming[ww,]
  discordant.table     <- discordant.table[!(discordant.table$l2r.timing %in% c('early', 'late') & discordant.table$repliscan.timing == 'ESLS'),]
  discordant <- round(nrow(discordant.table) / nrow(combined.classTiming) * 100)
  concordant <- 100 - discordant
  plot.data  <- data.frame(value = c(discordant, concordant), Type = c('discordant', 'concordant'))
  plot.data  <- plot.data %>% 
    arrange(desc(Type)) %>%
    mutate(ypos = cumsum(value) - 0.5*value )
  
  p <- ggplot(plot.data, aes(x = "", y = value, fill = Type)) +
    geom_bar(width = 1, stat = "identity", color = '#525252') + 
    coord_polar("y", start = 0) + 
    geom_text(aes(y = ypos, label = paste0(value, '%')), color = "black", size = 4) +
    scale_fill_manual(name = 'Type', values = c('#bababa', '#fdb863')) +
    xlab("") + ylab("") + ggtitle("Fraction of concordance between\nrepliscan and log2ratio approach") +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  pdf(paste0(PLOTDIR, PATIENT, '.comparison.classificationTiming.pdf'), width = 6, height = 4)
  plot(p)
  
  #create bar plot with different types of discordant cases
  discordant.table$key <- paste0('l2r_', discordant.table$l2r.timing, ':repliscan_', discordant.table$repliscan.timing)
  plot.data                  <- data.frame(value = round(as.numeric(table(discordant.table$key)) / nrow(discordant.table) * 100), 
                                           Difference = names(table(discordant.table$key)))
  plot.data  <- plot.data %>% 
    arrange(desc(Difference)) %>%
    mutate(ypos = cumsum(value) - 0.5*value )
  
  p <- ggplot(plot.data, aes(x = "", y = value, fill = Difference)) +
    geom_bar(width = 1, stat = "identity", color = '#525252') + 
    coord_polar("y", start = 0) + 
    geom_text(aes(y = ypos, label = paste0(value, '%')), color = "black", size = 4) +
    scale_fill_manual(name = 'Type', values = brewer.pal(n = nrow(plot.data), name = 'Set3')) +
    xlab("") + ylab("") + ggtitle("Fractions of different discordant types") +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  plot(p) 
  dev.off()
  









  






