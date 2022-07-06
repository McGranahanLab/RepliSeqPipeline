#!/usr/bin/env Rscript

################################################################################################################################
############                 Calculate combined replication timing signal and segmentation                         ############   
################################################################################################################################

## load docopt package from CRAN
suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN
suppressMessages(library(ggplot2)) 
suppressMessages(library(grid)) 
suppressMessages(library(gridExtra)) 
suppressMessages(library(data.table)) 
suppressMessages(library(GenomicRanges)) 

## configuration for docopt
doc <- "Usage: combineResults [-h] CONFIGFILE

-h --help                   show this help text"

## docopt parsing
opt           <- docopt(doc)
config.file   <- opt$CONFIGFILE

# cell.line   <- 'H1650'
# config.file <- paste0('/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/Repli-Seq/Lung/Release_w1000_20200203/', cell.line, '/commandHistory/', cell.line, '_Early.configFile.txt')

#read in paramters
config           <- read.table(config.file, stringsAsFactors = F, sep = '=')
colnames(config) <- c('variable', 'value')
for(i in 1:nrow(config)) {
  assign(config$variable[i], config$value[i])
}


### read in Log2-ratio and Repliscan results ###
l2r.classTiming         <- as.data.frame(fread(paste0(LOG2RATIODIR, PATIENT, '.timingClassification.txt')))
l2r.classTiming$timing_withTransition  <- sub('_early.to.late', '', l2r.classTiming$timing_withTransition)
l2r.classTiming$timing_withTransition  <- sub('_late.to.early', '', l2r.classTiming$timing_withTransition)
if(LOGREPLISCAN == TRUE) {
  repliscan.segmentation <- as.data.frame(fread(paste0(REPLISCANDIR, 'logFC_segmentation.gff3')))
} else {
  repliscan.segmentation <- as.data.frame(fread(paste0(REPLISCANDIR, 'ratio_segmentation.gff3')))
}
repliscan.classTiming           <- repliscan.segmentation[,c(1,4,5)]
colnames(repliscan.classTiming) <- c('chr', 'start', 'stop')
repliscan.classTiming$timing    <- matrix(unlist(strsplit(repliscan.segmentation[,9], '=')), ncol = 4, byrow = T)[,3]
repliscan.classTiming$timing    <- sub(';color', '', repliscan.classTiming$timing)
repliscan.classTiming$timing    <- sub('ES$', 'early', repliscan.classTiming$timing)
repliscan.classTiming$timing    <- sub('^LS', 'late', repliscan.classTiming$timing)
#repliscan.classTiming$timing    <- sub('ESLS', 'transition', repliscan.classTiming$timing)
repliscan.classTiming           <- repliscan.classTiming[repliscan.classTiming$chr %in% paste0('chr', c(1:22, 'X', 'Y')),]



### combine classTiming info from both approaches to get one final classification ###
# create combined table #
gr.l2r.classTiming       <- GRanges(seqnames = l2r.classTiming$chr, IRanges(start = l2r.classTiming$start, end = l2r.classTiming$stop))
gr.repliscan.classTiming <- GRanges(seqnames = repliscan.classTiming$chr, IRanges(start = repliscan.classTiming$start, end = repliscan.classTiming$stop))
overlap                  <- findOverlaps(gr.repliscan.classTiming, gr.l2r.classTiming)

combined.classTiming                  <- l2r.classTiming[,c(1,2,3,5,6)]
colnames(combined.classTiming)[4:5]   <- c('l2r.timing', "l2r.timing_withTransition")
combined.classTiming$repliscan.timing <- 'unknown'
combined.classTiming$repliscan.timing[subjectHits(overlap)] <- repliscan.classTiming$timing[queryHits(overlap)]
combined.classTiming$position.key <- paste(combined.classTiming$chr, combined.classTiming$start, combined.classTiming$stop, sep = ':')

# decide discordant regions without transition zones #
#1) if one algorithm classifies as unknown, then use classification of other algorithm
#2) repliscan classifies as ESLS and l2r as early or late, then use l2r classification
#4) classify non-matching regions (l2r early - repliscan late and l2r late - repliscan early) as discordant and treat them like unknowns
#  --> see denisty plot of values in regions that have been classified as early in one algorithm and late in the other one (and vice versa)
#  --> values in non-matching regions seem to be borderline regions and that is why there are discrepancies

combined.classTiming$combined.timing <- combined.classTiming$l2r.timing

ww <- which((combined.classTiming$l2r.timing == 'unknown') & (combined.classTiming$repliscan.timing %in% c('early', 'late'))) 
combined.classTiming$combined.timing[ww] <- combined.classTiming$repliscan.timing[ww]

ww <- which((combined.classTiming$l2r.timing == 'early' & combined.classTiming$repliscan.timing == 'late') | 
              (combined.classTiming$l2r.timing == 'late' & combined.classTiming$repliscan.timing == 'early'))
combined.classTiming$combined.timing[ww] <- 'discordant'


# decide discordant regions with transition zones #
#1) if one algorithm classifies as unknown, then use classification of other algorithm
#2) if repliscan classifies as transitions (ESLS) and l2r classifies as early or late, then use l2r classification
#  --> because repliscan classified as both, so also as early and late
#3) if l2r classifies as transition and repliscan as early or late, the use l2r classification
#  --> transition in l2r is classified from peak to valley or otherway round so their are definitely regions in between that are early or late
#4) classify non-matching regions (l2r early - repliscan late and l2r late - repliscan early) as discordant and treat them like unknowns
#  --> see denisty plot of values in regions that have been classified as early in one algorithm and late in the other one (and vice versa)
#  --> values in non-matching regions seem to be borderline regions and that is why there are discrepancies

combined.classTiming$combined.timing_withTransition <- combined.classTiming$l2r.timing_withTransition

ww <- which((combined.classTiming$l2r.timing_withTransition == 'unknown') & (combined.classTiming$repliscan.timing %in% c('early', 'late'))) 
combined.classTiming$combined.timing_withTransition[ww] <- combined.classTiming$repliscan.timing[ww]
ww <- which((combined.classTiming$repliscan.timing == 'unknown') & (combined.classTiming$l2r.timing_withTransition %in% c('early', 'late'))) 
combined.classTiming$combined.timing_withTransition[ww] <- combined.classTiming$l2r.timing_withTransition[ww]

ww <- which((combined.classTiming$l2r.timing_withTransition == 'early' & combined.classTiming$repliscan.timing == 'late') | 
            (combined.classTiming$l2r.timing_withTransition == 'late' & combined.classTiming$repliscan.timing == 'early'))
combined.classTiming$combined.timing_withTransition[ww] <- 'discordant'


#add smooth values to table
early.repliscan    <- as.data.frame(fread(paste0(REPLISCANDIR, 'ES_ratio.bedgraph'), col.names = c('chr', 'start', 'stop', 'value'), stringsAsFactors = F))
late.repliscan     <- as.data.frame(fread(paste0(REPLISCANDIR, 'LS_ratio.bedgraph'), col.names = c('chr', 'start', 'stop', 'value'), stringsAsFactors = F))
gr_early.repliscan <- makeGRangesFromDataFrame(early.repliscan)
gr_late.repliscan  <- makeGRangesFromDataFrame(late.repliscan)

gr_combined.classTiming <- makeGRangesFromDataFrame(combined.classTiming)
overlap.l2r             <- findOverlaps(gr_combined.classTiming, gr.l2r.classTiming, minoverlap = 5)
overlap.repliscan.early <- findOverlaps(gr_combined.classTiming, gr_early.repliscan, minoverlap = 5)
overlap.repliscan.late  <- findOverlaps(gr_combined.classTiming, gr_late.repliscan, minoverlap = 5)

combined.classTiming$l2r.value             <- NA
combined.classTiming$repliscan.early.value <- NA
combined.classTiming$repliscan.late.value  <- NA
combined.classTiming$l2r.value[queryHits(overlap.l2r)]                         <- l2r.classTiming$log2ratio[subjectHits(overlap.l2r)]
combined.classTiming$repliscan.early.value[queryHits(overlap.repliscan.early)] <- early.repliscan$value[subjectHits(overlap.repliscan.early)]
combined.classTiming$repliscan.late.value[queryHits(overlap.repliscan.late)]   <- late.repliscan$value[subjectHits(overlap.repliscan.late)]


# check density of values in non-matching regions #
ww <- which(combined.classTiming$combined.timing == 'discordant')
if(length(ww) != 0){
  discordant.regions      <- combined.classTiming[ww,]
  
  #l2r early and repliscan late
  pdf(paste0(PLOTDIR, PATIENT, '.density.discordantTiming.pdf'), width = 8, height = 6)
  if(sum(discordant.regions$l2r.timing == 'early' & discordant.regions$repliscan.timing == 'late') > 0){
    plot.data.l2r <- data.frame(key = paste(l2r.classTiming$chr[l2r.classTiming$timing == 'early'], l2r.classTiming$start[l2r.classTiming$timing == 'early'], l2r.classTiming$stop[l2r.classTiming$timing == 'early'], sep = ':'),
                                group = 'matching', 
                                value = l2r.classTiming$log2ratio[l2r.classTiming$timing == 'early'], 
                                stringsAsFactors = F)
    plot.data.l2r$group[plot.data.l2r$key %in% discordant.regions$position.key[discordant.regions$l2r.timing == 'early']] <- 'non-matching'
    plot.data.l2r$group <- factor(plot.data.l2r$group, levels = c('matching', 'non-matching'))
    
    p.l2r <- ggplot(plot.data.l2r) + 
      geom_density(aes(x = value, group = group, fill = group)) + 
      scale_fill_manual(name = 'Classification', values = c('#fdb863', '#8073ac'), labels = c('l2r_early:repliscan_early', 'l2r_early:repliscan_late')) +
      theme_bw() + 
      xlab('log2ratio') + ggtitle('Denisty of Log2ratio early and Repliscan late classified values')
    #p.l2r
    
    tmp      <- repliscan.classTiming[repliscan.classTiming$timing == 'late',]
    gr_tmp   <- makeGRangesFromDataFrame(tmp)
    overlap  <- findOverlaps(gr_early.repliscan, gr_tmp, minoverlap = 5)
    plot.data.repliscan.early <- early.repliscan[queryHits(overlap),]
    plot.data.repliscan.early$group <- 'matching'
    plot.data.repliscan.early$group[paste(plot.data.repliscan.early$chr, plot.data.repliscan.early$start, plot.data.repliscan.early$stop, sep = ':') %in% discordant.regions$position.key[discordant.regions$l2r.timing == 'early']] <- 'non-matching'
    plot.data.repliscan.early$group <- factor(plot.data.repliscan.early$group, levels = c('matching', 'non-matching'))
    
    p.repliscan.early <- ggplot(plot.data.repliscan.early) + 
      geom_density(aes(x = value, group = group, fill = group)) + 
      scale_fill_manual(name = 'Classification', values = c('#fdb863', '#8073ac'), labels = c('l2r_late:repliscan_late', 'l2r_early:repliscan_late')) +
      theme_bw() + 
      xlab('Repliscan early values') + ggtitle(" ")
    #p.repliscan.early
    
    tmp      <- repliscan.classTiming[repliscan.classTiming$timing == 'late',]
    gr_tmp   <- makeGRangesFromDataFrame(tmp)
    overlap  <- findOverlaps(gr_late.repliscan, gr_tmp, minoverlap = 5)
    plot.data.repliscan.late <- late.repliscan[queryHits(overlap),]
    plot.data.repliscan.late$group <- 'matching'
    plot.data.repliscan.late$group[paste(plot.data.repliscan.late$chr, plot.data.repliscan.late$start, plot.data.repliscan.late$stop, sep = ':') %in% discordant.regions$position.key[discordant.regions$l2r.timing == 'early']] <- 'non-matching'
    plot.data.repliscan.late$group <- factor(plot.data.repliscan.late$group, levels = c('matching', 'non-matching'))
    
    p.repliscan.late <- ggplot(plot.data.repliscan.late) + 
      geom_density(aes(x = value, group = group, fill = group)) + 
      scale_fill_manual(name = 'Classification', values = c('#fdb863', '#8073ac'), labels = c('l2r_late:repliscan_late', 'l2r_early:repliscan_late')) +
      theme_bw() + 
      xlab('Repliscan late values') + ggtitle(" ")
    #p.repliscan.late
    
    p.l2r             <- ggplotGrob(p.l2r + theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.5), "cm")))
    p.repliscan.early <- ggplotGrob(p.repliscan.early + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
    p.repliscan.late  <- ggplotGrob(p.repliscan.late + theme(plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm")))
    
    p.l2r$widths             <- unit.pmax(p.l2r$widths, p.repliscan.early$widths, p.repliscan.late$widths)
    p.repliscan.early$widths <- unit.pmax(p.l2r$widths, p.repliscan.early$widths, p.repliscan.late$widths)
    p.repliscan.late$widths  <- unit.pmax(p.l2r$widths, p.repliscan.early$widths, p.repliscan.late$widths)
    
    combined.list        <- list(p.l2r, p.repliscan.early, p.repliscan.late)
    combined.plot.layout <- matrix(c(rep(1, 1), rep(2,1), rep(3,1)), ncol = 1, byrow = T)
    combined.plot        <- arrangeGrob(grobs = combined.list, nrow = dim(combined.plot.layout)[1], ncol = 1, layout_matrix = combined.plot.layout)
    
    grid.draw(combined.plot)
  }
  
  
  
  #l2r late and repliscan early
  if(sum(discordant.regions$l2r.timing == 'late' & discordant.regions$repliscan.timing == 'early') > 0){
    plot.data.l2r <- data.frame(key = paste(l2r.classTiming$chr[l2r.classTiming$timing == 'late'], l2r.classTiming$start[l2r.classTiming$timing == 'late'], l2r.classTiming$stop[l2r.classTiming$timing == 'late'], sep = ':'),
                                group = 'matching', 
                                value = l2r.classTiming$log2ratio[l2r.classTiming$timing == 'late'], 
                                stringsAsFactors = F)
    plot.data.l2r$group[plot.data.l2r$key %in% discordant.regions$position.key[discordant.regions$l2r.timing == 'late']] <- 'non-matching'
    plot.data.l2r$group <- factor(plot.data.l2r$group, levels = c('matching', 'non-matching'))
    
    p.l2r <- ggplot(plot.data.l2r) + 
      geom_density(aes(x = value, group = group, fill = group)) + 
      scale_fill_manual(name = 'Classification', values = c('#fdb863', '#8073ac'), labels = c('l2r_late:repliscan_late', 'l2r_late:repliscan_early')) +
      theme_bw() + 
      xlab('log2ratio') + ggtitle('Denisty of Log2ratio late and Repliscan early classified values')
    #p.l2r
    
    tmp      <- repliscan.classTiming[repliscan.classTiming$timing == 'early',]
    gr_tmp   <- makeGRangesFromDataFrame(tmp)
    overlap  <- findOverlaps(gr_early.repliscan, gr_tmp, minoverlap = 5)
    plot.data.repliscan.early <- early.repliscan[queryHits(overlap),]
    plot.data.repliscan.early$group <- 'matching'
    plot.data.repliscan.early$group[paste(plot.data.repliscan.early$chr, plot.data.repliscan.early$start, plot.data.repliscan.early$stop, sep = ':') %in% discordant.regions$position.key[discordant.regions$l2r.timing == 'late']] <- 'non-matching'
    plot.data.repliscan.early$group <- factor(plot.data.repliscan.early$group, levels = c('matching', 'non-matching'))
    
    p.repliscan.early <- ggplot(plot.data.repliscan.early) + 
      geom_density(aes(x = value, group = group, fill = group)) + 
      scale_fill_manual(name = 'Classification', values = c('#fdb863', '#8073ac'), labels =  c('l2r_early:repliscan_early', 'l2r_late:repliscan_early')) +
      theme_bw() + 
      xlab('Repliscan early values') + ggtitle(" ")
    #p.repliscan.early
    
    tmp      <- repliscan.classTiming[repliscan.classTiming$timing == 'early',]
    gr_tmp   <- makeGRangesFromDataFrame(tmp)
    overlap  <- findOverlaps(gr_late.repliscan, gr_tmp, minoverlap = 5)
    plot.data.repliscan.late <- late.repliscan[queryHits(overlap),]
    plot.data.repliscan.late$group <- 'matching'
    plot.data.repliscan.late$group[paste(plot.data.repliscan.late$chr, plot.data.repliscan.late$start, plot.data.repliscan.late$stop, sep = ':') %in% discordant.regions$position.key[discordant.regions$l2r.timing == 'late']] <- 'non-matching'
    plot.data.repliscan.late$group <- factor(plot.data.repliscan.late$group, levels = c('matching', 'non-matching'))
    
    p.repliscan.late <- ggplot(plot.data.repliscan.late) + 
      geom_density(aes(x = value, group = group, fill = group)) + 
      scale_fill_manual(name = 'Classification', values = c('#fdb863', '#8073ac'), labels =  c('l2r_early:repliscan_early', 'l2r_late:repliscan_early')) +
      theme_bw() + 
      xlab('Repliscan late values') + ggtitle(" ")
    #p.repliscan.late
    
    p.l2r             <- ggplotGrob(p.l2r + theme(plot.margin = unit(c(1, 0.5, 0.1, 0.5), "cm")))
    p.repliscan.early <- ggplotGrob(p.repliscan.early + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
    p.repliscan.late  <- ggplotGrob(p.repliscan.late + theme(plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm")))
    
    p.l2r$widths             <- unit.pmax(p.l2r$widths, p.repliscan.early$widths, p.repliscan.late$widths)
    p.repliscan.early$widths <- unit.pmax(p.l2r$widths, p.repliscan.early$widths, p.repliscan.late$widths)
    p.repliscan.late$widths  <- unit.pmax(p.l2r$widths, p.repliscan.early$widths, p.repliscan.late$widths)
    
    combined.list        <- list(p.l2r, p.repliscan.early, p.repliscan.late)
    combined.plot.layout <- matrix(c(rep(1, 1), rep(2,1), rep(3,1)), ncol = 1, byrow = T)
    combined.plot        <- arrangeGrob(grobs = combined.list, nrow = dim(combined.plot.layout)[1], ncol = 1, layout_matrix = combined.plot.layout)
    
    grid.newpage()
    grid.draw(combined.plot)
  }
 dev.off()
}


# save combined timing table #
combined.classTiming <- combined.classTiming[,c('chr', 'start', 'stop', 'l2r.timing', 'l2r.timing_withTransition', 'repliscan.timing', 'l2r.value', 'repliscan.early.value', 'repliscan.late.value', 'combined.timing', 'combined.timing_withTransition')]
saveRDS(combined.classTiming, file = paste0(WORKDIR, PATIENT, '.replicationTiming.rds'))




### plot combined timing classififcation with log2ratio ###
output.dir <- paste0(PLOTDIR, '/combinedTiming/')
if(!dir.exists(output.dir)){
  dir.create(output.dir)
}

# without transition #
pdf(paste0(PLOTDIR, PATIENT, '.combinedTiming_log2ratio.pdf'), width = 22, height = 7)
for(chr in paste0('chr', c(1:22, 'X'))){
  sub             <- combined.classTiming[combined.classTiming$chr == chr,]
  sub$log2ratio   <- sub$l2r.value * 10
  sub$xmin        <- sub$start / 1000000
  sub$xmax        <- sub$stop / 1000000
  sub$ymin        <- floor(min(sub$log2ratio))
  sub$ymax        <- ceiling(max(sub$log2ratio))
  sub$combined.timing <- factor(sub$combined.timing, levels = c('early', 'late', 'discordant', 'unknown'))
  
  p <- ggplot() + 
    geom_rect(data = sub, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = combined.timing)) +
    scale_fill_manual(name = 'Timing', values = c('early'= '#de77ae', 'late' = '#7fbc41', 'discordant' = '#bababa', 'unknown' = '#f7f7f7'), labels = c('early', 'late', 'discordant', 'unknown')) + 
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = '#4d4d4d', size = 0.2) + 
    geom_path(data = sub, aes(x = xmin, y = log2ratio), size = 0.3) +
    theme_bw() + theme(legend.key = element_rect(colour = "black")) +
    ylab('smoothed Log2ratio') + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  # png(paste0(output.dir, PATIENT,'.', chr, '.combinedTiming_log2ratio.png'), width = 22, height = 7, units = 'cm', res = 300)
  plot(p)
  # dev.off()
}
dev.off()

# with transition zones #
pdf(paste0(PLOTDIR, PATIENT, '.combinedTiming_withTransition_log2ratio.pdf'), width = 22, height = 7)
for(chr in paste0('chr', c(1:22, 'X'))){
  sub             <- combined.classTiming[combined.classTiming$chr == chr,]
  sub$log2ratio   <- sub$l2r.value * 10
  sub$xmin        <- sub$start / 1000000
  sub$xmax        <- sub$stop / 1000000
  sub$ymin        <- floor(min(sub$log2ratio))
  sub$ymax        <- ceiling(max(sub$log2ratio))
  sub$combined.timing_withTransition <- factor(sub$combined.timing_withTransition, levels = c('early', 'late', 'transition', 'discordant', 'unknown'))
  
  p <- ggplot() + 
    geom_rect(data = sub, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = combined.timing_withTransition)) +
    scale_fill_manual(name = 'Timing', values = c('early'= '#de77ae', 'late' = '#7fbc41', 'transition' = '#fee0b6', 'discordant' = '#bababa', 'unknown' = '#f7f7f7'), labels = c('early', 'late', 'transition', 'discordant', 'unknown')) + 
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = '#4d4d4d', size = 0.2) + 
    geom_path(data = sub, aes(x = xmin, y = log2ratio), size = 0.3) +
    theme_bw() + theme(legend.key = element_rect(colour = "black")) +
    ylab('smoothed Log2ratio') + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  # png(paste0(output.dir, PATIENT,'.', chr, '.combinedTiming_withTransition_log2ratio.png'), width = 22, height = 7, units = 'cm', res = 300)
  plot(p)
  # dev.off()
}
dev.off()



