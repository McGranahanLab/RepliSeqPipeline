#!/usr/bin/env Rscript

## load docopt package from CRAN
options(scipen=9999)
options(stringsAsFactors = F)
suppressMessages(library(vroom))
suppressMessages(library(preprocessCore))
suppressMessages(library(dplyr))


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

if(SPANSIZE<1){
  stop("lspan must indicate distance in bp")
} 

#### Function ####
# # function to normalise log2ratio by covergae
# norm_log2ratio <- function(bgfile, coverageFiles, outname, chr.to.use = paste0('chr', c(1:22, 'X', 'Y'))){
#                                               
#   l2r <- as.data.frame(vroom(bgfile, col_names = F))
#   l2r <- l2r[l2r[,1] %in% chr.to.use,]
#   
#   #calculate weights
#   early.cov <- as.data.frame(vroom(coverageFiles[1], col_names = F))
#   early.cov <- early.cov[early.cov[,1] %in% chr.to.use,]
#   late.cov  <- as.data.frame(vroom(coverageFiles[2], col_names = F))
#   late.cov  <- late.cov[late.cov[,1] %in% chr.to.use,]
#   
#   ww <- match(paste(early.cov[,1], early.cov[,2], early.cov[,3], sep = ':'), paste(late.cov[,1], late.cov[,2], late.cov[,3], sep = ':'))
#   weights.table <- data.frame(key = paste(early.cov[,1], early.cov[,2], early.cov[,3], sep = ':'), 
#                               early.cov = early.cov[,4],
#                               late.cov = late.cov[ww,4], stringsAsFactors = F)
# 
#   weights.table$weights <- rowSums(weights.table[,2:3])
#   weights.table$weights <- (weights.table$weights - min(weights.table$weights)) / (max(weights.table$weights) -  min( weights.table$weights))
#   
#   #normalise log2ratio values
#   ww           <- match(paste(l2r[,1], l2r[,2], l2r[,3], sep = ':'), weights.table$key)
#   norm.l2r     <- l2r
#   norm.l2r[,4] <- norm.l2r[,4] * weights.table$weights[ww]
#   
#   #save results
#   write.table(norm.l2r, file=outname, col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE)
#   return(norm.l2r)
# }


# function to loess smooth 
loess.smooth_log2ratio <- function(bgfile, lspan, outname){
    curbg     <- as.data.frame(vroom(bgfile, col_names = F))
    chroms    <- unique(curbg[,1])
		numchroms <- length(chroms)
		all       <- split(curbg,curbg[,1])
		
		lscores <- lapply(1:numchroms,function(i){
			cur        <- all[[i]]
		  curchrom   <- cur[1,1]
			chromlspan <- lspan/(max(cur[,3])-min(cur[,2]))
			
			cura <- as.data.frame(lapply(4:ncol(cur), function(k){
				cur[,k] <- tryCatch({
						predict(loess(cur[,k]~cur[,2],span=chromlspan),cur[,2])
					},warning = function(war){
						out <- predict(loess(cur[,k]~cur[,2],span=chromlspan),cur[,2])
						return(out)
					},error = function(err){
						return(cur[,k])
					}
				)
			}))
			
			cura           <- cbind(cur[,1:3],cura)
			colnames(cura) <- colnames(cur)
			return(cura)
		})
		smoothbg <- Reduce(rbind,lscores)
		rownames(smoothbg) <- NULL
		
		#make sure that the smoothed and input values are not too different
		colnames(curbg)    <- c('chr', 'start', 'stop', 'input')
		colnames(smoothbg) <- c('chr', 'start', 'stop', 'smooth')
		combined_df        <- curbg %>% left_join(smoothbg, by = c('chr', 'start', 'stop'))
		combined_df$abs_diff   <- abs(combined_df$input - combined_df$smooth)
		ww <- which(combined_df$abs_diff > 3)
		if(length(ww) > 0){
		  smoothbg <- combined_df[-1*ww, c('chr', 'start', 'stop', 'smooth')]
		} else {
		  smoothbg <- combined_df[, c('chr', 'start', 'stop', 'smooth')]
		}
		
		#save results
		write.table(smoothbg, file=outname, col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE)
		rm(curbg)
		gc()
}

############    Main     #############
# # normalise by coverage (RPKM) #
# TYPES <- unlist(strsplit(TYPES, ' '))
# TYPES <- grep('[(]', TYPES, invert = T, value = T)
# TYPES <- grep('[)]', TYPES, invert = T, value = T)
# 
# bgfile        <- paste0(LOG2RATIODIR, PATIENT, '.l2r.bedGraph')
# rpkmFiles     <- c(paste0(COUNTDIR, PATIENT, '_', TYPES[length(TYPES)-1], '.filteredRPKM.bedGraph'), paste0(COUNTDIR, PATIENT, '_',TYPES[length(TYPES)],'.filteredRPKM.bedGraph'))
# output.file   <- paste0(LOG2RATIODIR, PATIENT, '.norm_l2r.bedGraph')
# norm.l2r      <- norm_log2ratio(bgfile, rpkmFiles, outname = output.file, chr.to.use = paste0('chr', c(1:22, 'X')))


# quantile normalisation to make different samples comparable #
if(NORMQUANTILREF != 'NA'){
  target <- as.data.frame(vroom(NORMQUANTILREF, col_names = c('chr', 'start', 'stop', 'target')))
  input  <- as.data.frame(vroom(paste0(LOG2RATIODIR, PATIENT, '.l2r.bedGraph'), col_names = c('chr', 'start', 'stop', 'l2r')))
  combined_matrix <- input %>% 
    full_join(target, by = c('chr', 'start', 'stop')) %>%
    filter(chr %in% paste0('chr', c(1:22, 'X'))) %>%
    mutate(chr = factor(chr, levels = paste0('chr', c(1:22, 'X')))) %>%
    group_by(chr) %>%
    arrange(start, .by_group = T) %>%
    as.data.frame()
  combined_matrix$chr <- as.character(combined_matrix$chr)
  
  input_matrix  <- apply(combined_matrix[,4:5], 2, as.numeric)
  qnorm_results <- normalize.quantiles.use.target(input_matrix[,1,drop = F], input_matrix[,2])
  qnorm_results <- cbind(combined_matrix[,1:3],qnorm_results)
  
  qnorm.l2r <- qnorm_results[!is.na(qnorm_results[,4]),]
  colnames(qnorm.l2r) <- NULL
  
  #save results
  output.file <- paste0(LOG2RATIODIR, PATIENT, '.qnorm_l2r.bedGraph')
  write.table(qnorm.l2r, file=output.file, col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE)
  
  bgfile <- paste0(LOG2RATIODIR, PATIENT, '.qnorm_l2r.bedGraph')
} else {
  bgfile <- paste0(LOG2RATIODIR, PATIENT, '.l2r.bedGraph')
}


# loess smooth #
output.file <- paste0(LOG2RATIODIR, PATIENT, '.loess', SPANSIZE, '.bedGraph')
loess.smooth_log2ratio(bgfile, lspan = as.numeric(SPANSIZE), outname = output.file)







