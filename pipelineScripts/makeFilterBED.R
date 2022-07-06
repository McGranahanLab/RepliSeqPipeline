#!/usr/bin/env Rscript

## load docopt package from CRAN
options(scipen=9999)
suppressMessages(library(vroom))

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
TYPES <-  unlist(strsplit(TYPES, ' '))
TYPES <- TYPES[-c(1, length(TYPES))]
TYPES <- TYPES[TYPES != 'G1']
MAXRPKM <- as.numeric(MAXRPKM)
MINRPKM <- as.numeric(MINRPKM)


#read in RPKM values
combined.rpkm <- c()
for(type in TYPES){
  rpkm.data     <- as.data.frame(vroom(paste0(COUNTDIR, PATIENT,'_', type, '.w', WINDOWSIZE, '.RPKM.bedGraph'), col_names = c('chr', 'start', 'stop', 'value')))
  combined.rpkm <- cbind(combined.rpkm, rpkm.data$value)
}
combined.rpkm           <- cbind(combined.rpkm, rowSums(combined.rpkm))
colnames(combined.rpkm) <- c(TYPES, 'sum')
combined.rpkm           <- cbind(rpkm.data[,1:3], combined.rpkm)


#find bins with RPKM values > MAXRPKM per type and sum(RPKM) < 0.1 to exclude them from the bedfile
if(!is.na(MAXRPKM)){
  exclude.maxRPKM <- lapply(TYPES, function(x){
    which(combined.rpkm[,x] > MAXRPKM)
  })
  exclude.maxRPKM <- unique(unlist(exclude.maxRPKM))
} else {
  exclude.maxRPKM <- NULL
}

exclude.0RPKM <- lapply(TYPES, function(x){
  which(combined.rpkm[,x] == 0)
})
exclude.0RPKM <- unique(unlist(exclude.0RPKM))

exclude.minRPKM <- which(combined.rpkm[,'sum'] < MINRPKM)
exclude         <- unique(c(exclude.maxRPKM, exclude.minRPKM, exclude.0RPKM))

bedFile <- combined.rpkm[-exclude, 1:3]


#save bed file
write.table(bedFile, file = paste0(COUNTDIR, PATIENT, '.filtered.windows.bed'), col.names = F, row.names = F, quote = F, sep = '\t')




