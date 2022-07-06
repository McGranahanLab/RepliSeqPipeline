##################################################################################################
#######                     Create a config file for each sample                           #######
##################################################################################################


### Command line arguments ###
cmdArgs        <- commandArgs(trailingOnly = TRUE)
designFile     <- cmdArgs[1]
releaseVersion <- cmdArgs[2]
work.dir       <- cmdArgs[3]
scratch.dir    <- cmdArgs[4]
script.dir     <- cmdArgs[5]
refGenome      <- cmdArgs[6]
chr.size.file  <- cmdArgs[7]
adapter        <- cmdArgs[8]
window.size    <- cmdArgs[9]
minRPKM        <- cmdArgs[10]
maxRPKM        <- cmdArgs[11]
span.size      <- cmdArgs[12]
qnorm.ref      <- cmdArgs[13]
threshold.l2   <- cmdArgs[14]
log.repliscan  <- cmdArgs[15]
nProcs         <- cmdArgs[16]
memRequest     <- cmdArgs[17]



### create config file for each sample ###
design      <- read.delim(designFile,header=T,sep="\t",stringsAsFactors=F,colClasses="character")

configFiles <- c()
for (sample in unique(design$Sample)) {
  patient    <- unlist(strsplit(sample, '_'))[1]
  
  #sample design
  sample.design <- design[design$Sample==sample,]
  
  #directories
  workDir        <- paste(work.dir, '/', releaseVersion,"/",sep='');
  scratchDir     <- paste(scratch.dir, '/' ,releaseVersion,"/",sep='');
  dataDir        <- "/camp/project/proj-tracerx-lung/tctProjects/mcgranahanLab/resources/farmFolder/tracerx/lung/data/";
  fastqDir       <- unlist(strsplit(sample.design$FQ1[1], '/'))
  fastqDir       <- paste(c(fastqDir[-length(fastqDir)], ''), collapse = '/')
  #fastqcDir      <- '/camp/lab/swantonc/inputs/babs-data-swanton/michelle.dietzen/robert.goldstone/DN18278/190328_K00102_0325_BH5KCGBBXY/fastqc/'
  
  #Create directory structure
  dirList <- list();
  dirList[["patientWorkDir"]]    <- patientWorkDir    <- paste(workDir, patient, "/",sep='');
  dirList[["patientScratchDir"]] <- patientScratchDir <- paste(scratchDir, patient ,"/",sep='');
  dirList[["cmdDir"]]            <- cmdDir            <- paste(patientWorkDir,"commandHistory/",sep='');
  dirList[["clippedFastqDir"]]   <- clippedFastqDir   <- paste(patientWorkDir,"ClippedFastq/",sep='');
  dirList[["fastqcDir"]]         <- fastqcDir         <- paste(patientWorkDir,"FastQC/",sep='');
  dirList[["bamFullDir"]]        <- bamFullDir        <- paste(patientWorkDir,"BAM/full/",sep='');
  dirList[["bamProcDir"]]        <- bamProcDir        <- paste(patientWorkDir,"BAM/processed/",sep='');
  dirList[["bamScratchDir"]]     <- bamScratchDir     <- paste(patientScratchDir,"BAM/",sep='');
  dirList[["countDir"]]          <- countDir          <- paste(patientWorkDir,"Counts/",sep='');
  dirList[["log2ratioDir"]]      <- log2ratioDir      <- paste(patientWorkDir,"Log2Ratio/",sep='');
  dirList[["plotDir"]]           <- plotDir           <- paste(patientWorkDir,"QC_Plots/",sep='');
  dirList[["repliScanDir"]]      <- repliScanDir      <- paste(patientWorkDir,"Repliscan/",sep='');
  
  dirCreated <- lapply(dirList,function(x) { if (!file.exists(x)) { dir.create(x, showWarnings=TRUE, recursive=TRUE, mode = "0775"); } });
  
  #fastq files
  nLanes <- nrow(sample.design)
  FQ1    <- sample.design$FQ1
  FQ2    <- sample.design$FQ2 #includes NAs for single-end reads
  
  #all types for patients
  allTypes <- unique(design$Type[design$Patient == patient])
  
  #create file
  configFile  <- paste0(cmdDir, sample, '.configFile.txt')
  
  sink(file=configFile);
  
  cat('PATIENT=', patient, '\n', sep='')
  cat('SAMPLE=', sample, '\n', sep='')
  cat('TYPE=', unique(sample.design$Type), '\n', sep='')
  cat('TYPES=( ', paste(allTypes, collapse = ' '), ' )\n', sep='')
  cat('RELEASE=', releaseVersion, '\n', sep='')
  cat('FQ1=( ', paste(FQ1, collapse = ' '), ' )\n', sep='')
  cat('FQ2=( ', paste(FQ2, collapse = ' '), ' )\n', sep='')
  cat('WORKDIR=', patientWorkDir, '\n', sep='')
  cat('SCRATCHDIR=', patientScratchDir, '\n', sep='')
  cat('CMDDIR=', cmdDir, '\n', sep='')
  cat('SCRIPTDIR=', script.dir, '\n', sep = '')
  cat('FASTQDIR=', fastqDir, '\n', sep='')
  cat('CLIPFASTQDIR=', clippedFastqDir, '\n', sep='')
  cat('FASTQCDIR=', fastqcDir, '\n', sep='')
  cat('BAMDIR=', patientWorkDir, 'BAM/', '\n', sep='')
  cat('COUNTDIR=', countDir, '\n', sep='')
  cat('LOG2RATIODIR=', log2ratioDir, '\n', sep='')
  cat('PLOTDIR=', plotDir, '\n', sep='')
  cat('REPLISCANDIR=', repliScanDir, '\n', sep='')
  cat('REFGENOME=', refGenome, '\n', sep = '')
  cat('CHROMSIZEFILE=', chr.size.file, '\n', sep = '')
  cat('ADAPTER=', adapter, '\n' ,sep ='')
  cat('WINDOWSIZE=', window.size, '\n', sep = '')
  cat('MINRPKM=', minRPKM, '\n', sep = '')
  cat('THRESHOLDL2R=', threshold.l2, '\n', sep = '')
  cat('MAXRPKM=', maxRPKM, '\n', sep = '')
  cat('SPANSIZE=', span.size, '\n', sep = '')
  cat('NORMQUANTILREF=', qnorm.ref, '\n', sep = '')
  cat('LOGREPLISCAN=', log.repliscan, '\n', sep = '')
  cat('NTHREATS=', nProcs, '\n', sep='')
  cat('MEM=', memRequest, '\n', sep='')
  
  sink()
  
  configFiles <- c(configFiles, configFile)

}

write.table(configFiles, file =  paste0(work.dir, '/', releaseVersion,"/configFile.list.txt"), sep = '\n', quote = F, row.names = F, col.names = F)




