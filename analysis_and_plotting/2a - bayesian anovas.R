### Do multiple regression on RSA values, cross-validated
### June 2023
### adam.milton.morgan@gmail.com

###
### Readme
###


###
### Setup
###

### Packages
library('beepr') # for playing beeps
library('dplyr') # for data organization, including bind_rows()
library('parallel')
library('foreach')
library('doParallel')
library('NMF')
library('BayesFactor') # for Bayes Factor
library('corrplot')
# library('lme4')

# Clean up
rm(list=ls())
cat("\014")
message("Begin multiple regressions on electrode RSA values (correlations across short time windows) on warped data. ",Sys.time())

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 4
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 5
    n.cores.to.use.nmf = 8
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 12
    n.cores.to.use.nmf = 18
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 10 # big RDMs
  n.cores.to.use.nmf = 32
}

### Lemmas
# Function to convert between samples, sample labels, and times
source(paste0(path,'/analysis/R/functions/time_convert.R'))
# Function for rolling average smoothing
source(paste0(path,'/analysis/R/functions/smoothing.R'))
# Function for plotting stacked trials in z (color) dimension
source(paste0(path,'/analysis/R/functions/stack_plot.R'))
# Close any old parallel backends
source(paste0(path,'/analysis/R/functions/unregister_dopar.R'))
# Subset to just good picture naming trials
source(paste0(path,'/analysis/R/functions/just_good_trial_functions.R'))
# Plot time series
source(paste0(path,'/analysis/R/functions/plot_time_series.R'))
# Add colored text to plots
source(paste0(path,'/analysis/R/functions/add_text_line_multiple_colors.R'))
# Create sigmoids
source(paste0(path,'/analysis/R/functions/sigmoid.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))
# Elementwise matrix apply
source(paste0(path,'/analysis/R/functions/elementwise_matrix_apply.R'))
# Squish Bayes Factors
source(paste0(path,'/analysis/R/functions/squish_bayes_factors.R'))
# Weighted mean and standard error functions
source(paste0(path,'/analysis/R/functions/weighted_summary_stats.R'))


### Set seed
set.seed(2024)


### Loop thru beta/high gamma data
# For cleanup
band.keep <- c(ls(), 'band.keep', 'band.loop')
# for(band.loop in c('high_gamma','beta')){ ## UNCOMMENT
band.loop = c('high_gamma','beta')[1] ## UNCOMMENT

# Clean up
rm(list = ls()[! ls() %in% band.keep])
gc()

# Save directory
output.path <- paste0(path, 'analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/',band.loop,'/')

### Load elec info
load(paste0(output.path, 'data/0a - definitions/elec info/elec info.RData'))

# Load data
message("Loading data... ",Sys.time())
load(paste0(path,
            'analysis/R/downsample data/warped data/',band.loop,'/warped ',band.loop,' data at 256 Hz.RData')) # loads "warped.data", "median.rt.samples", "patient.trial.info", and "sampling.rate"
message('...done! ', Sys.time())

### Clean up and get metadata
pn.data <- lapply(warped.data, function(x){x[['pn']]})
rm(warped.data)
pn.info <- lapply(patient.trial.info, function(x){x[['pn']]})
rm(patient.trial.info)

## Get rid of extremely early/late time samples (< 400ms before stim onset and > 500ms post articulation)
# Get "keep" sample ranges
pn.keep.sample.range <- c('start.time' = unname((-median.rt.samples['pn']) - time.convert(300, "times", "samples", 256)),
                          'end.time' = time.convert(1105, "times", "samples", 256))

# Subset
pn.data <- lapply(pn.data, function(x){
  lapply(x, function(y){
    current.samples <- time.convert(colnames(y), "sample.labels", "samples", 256)
    keep.cols <- which((current.samples >= pn.keep.sample.range[1]) & (current.samples <= pn.keep.sample.range[2]))
    return(y[,keep.cols])
  })
})

# Get new sample labels
pn.sample.labels <- colnames(pn.data[[1]][[1]])

# Good elecs:
patient.elecs <- lapply(pn.data, names)
patient.elecs <- lapply(patient.elecs, function(x){x[x %in% use.these.elecs]})
elecs.to.loop.thru <- data.frame('elec' = unlist(patient.elecs))
elecs.to.loop.thru$patient <- substr(elecs.to.loop.thru$elec, 1, 5)
rownames(elecs.to.loop.thru) <- elecs.to.loop.thru$elec
elecs.to.loop.thru$region_clinical <- elec.info[elecs.to.loop.thru$elec, 'region_clinical']

# Subset data to just use.these.elecs
for(patient in names(pn.data)){
  pn.data[[patient]] <- pn.data[[patient]][patient.elecs[[patient]]]
}; rm(patient)
gc()



###
### Bayesian anova
###

### Get Bayes factor for word information per electrode per time sample
## If already run, skip
word.aov.bf.path <- paste0(output.path,'data/2a - bayesian anovas/')
word.aov.bf.file <- 'anova Bayes Factors.RData'
if(file.exists(paste0(word.aov.bf.path, word.aov.bf.file))){
  load(paste0(word.aov.bf.path, word.aov.bf.file))
}else{
  
  ## Set up parallel processing
  # Close any old parallel backends
  unregister_dopar()
  # Set up parallel workers in case caret wants to use it
  cl <- makeCluster(n.cores.to.use + 2, type = "FORK")
  registerDoParallel(cl)
  
  word.aov.bfs <- 
    foreach(elec.loop = elecs.to.loop.thru$elec) %dopar% {
      # elec.loop = elecs.to.loop.thru$elec[1]
      
      # Set up
      patient <- elecs.to.loop.thru[elec.loop, 'patient']
      
      # Get Bayes Factors
      current.bfs <- sapply(pn.data[[patient]][[elec.loop]], function(x){
        current.data <- cbind('word' = pn.info[[patient]][,'word', drop = FALSE],
                              data.frame('ecog' = x))
        bf <- extractBF(anovaBF(ecog ~ word, data = current.data))$bf
        return(bf)
      })
      
      return(current.bfs)
    } # elec.loop
  ## Close parallel backend
  stopCluster(cl)
  unregister_dopar()
  
  message("Done w stats! ",Sys.time())
  beep()
  
  # Clean up
  names(word.aov.bfs) <- elecs.to.loop.thru$elec
  word.aov.bfs <- data.frame(bind_cols(word.aov.bfs))
  rownames(word.aov.bfs) <- colnames(pn.data[[1]][[1]])
  
  ### Save
  dir.create(word.aov.bf.path, showWarnings = FALSE, recursive = TRUE)
  save(word.aov.bfs,
       file = paste0(word.aov.bf.path, word.aov.bf.file))
  
} # if(file.exists(paste0(word.aov.bf.path, word.aov.bf.file))){}else{

# } # band.loop








# Finish!
message('Script completed successfully. ',Sys.time())







