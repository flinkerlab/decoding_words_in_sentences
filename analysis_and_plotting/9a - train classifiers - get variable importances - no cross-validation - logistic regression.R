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
library('ecp')
library('pheatmap')
library('caret')
library('nnet')


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
    n.cores.to.use = 15
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
# Classification
source(paste0(path,'/analysis/R/functions/classify_parallel_chunks.R'))


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

### Load data
# loads "median.rt.samples","patient.trial.info","sampling.rate","warped.data","warped.sample.labels"
load(paste0(path,'analysis/R/downsample data/warped data/high_gamma/warped high_gamma data at 256 Hz.RData')) 
patients <- names(warped.data)

## Clean up and get metadata
pn.data <- lapply(warped.data, function(x){x[['pn']]})
rm(warped.data)
pn.info <- lapply(patient.trial.info, function(x){x[['pn']]})
rm(patient.trial.info)

## Get rid of extremely early/late time samples (< 400ms before stim onset and > 500ms post articulation)
# Get "keep" sample ranges
pn.keep.sample.range <- c('start.time' = unname((-median.rt.samples['pn']) - time.convert(205, "times", "samples", 256)),
                          'end.time' = time.convert(605, "times", "samples", 256))

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


### Load AOV BFs for weighting PN data in windows
load(paste0(path,'analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/data/2a - bayesian anovas/anova Bayes Factors.RData'))


### Load region assignments (roi.elecs etc.)
load(paste0(output.path, 'data/0a - definitions/roi metadata/ROIs and temporal lobe splits.RData'))
roi.elecs <- roi.elecs[c('IFG','MFG','SMC','ATL','MTL','PTL','IPL')]
rois <- names(roi.elecs)


### Define windows - every 50ms
breakpoint.times <- rev(seq(250, time.convert(-median.rt.samples['pn'], "samples", "times", 256), by = -50))
breakpoint.samples <- time.convert(breakpoint.times, "times", "samples", 256)

# Define ranges
stage.ranges <- data.frame(
  'start.time' = breakpoint.times[-length(breakpoint.times)],
  'end.time' = breakpoint.times[-1],
  'start.sample' = breakpoint.samples[-length(breakpoint.samples)],
  'end.sample' = breakpoint.samples[-1] - 1)
stage.ranges$label <- 
  with(stage.ranges, 
       paste0('window_',
              gsub('-','neg',as.character(start.time)),
              '_to_',
              gsub('-','neg',as.character(end.time)),
              'ms'))
rownames(stage.ranges) <- stage.ranges$label

# Get sample labels per stage
stage.sample.labels <- list()
for(stage.loop in rownames(stage.ranges)){
  stage.sample.labels[[stage.loop]] <-
    time.convert(stage.ranges[stage.loop,'start.sample']:stage.ranges[stage.loop,'end.sample'], 
               "samples", "sample.labels", 256)
}; rm(stage.loop)


### Save stages and rois
save.stages.and.rois <- c('roi.elecs',
                          'rois',
                          'stage.ranges',
                          'stage.sample.labels')
save.stages.and.rois.dir <- 
  paste0(output.path,
         'data/9a - train classifiers - get variable importances - 50ms/stages and rois/')
dir.create(save.stages.and.rois.dir, showWarnings = FALSE, recursive = TRUE)
save(list = save.stages.and.rois,
     file = paste0(save.stages.and.rois.dir, 'stages and rois.RData'))


###
### Classification
###

### Define hyperparameters to try
grids <- list(
  "multinom" = expand.grid(decay = c(10^c(-4:0))),
  "nnet" = expand.grid(size = 3:6,
                       decay = c(10^c(-4:0))))
model.types <- names(grids)

### Loop thru models
for(model.loop in model.types){
  # model.loop = model.types[2]
  
  ### Loop thru patients
  trained.classifiers <- list()
  winning.hps <- list()
  variable.importances <- list()
  for(cluster.loop in rois){
    # cluster.loop = rois[1]
    
    ### Loop thru stages
    trained.classifiers[[cluster.loop]] <- list()
    winning.hps[[cluster.loop]] <- list()
    variable.importances[[cluster.loop]] <- list()
    for(stage.loop in names(stage.sample.labels)){
      # stage.loop = names(stage.sample.labels)[1]
      
      # For cleaning up
      patient.keep <- c(ls(), 'patient.keep', 'patient')
      
      trained.classifiers[[cluster.loop]][[stage.loop]] <- list()
      winning.hps[[cluster.loop]][[stage.loop]] <- list()
      variable.importances[[cluster.loop]][[stage.loop]] <- list()
      for(patient in patients){
        # patient = patients[1]
        
        # Clean up
        rm(list = ls()[! ls() %in% patient.keep])
        gc()
        
        current.elecs <- roi.elecs[[cluster.loop]][grepl(patient,roi.elecs[[cluster.loop]])]
        
        # Just those in word.aov.bfs (excludes visual and bad elecs etc.)
        current.elecs <- current.elecs[current.elecs %in% colnames(word.aov.bfs)]
        
        # Only perform classification for patient rois with at least 3 elecs 
        if(length(current.elecs) >= 3){
          
          current.sample.labels <- stage.sample.labels[[stage.loop]]
          
          # Weights for weighted average
          current.weights <-  word.aov.bfs[current.sample.labels, current.elecs]
          
          # Get the PN data for this patient's cluster-specific elecs this cluster-specific stage (time window)
          train.data <- lapply(pn.data[[patient]][current.elecs], 
                               function(x){x[,current.sample.labels]})
          
          # Get weighted average of time samples per elec
          for(elec.loop in current.elecs){
            # elec.loop = current.elecs[1]
            
            # Weight training data
            train.data[[elec.loop]] <- 
              apply(train.data[[elec.loop]], 1, function(x){
                weighted.average(x, w = current.weights[,elec.loop])})
            
          }; rm(elec.loop)
          
          # Combine
          train.data <- data.frame(bind_cols(train.data))
          
          # Add IV to train data
          train.data <- cbind(data.frame('word' = as.factor(pn.info[[patient]]$word)),
                              train.data)
          
          # Begin suppress output
          sink("/dev/null")
          
          # Tune parameters to find best fit
          tuning.classifier <-
            train(word ~ .,
                  data = train.data,
                  method = model.loop,
                  preProcess = c("center","scale"),
                  tuneGrid = grids[[model.loop]],
                  trControl = trainControl(method = "repeatedcv",
                                           number = 10,
                                           repeats = 3,
                                           allowParallel = TRUE),
                  verbose = FALSE,
                  maxit = 1000)
          
          # End suppress output
          sink()
          
          # Get best hyperparameters
          winning.hps[[cluster.loop]][[stage.loop]][[patient]] <- 
            tuning.classifier$bestTune
          
          # Train model with the winning hyperparmeter
          trained.classifiers[[cluster.loop]][[stage.loop]][[patient]] <- 
            train(word ~ .,
                  data = train.data,
                  method = model.loop,
                  preProcess = c("center","scale"),
                  tuneGrid = winning.hps[[cluster.loop]][[stage.loop]][[patient]],
                  trControl = trainControl(method = "none",
                                           allowParallel = FALSE),
                  verbose = FALSE,
                  maxit = 1000)
          
          # Reduce storage
          trained.classifiers[[cluster.loop]][[stage.loop]][[patient]]$trainingData <- NA
          
          # Get variable importances
          variable.importances[[cluster.loop]][[stage.loop]][[patient]] <- 
            data.frame(t(varImp(trained.classifiers[[cluster.loop]][[stage.loop]][[patient]], scale = FALSE)$importance))
          
          message('Completed: ',model.loop,' - ',band.loop,' - ',patient,' - ',cluster.loop,' - ',stage.loop,'. ',Sys.time())
          
        } # if(length(current.elecs) > 3){
      }; rm(patient)
    }; rm(stage.loop)
  }; rm(cluster.loop)
  
  ### Save directories
  save.classifier.dir <- 
    paste0(output.path,
           'data/9a - train classifiers - get variable importances - 50ms/',
           '/weighted by aov.bfs/',
           model.loop,'/')
  dir.create(save.classifier.dir, showWarnings = FALSE, recursive = TRUE)
  save(list = c('trained.classifiers','winning.hps'),
       file = paste0(save.classifier.dir, 'trained.classifiers.RData'))
  save(variable.importances,
       file = paste0(save.classifier.dir, 'variable.importances.RData'))
  
}; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


