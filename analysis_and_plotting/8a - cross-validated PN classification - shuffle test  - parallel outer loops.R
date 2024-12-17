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
    n.cores.to.use = 7
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
# Downsample data
downsample.trial <- function(trial.data,
                             input.rate = 512,
                             output.rate = 256){
  output <- signal::resample(trial.data,
                             p = output.rate,
                             q = input.rate)
  return(output)
} # downsample.trial()


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


###
### Define function for classification
###

predict.time.series <- 
  function(.trained.classifiers, # trained.classifiers
           .train.folds, # train.folds
           test.data, # pn.data[[patient]][current.elecs]
           prediction.seeds = NA,
           iv.labels, # pn.info[[patient]]$word
           sample.labels.to.test, # pn.sample.labels
           words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
           shuffle.test.labels = FALSE,
           save.data.dir = save.classification.data.dir){  
    
    
    ### Inputs
    # .trained.classifiers: list of pre-trained classifiers, one for each train fold
    # .train.folds: training rows per cross-validation fold
    # test.data: list of dataframes (trial X elec), one per sample
    # prediction.seeds: vector of seed values, one per prediction run (so 1 for real data, 1000 for shuffled)
    # iv.labels: the Y values (true trial labels, i.e., what word they said at time 0)
    # sample.labels.to.test: run prediction at these samples
    # words: the unique Y values
    # shuffle.test.labels: real (FALSE) or shuffled (TRUE) loop
    # save.data.dir: root directory to save files
    
    
    # ### Predict!
    # ## Set up parallel processing
    # # Close any old parallel backends
    # unregister_dopar()
    # # Set up parallel workers in case caret wants to use it
    # cl <- makeCluster(n.cores.to.use, type = "FORK")
    # registerDoParallel(cl)
    
    # current.predictions <- 
    #   foreach(shuffle.loop = prediction.seeds) %dopar% {
    #     # shuffle.loop = prediction.seeds[1]
    
    for(shuffle.loop in prediction.seeds){
      # shuffle.loop = prediction.seeds[1]
      
      # Set seed
      set.seed(seed = shuffle.loop)
      
      current.predictions <- list()
      for(sample.loop in sample.labels.to.test){
        # sample.loop = sample.labels.to.test[1]
        
        # Loop thru RCVs
        current.predictions[[sample.loop]] <- list()
        for(rcv.loop in 1:length(.train.folds)){
          # rcv.loop = 1
          
          # Define test data
          current.test.data <- 
            cbind(data.frame('word' = iv.labels[-.train.folds[[rcv.loop]]]),
                  data.frame(bind_cols(lapply(test.data, function(x){
                    x[-.train.folds[[rcv.loop]], sample.loop]}))))
          
          # Shuffle?
          if(shuffle.test.labels){
            current.test.data$word <- as.factor(sample(as.character(current.test.data$word),
                                                       size = nrow(current.test.data), 
                                                       replace = FALSE))
          } # if(shuffle.test.labels){
          
          # Get prediction probabilities
          current.predictions[[sample.loop]][[rcv.loop]] <-
            predict(trained.classifiers[[rcv.loop]], 
                    newdata = current.test.data,
                    type = "prob")
          
          # Get winners
          current.predictions[[sample.loop]][[rcv.loop]]$winner <-
            predict(trained.classifiers[[rcv.loop]],
                    newdata = current.test.data,
                    type = "raw")
          current.predictions[[sample.loop]][[rcv.loop]]$actual <-
            current.test.data$word
          
        }; rm(rcv.loop)
      }; rm(sample.loop)
      
      
      ### Get summary stats
      # Trial accuracies
      current.predictions <- lapply(current.predictions, function(x){
        x <- bind_rows(x)
        x$accuracy <- with(x, as.numeric(winner == actual))
        return(x)
      })
      
      # Mean accuracy
      classifier.accuracies <-
        sapply(current.predictions, function(x){mean(x$accuracy)})
      
      # Prediction probabilities
      predictions <- list()
      for(sample.loop in names(current.predictions)){
        # sample.loop = names(current.predictions)[1]
        predictions[[sample.loop]] <- 
          reshape2::melt(current.predictions[[sample.loop]],
                         id.vars = 'actual',
                         measure.vars = words,
                         variable.name = 'model',
                         value.name = 'probability')
        predictions[[sample.loop]] <- 
          data.frame(with(predictions[[sample.loop]], 
                          tapply(probability, 
                                 list(actual, model), 
                                 mean)))
        predictions[[sample.loop]]$actual <- rownames(predictions[[sample.loop]])
        predictions[[sample.loop]]$sample.label <- sample.loop
        predictions[[sample.loop]] <- predictions[[sample.loop]][,c('sample.label','actual',words)]
        rownames(predictions[[sample.loop]]) <- NULL  
      }; rm(sample.loop)
      rm(current.predictions)
      
      predictions <- bind_rows(predictions)
      
      ### Save
      save.accuracies.dir <- 
        paste0(save.data.dir,
               'accuracies/',
               ifelse(shuffle.test.labels, 'shuffled/', 'real/'))
      save.predictions.dir <- 
        paste0(save.data.dir,
               'predictions/',
               ifelse(shuffle.test.labels, 'shuffled/', 'real/'))
      dir.create(save.accuracies.dir, showWarnings = FALSE, recursive = TRUE)
      dir.create(save.predictions.dir, showWarnings = FALSE, recursive = TRUE)
      
      # Save
      save(classifier.accuracies,
           file = paste0(save.accuracies.dir,'classifier.accuracies',
                         ifelse(is.na(shuffle.loop),'',paste0('_seed=', shuffle.loop)),
                         '.RData'))
      save(predictions,
           file = paste0(save.predictions.dir,'predictions',
                         ifelse(is.na(shuffle.loop),'',paste0('_seed=', shuffle.loop)),
                         '.RData'))
      
    }#; rm(shuffle.loop)
    
    #   } # foreach
    # 
    # ## Close parallel backend
    # stopCluster(cl)
    # unregister_dopar()
    
  } # predict.time.series()


### Load elec info
load(paste0(output.path, 'data/0a - definitions/elec info/elec info.RData'))


### Load ROIs and stages
load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/stages and rois/stages and rois.RData'))
roi.elecs <- lapply(roi.elecs, function(.elec){
  .elec[.elec %in% use.these.elecs]
})

### Load RTs
load(paste0(path, 'analysis/R/downsample data/warped data/high_gamma/median.rt.samples 256 Hz.RData'))
median.rt.times <- time.convert(median.rt.samples, "samples", "times", 256)
rm(median.rt.samples) # risky to keep since dealing with 256 and 512 Hz datasets


# ### Load Bayes factors for word information per electrode per time sample
# load(paste0(path,'analysis/R/track words during sentences/output/high_gamma/data/2a - bayesian anovas/anova Bayes Factors.RData')) # loads "word.aov.bfs"




###
### Classification
###

### Loop thru models
for(model.loop in rev(c('multinom','nnet'))){ ## UNCOMMENT
  # model.loop = c('multinom','nnet')[2] ## UNCOMMENT
  
  # Select some reasonable arbitrary values for model hyperparameters (not cross-validating for picture naming; not enough data)
  model.hps <- list(
    'multinom' = data.frame('decay' = .01),
    'nnet' = data.frame('size' = 4, 'decay' = .01)
  )[[model.loop]]
  
  ### Loop thru patients
  for(roi.loop in rois){ ## UNCOMMENT
    # roi.loop = rois[1] ## UNCOMMENT
    
    message("Beginning [",band.loop,", ",model.loop,", ",roi.loop,"] loop (with patients and training stages in parallel). ",Sys.time())
    
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
    
    # Just elecs in this ROI
    pn.data <- lapply(pn.data, function(x){
      x <- x[names(x)[names(x) %in% roi.elecs[[roi.loop]]]]
    })
    
    
    ### Predict!
    ## Set up parallel processing
    # Close any old parallel backends
    unregister_dopar()
    # Set up parallel workers in case caret wants to use it
    cl <- makeCluster(n.cores.to.use, type = "FORK")
    registerDoParallel(cl)
    
    # for(patient in rev(patients)){ ## UNCOMMENT
    #   # patient = patients[1] ## UNCOMMENT
    
    foreach(patient = rev(patients)) %:% 
      # patient = patients[1] ## UNCOMMENT
      
      foreach(stage.loop = names(stage.sample.labels)) %dopar% {
        # stage.loop = names(stage.sample.labels)[1] ## UNCOMMENT
        
        current.elecs <- roi.elecs[[roi.loop]][grepl(patient,roi.elecs[[roi.loop]])]
        
        # Only perform classification for patient rois with at least 3 elecs 
        if(length(current.elecs) >= 3){
          
          # for(stage.loop in names(stage.sample.labels)){ ## UNCOMMENT
          #   # stage.loop = names(stage.sample.labels)[1] ## UNCOMMENT
          
          
          current.sample.labels <- stage.sample.labels[[stage.loop]]
          
          ### Collapse across samples to get one training matrix per window per roi
          ## Progress update
          message('Beginning: ',band.loop,' - ',model.loop,' - ',patient,' - ',roi.loop,' - ',stage.loop,'. ',Sys.time())
          
          ## Define save directory
          save.classification.data.dir <- 
            paste0(output.path,
                   'data/8a - cross-validated PN classification - 50ms/',
                   model.loop,'/',
                   roi.loop,'/',
                   stage.loop,'/',
                   patient,
                   '/training samples averaged - unweighted/')
          
          
          ##
          ## Figure out if there's anything left to run for this loop
          ##
          
          ### Real data
          # Randomization seed
          real.seed <- 18
          # Accuracies and predictions filenames
          real.accs.filename <- 
            paste0(save.classification.data.dir,
                   'accuracies/real/classifier.accuracies_seed=',real.seed,'.RData')
          real.preds.filename <-
            paste0(save.classification.data.dir,
                   'predictions/real/predictions_seed=',real.seed,'.RData')
          
          ### Shuffled data
          
          ## Set up loops
          # Target number
          n.shuffle.loops <- 50
          
          # Which loops already complete?
          already.done.shuff.acc.seeds <- 
            gsub('classifier.accuracies_seed=','',
                 list.files(paste0(save.classification.data.dir,'accuracies/shuffled/')))
          already.done.shuff.acc.seeds <- 
            as.numeric(gsub('.RData','',already.done.shuff.acc.seeds))
          already.done.shuff.pred.seeds <- 
            gsub('predictions_seed=','',
                 list.files(paste0(save.classification.data.dir,'predictions/shuffled/')))
          already.done.shuff.pred.seeds <- 
            as.numeric(gsub('.RData','',already.done.shuff.pred.seeds))
          already.done.seeds <- 
            intersect(already.done.shuff.acc.seeds, already.done.shuff.pred.seeds)
          
          # Shuffles to do
          seeds.to.do <- ((1:n.shuffle.loops) * 100) + 18
          seeds.to.do <- seeds.to.do[! seeds.to.do %in% already.done.seeds]
          
          
          ### Loop if there's anything left to run
          if((! (file.exists(real.accs.filename) &
                 file.exists(real.preds.filename))) |
             (length(seeds.to.do) > 0)){
            
            ### Create train data
            # Get the PN data for this patient's roi-specific elecs this roi-specific stage (time window)
            current.raw.train.data <- lapply(pn.data[[patient]][current.elecs], 
                                             function(x){x[,current.sample.labels]})
            
            # Get weighted average of time samples per elec
            for(elec.loop in current.elecs){
              # elec.loop = current.elecs[1]
              
              # Weight training data
              current.raw.train.data[[elec.loop]] <- 
                apply(current.raw.train.data[[elec.loop]], 1, mean, na.rm = TRUE)
              
            }; rm(elec.loop)
            
            # Combine
            current.raw.train.data <- data.frame(bind_cols(current.raw.train.data))
            
            # Add IV to train data
            current.raw.train.data <- cbind(data.frame('word' = pn.info[[patient]]$word),
                                            current.raw.train.data)
            
            
            ##
            ## Train classifiers on real data
            ##
            
            ## Set up repeated k-fold cross-validation
            train.folds <- createMultiFolds(y = as.character(current.raw.train.data$word),
                                            k = 10,
                                            times = 2) # just one-fold validation is fine
            
            # Loop thru RCVs
            trained.classifiers <- list()
            for(rcv.loop in 1:length(train.folds)){
              # rcv.loop = 1
              
              # Define train data
              current.train.data <- current.raw.train.data[train.folds[[rcv.loop]],]
              
              # Begin suppress output
              sink("/dev/null") 
              
              # Train
              trained.classifiers[[rcv.loop]] <-
                train(word ~ .,
                      data = current.train.data,
                      method = model.loop,
                      preProcess = c("center","scale"),
                      tuneGrid = model.hps,
                      trControl = trainControl(method = "none",
                                               allowParallel = FALSE),
                      verbose = FALSE,
                      maxit = 1000)
              
              # End suppress output
              sink()
              
              # Remove unnecessary memory burden
              trained.classifiers[[rcv.loop]]$trainingData <- NA
              
            }; rm(rcv.loop)
            
            
            ##
            ## Real data
            ##
            
            # If not already run, run real data!
            if(! (file.exists(real.accs.filename) &
                  file.exists(real.preds.filename))){
              
              # Predict!
              predict.time.series(.trained.classifiers = trained.classifiers,
                                  .train.folds = train.folds,
                                  test.data = pn.data[[patient]][current.elecs],
                                  prediction.seeds = real.seed,
                                  iv.labels = pn.info[[patient]]$word,
                                  sample.labels.to.test = pn.sample.labels,
                                  words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                  shuffle.test.labels = FALSE,
                                  save.data.dir = save.classification.data.dir)
              
            } # if(! file.exists())
            
            
            ### Shuffled data
            ## If there are any shuffles to do
            if(length(seeds.to.do) > 0){
              
              predict.time.series(.trained.classifiers = trained.classifiers,
                                  .train.folds = train.folds,
                                  test.data = pn.data[[patient]][current.elecs],
                                  prediction.seeds = seeds.to.do,
                                  iv.labels = pn.info[[patient]]$word,
                                  sample.labels.to.test = pn.sample.labels,
                                  words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                  shuffle.test.labels = TRUE,
                                  save.data.dir = save.classification.data.dir)
              
            } # if(length(seeds.to.do) > 0){
            
          } # Loop if there's anything left to run
          
          # }; rm(weight.loop)
          # }; rm(stage.loop)
        } # if(length(current.elecs) > 3){
      } # foreach patients %:% stages
    
    # ## Close parallel backend
    stopCluster(cl)
    unregister_dopar()
    
    # }; rm(patient)
  }; rm(roi.loop)
  
}; rm(model.loop)

# }; rm(band.loop)












# Finish!
message('Script completed successfully. ',Sys.time())

