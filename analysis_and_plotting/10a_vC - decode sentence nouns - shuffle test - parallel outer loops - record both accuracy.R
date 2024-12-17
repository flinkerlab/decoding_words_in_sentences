### Decode sentence nouns
### August 2024
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
    n.cores.to.use = 10 # for 10a_vB: max 11
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
# Reverse hierarchy of list of dataframes
source(paste0(path,'/analysis/R/functions/reverse.data.frame.list.hierarchy.R'))
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

### Load elec info
load(paste0(output.path, 'data/0a - definitions/elec info/elec info.RData'))


### Load ROIs and stages
load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/stages and rois/stages and rois.RData'))


### Load RTs
load(paste0(path, 'analysis/R/downsample data/warped data/high_gamma/median.rt.samples 256 Hz.RData'))
median.rt.times <- time.convert(median.rt.samples, "samples", "times", 256)
rm(median.rt.samples) # risky to keep since dealing with 256 and 512 Hz datasets


# Words
words <- c('chicken','dog','dracula','frankenstein','ninja','nurse')


###
### Classification
###

### Loop thru models
# For clean up
model.keep <- c(ls(), 'model.keep', 'model.loop')

model.types <- list.files(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/'))
for(model.loop in model.types[1]){ ## UNCOMMENT
  # model.loop = rev(model.types)[1] ## UNCOMMENT
  
  # Clean up
  rm(list = ls()[! ls() %in% model.keep])
  gc()
  
  ### Load classifiers
  load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/',
              model.loop,
              '/trained.classifiers.RData'))
  
  ### Patients
  patients <- unique(sapply(strsplit(elec.info$patient_elec, split = '_', fixed = TRUE), function(x){x[1]}))
  
  
  ### Decode!
  ## Set up parallel processing
  # Close any old parallel backends
  unregister_dopar()
  # Set up parallel workers in case caret wants to use it
  cl <- makeCluster(n.cores.to.use, type = "FORK")
  registerDoParallel(cl)
  
  foreach(roi.loop = rois) %:% 
    # roi.loop = rois[2]
    
    foreach(patient = patients) %dopar% {
      # patient = patients[1]
      
      # Only run if this patient has >= 3 elecs in this ROI (i.e., a classifier was trained)
      current.cluster.patients <- unique(unname(unlist(lapply(trained.classifiers[[roi.loop]], names))))
      if(patient %in% current.cluster.patients){
        
        # Status update
        message('Beginning prediction loop: ',band.loop,', ',model.loop,', ',roi.loop,', ',patient,'. ',Sys.time())
        start.time <- Sys.time()
        
        # This patient's elecs in this cluster
        current.elecs <- roi.elecs[[roi.loop]][grepl(patient,roi.elecs[[roi.loop]])]
        
        ### Load patient sample data
        load(paste0(path,
                    'data/',
                    patient, 
                    '/data/epoched/elec_data_without_bad_trials/multiband/',
                    'locked_to_production_onset',
                    '/z_scored/', 
                    patient,' elec data.RData')) # reads in list "all.data"
        
        # Trial data
        trial.info <- read.csv(paste0(path,
                                      '/data/',
                                      patient,
                                      '/data/epoched/',
                                      patient,
                                      '_trial_labels_without_bad_trials.csv'))
        
        # Subset to just sentence and list production data
        data <- lapply(all.data, function(x){
          x <- rbind(just.sentence.noun.verb.data(x),
                     just.list.noun.data(x))
          # Remove linguistic metadata columns from data (identical to trial.info)
          x <- x[,-which(colnames(x) %in% colnames(trial.info))]
          return(x)
        })
        trial.info <- rbind(just.sentence.noun.verb.data(trial.info),
                            just.list.noun.data(trial.info))
        
        # Free up some memory
        rm(all.data)
        gc()
        
        ### Downsample data
        data <- lapply(data, function(z){
          output <- data.frame(t(apply(z, 1, downsample.trial)))
          # Get rid of times outside original window
          output <- output[, colnames(output) %in% colnames(z)]
          # Convert colnames to new sampling rate
          colnames(output) <- time.convert(colnames(output), 
                                           input.type = "sample.labels", 
                                           output.type = "sample.labels", 
                                           input.sampling.rate = 512, 
                                           output.sampling.rate = 256)
          return(output)
        })
        
        ### Limit to just samples of interest
        decode.these.samples.256 <- 
          time.convert(-median.rt.times['sp'] - 100, "times", "samples", 256):
          time.convert(600, "times", "samples", 256)
        decdode.these.sample.labels.256 <-
          time.convert(decode.these.samples.256, "samples", "sample.labels", 256)
        data <- lapply(data, function(x){
          x[,decdode.these.sample.labels.256]
        })
        
        # Convert from list of elecs to list of samples
        data <- reverse.data.frame.list.hierarchy(data)
        
        
        ### Loop thru stages
        # For clean up
        stage.keep <- c(ls(), 'stage.keep', 'stage.loop')
        
        # Loop!
        for(stage.loop in names(stage.sample.labels)){ ## UNCOMMENT
          # stage.loop = names(stage.sample.labels)[1] ## UNCOMMENT
          
          # Clean up
          rm(list = ls()[! ls() %in% stage.keep])
          gc()
          
          
          ### Get row indices for trials to average for each representation of interest
          # Add test class column
          trial.data <- trial.info
          trial.data$test.class <- NA
          trial.data$test.class[which(trial.data$case == 's' & trial.data$pos == 'n' & trial.data$verb_voice == 'active')] <- 'subject.noun.active'
          trial.data$test.class[which(trial.data$case == 's' & trial.data$pos == 'n' & trial.data$verb_voice == 'passive')] <- 'subject.noun.passive'
          trial.data$test.class[which(trial.data$case == 's' & trial.data$pos == 'd' & trial.data$verb_voice == 'active')] <- 'subject.det.active'
          trial.data$test.class[which(trial.data$case == 's' & trial.data$pos == 'd' & trial.data$verb_voice == 'passive')] <- 'subject.det.passive'
          trial.data$test.class[which(trial.data$case == 'do' & trial.data$pos == 'n')] <- 'object.noun.active'
          trial.data$test.class[which(trial.data$case == 'bo' & trial.data$pos == 'n')] <- 'object.noun.passive'
          trial.data$test.class[which(trial.data$case == 'do' & trial.data$pos == 'd')] <- 'object.det.active'
          trial.data$test.class[which(trial.data$case == 'bo' & trial.data$pos == 'd')] <- 'object.det.passive'
          trial.data$test.class[which(trial.data$pos == 'v' & trial.data$verb_voice == 'active')] <- 'verb.active'
          trial.data$test.class[which(trial.data$pos == 'v' & trial.data$verb_voice == 'passive')] <- 'verb.passive'
          trial.data$test.class[which(trial.data$pos == 'a' & trial.data$verb_voice == 'active')] <- 'aux.active'
          trial.data$test.class[which(trial.data$pos == 'a' & trial.data$verb_voice == 'passive')] <- 'aux.passive'
          trial.data$test.class[which(trial.data$pos == 'p')] <- 'passive.be'
          trial.data$test.class[which(trial.data$pos == 'b')] <- 'passive.by'
          trial.data$test.class[which(trial.data$case == '1')] <- 'list1'
          trial.data$test.class[which(trial.data$case == '2')] <- 'list2'
          
          # Factorize
          trial.data$test.class <- as.factor(trial.data$test.class)
          
          ## Groups of trials you'll want to look at later
          lx.classes <- list()
          all.test.classes <- c('subject.noun.active',
                                'verb.active',
                                'object.noun.active',
                                'subject.noun.passive',
                                'verb.passive',
                                'object.noun.passive',
                                'list1',
                                'list2')
          
          for(i in all.test.classes){
            lx.classes[[i]] <- i
          }; rm(i)
          lx.classes[['subject.noun']] <- c('subject.noun.active','subject.noun.passive')
          # lx.classes[['subject.det']] <- c('subject.det.active','subject.det.passive')
          lx.classes[['object.noun']] <- c('object.noun.active','object.noun.passive')
          # lx.classes[['object.det']] <- c('object.det.active','object.det.passive')
          lx.classes[['sentence.noun']] <- c('subject.noun.active','subject.noun.passive','object.noun.active','object.noun.passive')
          # lx.classes[['det']] <- c('subject.det.active','subject.det.passive','object.det.active','object.det.passive')
          lx.classes[['verb']] <- c('verb.active','verb.passive')
          # lx.classes[['aux']] <- c('aux.active','aux.passive')
          # if(any(grepl("list", trial.data$test.class))){
          lx.classes[['list.noun']] <- c('list1','list2')
          # }
          
          # Corresponding rows in test data
          test.trial.rows.by.lx.class <- list()
          for(word.loop in c('noun1','noun2')){
            # word.loop = c('noun1','noun2')[1]
            
            test.trial.rows.by.lx.class[[word.loop]] <- list()
            for(lx.class.loop in names(lx.classes)){ 
              # (lx.class.loop = names(lx.classes)[3])
              
              test.trial.rows.by.lx.class[[word.loop]][[lx.class.loop]] <-
                which(trial.data[! is.na(trial.data[,word.loop]),]$test.class %in% lx.classes[[lx.class.loop]])
            }; rm(lx.class.loop)
          }; rm(word.loop)
          
          # Add 'both' rows
          test.trial.rows.by.lx.class[['both']] <- list()
          for(lx.class.loop in names(lx.classes)){ 
            # (lx.class.loop = names(lx.classes)[3])
            
            test.trial.rows.by.lx.class[['both']][[lx.class.loop]] <-
              which(trial.data[(! is.na(trial.data[,'noun1'])) & (! is.na(trial.data[,'noun2'])),]$test.class %in% lx.classes[[lx.class.loop]])
          }; rm(lx.class.loop)
          
          # Clean up
          # Open up more memory
          rm(trial.data, all.test.classes)
          
          
          ##
          ## Classification function (use exact same function for real and shuffle loops)
          ##
          
          # Define prediction function
          use.cluster.classifiers.to.predict <- 
            function(all.test.data,
                     test.trial.info,
                     current.patient,
                     trained.classifier,
                     shuffle.test = FALSE,
                     .seed = NA,
                     .save = TRUE,
                     .save.path = NA,
                     .save.name = NA,
                     .output = TRUE){
              
              if(.save){if(is.na(.save.path) | is.na(.save.name)){
                message("Error: .save is TRUE, but .save.path and/or .save.name are not specified.")
                break
              }}
              
              # Set random seed
              set.seed(seed = .seed)
              
              # Shuffle test labels if shuffle.loop
              if(shuffle.test){
                test.trial.info[,'noun1'] <- sample(test.trial.info[,'noun1'])
                # Shuffle independently? Can't decide... for now go for it
                test.trial.info[,'noun2'] <- sample(test.trial.info[,'noun2'])
              } # if(shuffle.test)
              
              # Storage template
              accuracies.by.lx.class.template <- data.frame('sample.label' = names(all.test.data))
              
              ### Loop thru first and second noun
              accuracies.by.lx.class <- list()
              for(word.loop in c('noun1','noun2')){
                # word.loop = c('noun1','noun2')[1]
                
                current.results <- 
                  data.frame(matrix(nrow = nrow(test.trial.info[! is.na(test.trial.info[,word.loop]),]),
                                    ncol = length(all.test.data)))
                colnames(current.results) <- names(all.test.data)
                
                ### Loop thru samples
                for(sample.loop in 1:length(all.test.data)){
                  # sample.loop = 1
                  
                  current.sample <- names(all.test.data)[sample.loop]
                  current.data <- cbind(data.frame('word' = test.trial.info[,word.loop]),
                                        all.test.data[[current.sample]][, current.elecs])
                  # Get rid of NAs
                  current.data <- droplevels(subset(current.data, ! is.na(word)))
                  
                  # Predict: Classify
                  current.predictions <- predict(trained.classifier, 
                                                 newdata = current.data)
                  current.results[,current.sample] <- 
                    as.numeric(as.character(current.data$word) == as.character(current.predictions))
                  
                  # Clean up
                  rm(current.predictions, current.sample, current.data)
                  
                }; rm(sample.loop)
                
                # Get mean accuracies for each linguistic class
                accuracies.by.lx.class[[word.loop]] <- accuracies.by.lx.class.template
                for(lx.loop in names(lx.classes)){ # (lx.loop = names(lx.classes)[1])
                  accuracies.by.lx.class[[word.loop]][,lx.loop] <- 
                    colMeans(current.results[test.trial.rows.by.lx.class[[word.loop]][[lx.loop]],])
                }; rm(lx.loop)
                rm(current.results)
                
              }; rm(word.loop)
              
              ### Now do the same but counting either noun1 or noun2 as correct
              current.results <- 
                data.frame(matrix(nrow = nrow(test.trial.info[(! is.na(test.trial.info$noun1)) & (! is.na(test.trial.info$noun2)),]),
                                  ncol = length(all.test.data)))
              colnames(current.results) <- names(all.test.data)
              
              ## Loop thru samples
              for(sample.loop in 1:length(all.test.data)){
                # sample.loop = 1
                
                current.sample <- names(all.test.data)[sample.loop]
                current.labels <- data.frame('noun1' = test.trial.info[,'noun1'],
                                             'noun2' = test.trial.info[,'noun2'])
                current.data <- all.test.data[[current.sample]][, current.elecs]
                # Get rid of NAs
                current.data <- current.data[(! is.na(current.labels$noun1)) & (! is.na(current.labels$noun2)),]
                current.labels <- current.labels[(! is.na(current.labels$noun1)) & (! is.na(current.labels$noun2)),]
                # Predict
                current.predictions <- predict(trained.classifier, 
                                               newdata = current.data)
                for(trial.loop in 1:nrow(current.data)){
                  current.results[trial.loop, current.sample] <- 
                    as.numeric(as.character(current.predictions[trial.loop]) %in% 
                                 as.character(unlist(current.labels[trial.loop,c('noun1','noun2')])))
                }; rm(trial.loop)
                rm(current.predictions, current.sample, current.data, current.labels)
                
              }; rm(sample.loop)
              
              # Get mean accuracies for each linguistic class
              accuracies.by.lx.class[['both']] <- accuracies.by.lx.class.template
              for(lx.loop in names(lx.classes)){ # (lx.loop = names(lx.classes)[1])
                accuracies.by.lx.class[['both']][,lx.loop] <- 
                  colMeans(current.results[test.trial.rows.by.lx.class[['both']][[lx.loop]],])
              }; rm(lx.loop)
              rm(current.results)
              
              
              ### Now get prediction probabilities
              probabilities.by.lx.class <- list()
                
              # Get data rows where neither noun1 nor noun2 are NAs
              keep.rows <- which((! is.na(test.trial.info$noun1)) & (! is.na(test.trial.info$noun2)))
              probs.template <- data.frame(matrix(nrow = length(keep.rows),
                                                  ncol = length(words),
                                                  dimnames = list(NULL,
                                                                  c('noun1','noun2','other1','other2','other3','other4'))))
              
                ### Loop thru samples
                for(sample.loop in 1:length(all.test.data)){
                  # sample.loop = 1
                  
                  current.sample <- names(all.test.data)[sample.loop]
                  current.data <- cbind(data.frame('word' = NA),
                                        all.test.data[[current.sample]][keep.rows, current.elecs])
                  
                  # Predict: Classify
                  current.probs <- predict(trained.classifier, 
                                                 newdata = current.data,
                                                type = "prob")
                  temp <- probs.template
                  for(row in 1:nrow(temp)){
                    # row = 1
                    
                    # Get all nouns and label them
                    current.noun1 <- test.trial.info[keep.rows[row],'noun1']
                    current.noun2 <- test.trial.info[keep.rows[row],'noun2']
                    current.other.words <- words[! words %in% c(current.noun1,current.noun2)]
                    current.noun.labels <- 
                      c('noun1' = current.noun1,
                      'noun2' = current.noun2,
                      'other1' = current.other.words[1],
                      'other2' = current.other.words[2],
                      'other3' = current.other.words[3],
                      'other4' = current.other.words[4])
                    
                    for(col in colnames(probs.template)){
                      # col = colnames(probs.template)[1]
                      temp[row, col] <- current.probs[row, current.noun.labels[col]]
                    }; rm(col)
                    
                    # Clean up
                    rm(current.noun1, current.noun2, current.other.words, current.noun.labels)
                    
                  }; rm(row)
                  
                  # Get mean probs for each linguistic class
                  probabilities.by.lx.class <- accuracies.by.lx.class.template
                  for(lx.loop in names(lx.classes)){ # (lx.loop = names(lx.classes)[1])
                    accuracies.by.lx.class[[word.loop]][,lx.loop] <- 
                      colMeans(current.results[test.trial.rows.by.lx.class[[word.loop]][[lx.loop]],])
                  }; rm(lx.loop)
                  
                  
                  # Clean up
                  rm(current.sample, current.data, temp)
                  
                }; rm(sample.loop)
                
                
                
              
              ### Save
              if(.save){
                
                # Save classification accuracies 
                dir.create(paste0(.save.path, 'accuracies/'), showWarnings = FALSE, recursive = TRUE)
                save(accuracies.by.lx.class,
                     file = paste0(.save.path, 'accuracies/', .save.name))
                
                # Save prediction probabilities
                dir.create(paste0(.save.path, 'prediction probabilities/'), showWarnings = FALSE, recursive = TRUE)
                save(probabilities.by.lx.class,
                     file = paste0(.save.path, 'prediction probabilities/', .save.name))
                
              } # if(save){
              
              if(.output){
                return(accuracies.by.lx.class)
              } # if(.output){
              
            } # use.cluster.classifiers.to.predict()
          
          
          
          ###
          ### Predict
          ### 
          
          ##
          ## Real data
          ## 
          
          save.path.real.data <- 
            paste0(output.path,
                   '/data/10a_vC - decode sentence nouns - shuffle test - 50ms/',
                   model.loop,'/',
                   'real data/',
                   roi.loop,'/',
                   stage.loop,'/',
                   patient,'/')
          save.name.real.data <- paste0(roi.loop,' - ',stage.loop,' - ',patient,'.RData')
          if(! file.exists(paste0(save.path.real.data, save.name.real.data))){
            
            # Time
            message("Beginning classification of real data for ",model.loop,", ",roi.loop,", ",stage.loop,", & ",patient," at ",start.time <- Sys.time())
            real.results <- use.cluster.classifiers.to.predict(
              all.test.data = data,
              test.trial.info = trial.info,
              current.patient = patient,
              trained.classifier = trained.classifiers[[roi.loop]][[stage.loop]][[patient]],
              .save.path = save.path.real.data,
              .save.name = save.name.real.data,
              shuffle.test = FALSE,
              .output = FALSE,
              .seed = 18)
            
          } # if(file.exists(paste0(save.path.real.data, save.name.real.data))){
          
          
          ##
          ## Shuffled data
          ## 
          
          save.path.shuffled.data <- 
            paste0(output.path,
                   '/data/10a_vC - decode sentence nouns - shuffle test - 50ms/',
                   model.loop,'/',
                   'shuffled data/',
                   roi.loop,'/',
                   stage.loop,'/',
                   patient,'/')
          dir.create(save.path.shuffled.data, showWarnings = FALSE, recursive = TRUE)
          save.name.shuffled.data.prefix <- paste0(roi.loop,' - ',stage.loop,' - ',patient,' - seed=')
          
          ## How many shuffle.loops?
          n.shuffles <- 1000
          shuffle.seeds.to.run <- (1:n.shuffles) * 100 + 18
          
          ## Get runs already done
          existing.shuffle.files <- list.files(save.path.shuffled.data)
          existing.shuffle.seeds <- 
            as.numeric(gsub(save.name.shuffled.data.prefix, "", 
                            gsub(".RData", "", existing.shuffle.files, fixed = TRUE), 
                            fixed = TRUE))
          shuffle.seeds.to.run <- shuffle.seeds.to.run[! shuffle.seeds.to.run %in% existing.shuffle.seeds]
          
          ### If there are any shuffle seeds to run, run them!
          if(length(shuffle.seeds.to.run > 0)){
            
            # Time
            message("Beginning classification of shuffled data for ",roi.loop,", ",stage.loop,", & ",patient," at ",start.time <- Sys.time())
            
            # ## Set up parallel processing
            # # Close any old parallel backends
            # unregister_dopar()
            # # Set up parallel workers in case caret wants to use it
            # cl <- makeCluster(n.cores.to.use, type = "FORK")
            # registerDoParallel(cl)
            # 
            # foreach(seed.loop = shuffle.seeds.to.run,
            #         .options.nws = list(chunkSize = 15)) %dopar% {
            
            for(seed.loop in shuffle.seeds.to.run){
              # seed.loop = shuffle.seeds.to.run[1]
              
              use.cluster.classifiers.to.predict(
                all.test.data = data,
                test.trial.info = trial.info,
                current.patient = patient,
                trained.classifier = trained.classifiers[[roi.loop]][[stage.loop]][[patient]],
                shuffle.test = TRUE,
                .output = FALSE,
                .save.path = save.path.shuffled.data,
                .save.name = paste0(save.name.shuffled.data.prefix,seed.loop,'.RData'),
                .seed = seed.loop)
              
            }; rm(seed.loop)
            # } # foreach
            
            # # End parallel processing
            # stopCluster(cl)
            # unregister_dopar()
            
          } # if(length(shuffle.seeds.to.run > 0)){
          
        }; rm(stage.loop)
      } # if(patient %in% current.cluster.patients){ 
      
      ### Progress update
      progress.update.path <- 
        paste0(output.path,
               'data/10a_vC - decode sentence nouns - shuffle test - 50ms/progress updates - completed runs (can delete)/',
               model.loop,'/')
      dir.create(progress.update.path, showWarnings = FALSE, recursive = TRUE)
      write.csv(data.frame('model.loop' = model.loop,
                           'roi.loop' = roi.loop,
                           'patient' = patient,
                           'n.shuffles.completed' = ifelse(patient %in% current.cluster.patients, n.shuffles, 'not run - not enough elecs this patient this region'),
                           'duration' = ifelse(patient %in% current.cluster.patients, Sys.time() - start.time, "N/A"),
                           'start.time' = ifelse(patient %in% current.cluster.patients, start.time, "N/A"),
                           'completed.time' = Sys.time()),
                file = paste0(progress.update.path, roi.loop,' - ',patient,'.csv'),
                row.names = FALSE, quote = FALSE)
      
      
    } # foreach roi.loop %:% patient
  
  # End parallel processing
  stopCluster(cl)
  unregister_dopar()

  # }#; rm(patient)
  #}#; rm(roi.loop)
  
}; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


