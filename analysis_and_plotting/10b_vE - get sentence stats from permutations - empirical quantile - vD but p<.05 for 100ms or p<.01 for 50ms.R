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
    n.cores.to.use = 10
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
# Get sig windows multiple thresholds
source(paste0(path,'/analysis/R/functions/get_significant_windows_multiple_thresholds.R'))
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
# Inverse get significant windows
source(paste0(path,'/analysis/R/functions/get_significant_windows_inverse.R'))

# Define function for getting sig windows at p<.05 for 100ms or p<.01 for 50ms
get.sig.05_100ms.or.01_50ms.vals <- function(..acc, 
                                             ..05.thresh,
                                             ..01.thresh,
                                             ..sample.labels){
  # Get significant windows for each criterion
  sw.05_100 <- get.significant.windows(sig.vals = ..acc > ..05.thresh,
                                       .sample.labels = ..sample.labels,
                                       output.class = 'data.frame',
                                       include.duration = TRUE,
                                       .exclude.sig.durations.under.ms = 100,
                                       .sampling.rate = 256)
  sw.01_50 <- get.significant.windows(sig.vals = ..acc > ..01.thresh,
                                      .sample.labels = ..sample.labels,
                                      output.class = 'data.frame',
                                      include.duration = TRUE,
                                      .exclude.sig.durations.under.ms = 50,
                                      .sampling.rate = 256)
  # Get vectors of 1s and 0s per time sample for each criterion
  sig.vector.05_100 <- get.significant.windows.inverse(sig.windows = sw.05_100,
                                                       .sample.labels = ..sample.labels,
                                                       .sampling.rate = 256)
  sig.vector.01_50 <- get.significant.windows.inverse(sig.windows = sw.01_50,
                                                      .sample.labels = ..sample.labels,
                                                      .sampling.rate = 256)
  # Merge
  sig.vector.out <- (sig.vector.05_100 + sig.vector.01_50) > 0
  
  # Get sig windows!
  output <- get.significant.windows(sig.vals = sig.vector.out,
                                    .sample.labels = ..sample.labels,
                                    output.class = 'data.frame',
                                    include.duration = TRUE,
                                    .sampling.rate = 256)
  
  # Output
  return(output)
} # get.sig.05_100ms.or.01_50ms.vals


### Set seed
set.seed(2024)

### Get n.trials.per.patient
print("Attaching trial info...")
attach(paste0(path,'analysis/R/downsample data/warped data/high_gamma/warped high_gamma data at 256 Hz.RData'))
n.trials.per.patient <- sapply(patient.trial.info, function(x){nrow(x[['pn']])})
detach()
print("...detached!")



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


### Loop thru models
model.types <- list.files(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/'))
for(model.loop in model.types){ ## UNCOMMENT
  # model.loop = model.types[2] ## UNCOMMENT
  
  message('Beginning ',model.loop,' loop. ',Sys.time())
  
  ### Load classifiers
  load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/',model.loop,'/trained.classifiers.RData'))
  patients.per.roi <- lapply(trained.classifiers, function(x){names(x[[1]])})
  rm(trained.classifiers) # not the most efficient way to do this -- 100MB file
  
  
  ### Loop thru rois
  # for(roi.loop in rois){ ## UNCOMMENT
  #   # roi.loop = rois[1] ## UNCOMMENT
  
  ## Set up parallel processing
  # Close any old parallel backends
  unregister_dopar()
  # Set up parallel workers in case caret wants to use it
  cl <- makeCluster(n.cores.to.use, type = "FORK")
  registerDoParallel(cl)
  
  foreach(roi.loop = rev(rois)) %:% 
    # roi.loop = rois[2] ## UNCOMMENT
    
    # ### Loop thru stages
    # # For clean up
    # stage.keep <- c(ls(), 'stage.keep', 'stage.loop')
    
    # Loop!
    # for(stage.loop in rev(names(stage.sample.labels))){ ## UNCOMMENT
    #   # stage.loop = names(stage.sample.labels)[1] ## UNCOMMENT
    
    # # Clean up
    # rm(list = ls()[! ls() %in% stage.keep])
    # gc()
    
    foreach(stage.loop = rev(names(stage.sample.labels))) %dopar% { ## UNCOMMENT
      # stage.loop = rev(names(stage.sample.labels))[1] ## UNCOMMENT
        
        # Get patients with >=3 elecs this roi
        current.patients <- patients.per.roi[[roi.loop]]
        
        ### Save data dir
        save.data.dir <- paste0(output.path, 'data/10b_vE - get sentence stats from permutations - 50ms/',
                                model.loop,'/',
                                roi.loop,'/',
                                stage.loop,'/')
        # save.1000.shuffles.name <- list.files(save.data.dir)#'patient and combined sentence stats - 1000 shuffles.RData'
        # save.1000.shuffles.name <- save.1000.shuffles.name[grepl("patient and combined sentence stats", save.1000.shuffles.name)]
        save.1000.shuffles.name <- 'patient and combined sentence stats - 1000 shuffles.RData'
        
        ### If data already saved, load and skip to plotting
        if(file.exists(paste0(save.data.dir, save.1000.shuffles.name))){
          
          message('Stats for ',model.loop,' - ',roi.loop, ' - ', stage.loop,' already calculated! Loading data for plotting.')
          
          ## Load data
          load(paste0(save.data.dir, save.1000.shuffles.name))
          
        }else{
          
          # Storage
          real.accs <- list()
          shuff.accs <- list()
          patient.sl.stats <- list() # sentence list stats
          
          ### Loop thru patients
          for(patient in current.patients){ ## UNCOMMENT
            # patient = current.patients[1]  ## UNCOMMENT
            
            message('Starting: ',model.loop,' - ',band.loop,' - ',roi.loop,' - ',stage.loop,' - ',patient,'. ',Sys.time())
            
            real.accs.dir <- 
              paste0(output.path,
                     'data/10a - decode sentence nouns - shuffle test - 50ms/',
                     model.loop,
                     '/real data/',
                     roi.loop,'/',
                     stage.loop,'/',
                     patient,'/')
            shuffled.accs.dir <- 
              paste0(output.path,
                     'data/10a - decode sentence nouns - shuffle test - 50ms/',
                     model.loop,
                     '/shuffled data/',
                     roi.loop,'/',
                     stage.loop,'/',
                     patient,'/')
            
            
            ### Load real data
            print("Attaching real accuracies...")
            attach(paste0(real.accs.dir,
                          roi.loop,' - ',
                          stage.loop,' - ',
                          patient,
                          '.RData'))
            real.accs[[patient]] <- accuracies.by.lx.class
            detach()
            print("...detached!")
            
            # Clean up
            real.accs[[patient]] <- lapply(real.accs[[patient]], function(x){
              rownames(x) <- x$sample.label
              x <- x[,-which(colnames(x) == 'sample.label')]
              return(x)
            })
            
            
            ### Load shuffled accuracies
            shuff.acc.files <- 
              list.files(shuffled.accs.dir)
            
            shuff.accs[[patient]] <- list()
            print("Attaching shuffle accuracies...")
            for(shuff.acc.file in shuff.acc.files){
              attach(paste0(shuffled.accs.dir, shuff.acc.file))
              shuff.accs[[patient]][[shuff.acc.file]] <- accuracies.by.lx.class
              detach()
            }; rm(shuff.acc.file, shuff.acc.files)
            print("...detached!")
            
            # Clean up
            shuff.accs[[patient]] <- lapply(shuff.accs[[patient]], function(y){
              lapply(y, function(x){
                rownames(x) <- x$sample.label
                x <- x[,-which(colnames(x) == 'sample.label')]
                return(x)  
              })})
            
            
            ### Exclude patient NY799 from passive analyses (only completed half of one trial which is doing weird things to classification accuracies)
            if(patient == "NY799"){
              
              passive.columns <- c('subject.noun.passive', 'object.noun.passive', 'verb.passive')
              
              for(..noun.loop in c('noun1', 'noun2')){
                real.accs[[patient]][[..noun.loop]][,passive.columns] <- NA
                shuff.accs[[patient]] <- lapply(shuff.accs[[patient]], function(x){
                  x[[..noun.loop]][,passive.columns] <- NA
                  return(x)
                })
              };rm(..noun.loop)
              
              rm(passive.columns)
            } # if(patient == "NY799"){
            
            
            ### Smooth
            half.n.smoothing.samples.sl <- time.convert(50, "times", "samples", 256) # 100ms window (in half)
            # real.accs[[patient]] <- lapply(real.accs[[patient]], function(y){
            #   data.frame(apply(y, 2, function(x){
            #     smoothing(x,
            #               n.samples.pre = half.n.smoothing.samples.sl,
            #               na.pad = FALSE)
            #   }))})
            # 
            # shuff.accs[[patient]] <- lapply(shuff.accs[[patient]], function(z){
            #   lapply(z, function(y){
            #     data.frame(apply(y, 2, function(x){
            #       smoothing(x,
            #                 n.samples.pre = half.n.smoothing.samples.sl,
            #                 na.pad = FALSE)
            #     }))})})
            
            ### Load trial information (to get number of trials per role for binomial distribution/Bayesian analyses)
            trial.info <- read.csv(paste0(path,'/data/',patient,'/data/epoched/',patient,'_trial_labels_without_bad_trials.csv'))
            sp.info <- just.sentence.noun.verb.data(trial.info)
            sp.info[is.na(sp.info$case), "case"] <- 'none'
            lp.info <- just.list.noun.data(trial.info)
            
            ## Get number of trials per linguistic class
            n.trials.per.role <- list()
            
            # Active sentences
            n.trials.per.role[['subject.noun.active']] <-
              nrow(sp.info[(sp.info$verb_voice == "active") & (sp.info$case == "s"),])
            n.trials.per.role[['verb.active']] <-
              nrow(sp.info[(sp.info$verb_voice == "active") & (sp.info$pos == "v"),])
            n.trials.per.role[['object.noun.active']] <-
              nrow(sp.info[(sp.info$verb_voice == "active") & (sp.info$case == "do"),])
            
            # Passive sentences
            n.trials.per.role[['subject.noun.passive']] <-
              nrow(sp.info[(sp.info$verb_voice == "passive") & (sp.info$case == "s"),])
            n.trials.per.role[['verb.passive']] <-
              nrow(sp.info[(sp.info$verb_voice == "passive") & (sp.info$pos == "v"),])
            n.trials.per.role[['object.noun.passive']] <-
              nrow(sp.info[(sp.info$verb_voice == "passive") & (sp.info$case == "bo"),])
            
            # Sentences (collapsing across active/passive)
            n.trials.per.role[['subject.noun']] <-
              nrow(sp.info[sp.info$case == "s",])
            n.trials.per.role[['verb']] <-
              nrow(sp.info[sp.info$pos == "v",])
            n.trials.per.role[['object.noun']] <-
              nrow(sp.info[sp.info$case %in% c("do", "bo"),])
            n.trials.per.role[['sentence.noun']] <-
              nrow(sp.info[sp.info$pos == "n",])
            
            # Lists
            n.trials.per.role[['list1']] <-
              nrow(lp.info[lp.info$case == "1",])
            n.trials.per.role[['list2']] <-
              nrow(lp.info[lp.info$case == "2",])
            n.trials.per.role[['list.noun']] <-
              nrow(lp.info[lp.info$case %in% c("1", "2"),])
            
            
            ### Stats for sentences and lists ("sl")
            patient.sl.stats[[patient]] <- list()
            for(lx.loop in colnames(real.accs[[patient]][[1]])){
              # lx.loop = colnames(real.accs[[patient]][[1]])[1]
              
              patient.sl.stats[[patient]][[lx.loop]] <- list()
              for(word.loop in names(real.accs[[patient]])){
                # word.loop = names(real.accs[[patient]])[1]
                
                ## Summary stats
                temp <- data.frame('sample.label' = rownames(real.accs[[patient]][[word.loop]]),
                                   'accuracy' = real.accs[[patient]][[word.loop]][,lx.loop],
                                   'shuffle.mean' = apply(bind_rows(lapply(shuff.accs[[patient]], function(x){
                                     x[[word.loop]][,lx.loop]})), 1, mean),
                                   'shuffle.95ci.upper' = apply(bind_rows(lapply(shuff.accs[[patient]], function(x){
                                     x[[word.loop]][,lx.loop]})), 1, quantile, .95, na.rm = TRUE),
                                   'shuffle.95ci.lower' = apply(bind_rows(lapply(shuff.accs[[patient]], function(x){
                                     x[[word.loop]][,lx.loop]})), 1, quantile, .05, na.rm = TRUE),
                                   'shuffle.99ci.upper' = apply(bind_rows(lapply(shuff.accs[[patient]], function(x){
                                     x[[word.loop]][,lx.loop]})), 1, quantile, .99, na.rm = TRUE),
                                   'shuffle.975ci.lower' = apply(bind_rows(lapply(shuff.accs[[patient]], function(x){
                                     x[[word.loop]][,lx.loop]})), 1, quantile, .025, na.rm = TRUE),
                                   'shuffle.975ci.upper' = apply(bind_rows(lapply(shuff.accs[[patient]], function(x){
                                     x[[word.loop]][,lx.loop]})), 1, quantile, .975, na.rm = TRUE))
                
                ## Add Bayes Factor
                # Add column
                temp$bayes.factor.log10 <- NA
                if(! any(c(is.na(temp$shuffle.mean), is.na(temp$accuracy)))){
                  # Null hypothesis (conservative: chance is defined as whichever's higher: theoretical (1/6) or empirical (mean of shuffle accuracies))
                  p0 <- sapply(temp$shuffle.mean, function(x){max(c(x, 1/6))})
                  n.trials <- n.trials.per.role[[lx.loop]]
                  n.correct.predictions <- round(temp$accuracy * n.trials)
                  for(row.loop in 1:nrow(temp)){
                    # row.loop = 1
                    
                    # Get BF using default priors (Jeffreys prior: Beta(.5, .5))
                    temp$bayes.factor.log10[row.loop] <- 
                      log10(extractBF(proportionBF(y = n.correct.predictions[row.loop], 
                                                   N = n.trials, 
                                                   p = p0[row.loop],
                                                   nullInterval = c(p0, 1)))[1,"bf"]) # test: y greater than p?
                  }; rm(row.loop)
                  rm(p0, n.trials, n.correct.predictions)
                } # If any non-NA data
                
                ## Smooth
                for(smooth.loop in c('accuracy','shuffle.mean','shuffle.95ci.upper','shuffle.95ci.lower','shuffle.99ci.upper','shuffle.975ci.lower','shuffle.975ci.upper','bayes.factor.log10')){
                  temp[,smooth.loop] <- smoothing(temp[,smooth.loop], n.samples.pre = half.n.smoothing.samples.sl)  
                }; rm(smooth.loop)
                
                # Remove NAs introduced by smoothing
                n.rows <- nrow(temp)
                temp <- temp[-c((1:half.n.smoothing.samples.sl),
                                ((n.rows - half.n.smoothing.samples.sl + 1):n.rows)),]
                rm(n.rows)
                
                # # Analystical stats
                # temp$z.score <- 
                #   with(temp, 
                #        (accuracy - shuffle.mean) / shuffle.sd)
                # temp$p.value <- 
                #   1 - pnorm(temp$z.score) # one-tailed
                
                # Store
                patient.sl.stats[[patient]][[lx.loop]][[word.loop]] <- temp
                rm(temp)
                
              }; rm(word.loop)
            }; rm(lx.loop)
            
          }; rm(patient)
          
          n.shuffles.done <- min(sapply(shuff.accs, length))
          print(paste0(roi.loop, ' - ',stage.loop, ': All patients have at least ',n.shuffles.done,' shuffles completed.'))
          
          ### Get sig windows
          patient.sl.sig.windows <- 
            lapply(patient.sl.stats, function(z){
              lapply(z, function(y){
                lapply(y, function(x){
                  get.sig.05_100ms.or.01_50ms.vals(..acc = x$accuracy,
                                                   ..05.thresh = x$shuffle.95ci.upper,
                                                   ..01.thresh = x$shuffle.99ci.upper,
                                                   ..sample.labels = x$sample.label)})})})
          
          ### Get Bayesian "sig" windows
          # Thresholded at BF > 3
          bf3.threshold <- log10(3)
          patient.sl.BF.over.3.windows <-
            lapply(patient.sl.stats, function(z){
              lapply(z, function(y){
                lapply(y, function(x){
                  get.significant.windows(x$bayes.factor.log10 > bf3.threshold,
                                          .sample.labels = x$sample.label,
                                          output.class = 'data.frame',
                                          .exclude.sig.durations.under.ms = 100,
                                          .sampling.rate = 256,
                                          include.duration = TRUE)})})})
          
          # Thresholded at BF > 3
          bf10.threshold <- log10(10)
          patient.sl.BF.over.10.windows <-
            lapply(patient.sl.stats, function(z){
              lapply(z, function(y){
                lapply(y, function(x){
                  get.significant.windows(x$bayes.factor.log10 > bf10.threshold,
                                          .sample.labels = x$sample.label,
                                          output.class = 'data.frame',
                                          .exclude.sig.durations.under.ms = 100,
                                          .sampling.rate = 256,
                                          include.duration = TRUE)})})})
          
          
          ##
          ## Combine across patients
          ##
          
          message("Starting stats on combined patients. ", Sys.time())
          
          ### Real data
          wm <- function(z){weighted.average(z, w = n.trials.per.patient[names(real.accs)])}
          combined.accs <- list()
          for(word.loop in names(real.accs[[1]])){
            # word.loop = names(real.accs[[1]])[1]
            
            combined.accs[[word.loop]] <- elementwise.matrix.apply(
              lapply(real.accs, function(x){x[[word.loop]]}),
              .function = 'wm')
            
          }; rm(word.loop)
          
          
          ### Shuffled data
          wm <- function(z){weighted.average(z, w = n.trials.per.patient[names(shuff.accs)])}
          combined.shuffs <- list()
          for(shuffle.loop in 1:min(sapply(shuff.accs, length))){
            # shuffle.loop = 1
            
            combined.shuffs[[shuffle.loop]] <- list()
            for(word.loop in names(shuff.accs[[1]][[1]])){
              # word.loop = names(shuff.accs[[1]][[1]])[1]
              
              combined.shuffs[[shuffle.loop]][[word.loop]] <- 
                data.frame(elementwise.matrix.apply(
                  lapply(shuff.accs, function(x){
                    x[[shuffle.loop]][[word.loop]]}),
                  .function = 'wm'))
              
            }; rm(word.loop)
          }; rm(shuffle.loop)
          
          
          ### Stats for sentences and lists ("sl")
          combined.sl.stats <- list()
          for(lx.loop in colnames(combined.accs[[1]])){
            # lx.loop = colnames(combined.accs[[1]])[1]
            
            combined.sl.stats[[lx.loop]] <- list()
            for(word.loop in names(combined.accs)){
              # word.loop = names(combined.accs)[1]
              
              # Summary stats
              temp <- data.frame('sample.label' = rownames(combined.accs[[word.loop]]),
                                 'accuracy' = combined.accs[[word.loop]][,lx.loop],
                                 'shuffle.mean' = apply(do.call(cbind, lapply(combined.shuffs, function(x){
                                   x[[word.loop]][,lx.loop]})), 1, mean),
                                 'shuffle.95ci.upper' = apply(do.call(cbind, lapply(combined.shuffs, function(x){
                                   x[[word.loop]][,lx.loop]})), 1, quantile, .95, na.rm = TRUE),
                                 'shuffle.95ci.lower' = apply(do.call(cbind, lapply(combined.shuffs, function(x){
                                   x[[word.loop]][,lx.loop]})), 1, quantile, .05, na.rm = TRUE),
                                 'shuffle.99ci.upper' = apply(do.call(cbind, lapply(combined.shuffs, function(x){
                                   x[[word.loop]][,lx.loop]})), 1, quantile, .99, na.rm = TRUE),
                                 'shuffle.975ci.lower' = apply(do.call(cbind, lapply(combined.shuffs, function(x){
                                   x[[word.loop]][,lx.loop]})), 1, quantile, .025, na.rm = TRUE),
                                 'shuffle.975ci.upper' = apply(do.call(cbind, lapply(combined.shuffs, function(x){
                                   x[[word.loop]][,lx.loop]})), 1, quantile, .975, na.rm = TRUE))
              
              
              # Smooth
              for(smooth.loop in c('accuracy','shuffle.mean','shuffle.95ci.upper','shuffle.95ci.lower','shuffle.99ci.upper','shuffle.975ci.lower','shuffle.975ci.upper')){
                temp[,smooth.loop] <- smoothing(temp[,smooth.loop], n.samples.pre = half.n.smoothing.samples.sl)  
              }; rm(smooth.loop)
              
              # Remove NAs introduced by smoothing
              n.rows <- nrow(temp)
              temp <- temp[-c((1:half.n.smoothing.samples.sl),
                              ((n.rows - half.n.smoothing.samples.sl + 1):n.rows)),]
              rm(n.rows)
              
              # # Analytical stats
              # temp$z.score <- 
              #   with(temp, 
              #        (accuracy - shuffle.mean) / shuffle.sd)
              # temp$p.value <- 
              #   1 - pnorm(temp$z.score) # one-tailed
              
              # Store
              combined.sl.stats[[lx.loop]][[word.loop]] <- temp
              rm(temp)
              
            }; rm(word.loop)
          }; rm(lx.loop)
          
          
          ### Get sig windows
          combined.sl.sig.windows <- 
            lapply(combined.sl.stats, function(x){
              lapply(x, function(y){
                get.sig.05_100ms.or.01_50ms.vals(..acc = y$accuracy,
                                                 ..05.thresh = y$shuffle.95ci.upper,
                                                 ..01.thresh = y$shuffle.99ci.upper,
                                                 ..sample.labels = y$sample.label)
              })})
          
          
          ### Save
          save.these <- c('patient.sl.stats',
                          'patient.sl.sig.windows',
                          'patient.sl.BF.over.3.windows',
                          'patient.sl.BF.over.10.windows',
                          'combined.sl.stats',
                          'combined.sl.sig.windows',
                          'half.n.smoothing.samples.sl')
          dir.create(save.data.dir, showWarnings = FALSE, recursive = TRUE)
          save(list = save.these,
               file = paste0(save.data.dir, 'patient and combined sentence stats - ',n.shuffles.done,' shuffles.RData'))
          
        } # if stats already saved
        
        
        ##
        ## Quick and dirty plot of each word by patient
        ##
        
        message("Plotting - quick and dirty.", Sys.time())
        
        # Save dir
        save.patient.plot.dir <- paste0(output.path, 'figures/10b_vE - get sentence stats from permutations - 50ms/',
                                        'quick and dirty time series/single words by patient/')
        dir.create(save.patient.plot.dir, showWarnings = FALSE, recursive = TRUE)
        
        # Data
        for(lx.loop in names(patient.sl.stats[[1]])){
          # lx.loop = names(patient.sl.stats[[1]])[1]
          
          # Set up save image
          pdf(paste0(save.patient.plot.dir, 'accuracy - ',roi.loop,' - ',stage.loop,' - ',lx.loop,' - ',model.loop,'.pdf'),
              width = 14, height = 6)
          par(mfrow = c(1,2))
          
          for(word.loop in names(patient.sl.stats[[1]][[lx.loop]])){
            # word.loop = names(patient.sl.stats[[1]][[lx.loop]])[1]
            
            # Figure data
            ys <- lapply(patient.sl.stats, function(x){x[[lx.loop]][[word.loop]]$accuracy})
            ys[['combined']] <- combined.sl.stats[[lx.loop]][[word.loop]]$accuracy
            sig.windows <- lapply(patient.sl.sig.windows, function(x){x[[lx.loop]][[word.loop]]})
            sig.windows[['combined']] <- combined.sl.sig.windows[[lx.loop]][[word.loop]]
            
            # Plot
            plot.time.series(.y.values = ys,
                             .x.values = combined.sl.stats[[lx.loop]][[word.loop]]$sample.label,
                             .sampling.rate = 256,
                             .colors = c(grey.colors(length(patient.sl.stats)), 'magenta'),
                             .y.lwd = c(rep(3, times = length(patient.sl.stats)), 6),
                             .y.limits.min.at.least = .08,
                             .y.limits.max.at.least = .4,
                             .y.ticks = c(.1,.2,.3,.4),
                             .y.label = 'accuracy',
                             .shuffle.dist.mean = combined.sl.stats[[lx.loop]][[word.loop]]$shuffle.mean,
                             .shuffle.dist.error.bars.upper = combined.sl.stats[[lx.loop]][[word.loop]]$shuffle.975ci.upper - combined.sl.stats[[lx.loop]][[word.loop]]$shuffle.mean,
                             .shuffle.dist.error.bars.lower = combined.sl.stats[[lx.loop]][[word.loop]]$shuffle.mean - combined.sl.stats[[lx.loop]][[word.loop]]$shuffle.975ci.lower,
                             .shuffle.dist.color = adjust.transparency('magenta', alpha = .3),
                             # .polygons.x = stage.ranges[[stage.loop]],
                             .polygons.x = unlist(stage.ranges[stage.loop,c('start.time','end.time')]),
                             .polygons.color = 'white',
                             .horizontal.line.at = 1/6,
                             .sig.windows = sig.windows,
                             .sig.color = c(grey.colors(length(patient.sl.stats)), 'magenta'),
                             .center.sig.bars.vertically = FALSE,
                             .title = paste0(roi.loop,' - ',stage.loop,'\n',lx.loop,' - ',word.loop))
            
          }; rm(word.loop)
          
          # Save image
          dev.off()
          
        }; rm(lx.loop)
        
        
        ### Quick and dirty plot of noun1 and noun2 averaging across patients
        
        # What data to plot: sentences? lists?
        lx.comparisons <- list(
          'sentences' = c('subject.noun','verb','object.noun'),
          'actives' = c('subject.noun.active','verb.active','object.noun.active'),
          'passives' = c('subject.noun.passive','verb.passive','object.noun.passive'),
          'lists' = c('list1','list2','list.noun')
        )
        
        for(comparison.loop in names(lx.comparisons)){
          # comparison.loop = names(lx.comparisons)[1]
          
          # Save dir
          save.nouns.1.2.plot.dir <- 
            paste0(output.path, 
                   'figures/10b_vE - get sentence stats from permutations - 50ms/',
                   'quick and dirty time series/noun1 and noun2 averaged across patients/',
                   comparison.loop,'/')
          dir.create(save.nouns.1.2.plot.dir, showWarnings = FALSE, recursive = TRUE)
          
          
          # Save!
          pdf(paste0(save.nouns.1.2.plot.dir, 'accuracy - ',comparison.loop,' - ',roi.loop,' - ',stage.loop,' - ',model.loop,'.pdf'),
              width = 7 * 3, 
              height = 6)
          
          par(mfrow = c(1, 3))
          
          for(lx.loop in lx.comparisons[[comparison.loop]]){
            # lx.loop = lx.comparisons[[comparison.loop]][1]
            
            # Figure data
            ys <- lapply(combined.sl.stats[[lx.loop]], function(x){x$accuracy})
            sig.windows <- combined.sl.sig.windows[[lx.loop]]
            
            plot.time.series(.y.values = ys,
                             .x.values = lapply(combined.sl.stats[[lx.loop]], function(x){x$sample.label}),
                             .sampling.rate = 256,
                             .colors = c('red', 'blue'),
                             .y.lwd = c(rep(3, times = length(patient.sl.stats)), 6),
                             .y.limits.min.at.least = .08,
                             .y.limits.max.at.least = .4,
                             .y.ticks = c(.1,.2,.3,.4),
                             .y.label = 'accuracy',
                             .shuffle.dist.mean = combined.sl.stats[[lx.loop]]$noun1$shuffle.mean,
                             .shuffle.dist.error.bars.upper = combined.sl.stats[[lx.loop]]$noun1$shuffle.975ci.upper - combined.sl.stats[[lx.loop]]$noun1$shuffle.mean,
                             .shuffle.dist.error.bars.lower = combined.sl.stats[[lx.loop]]$noun1$shuffle.mean - combined.sl.stats[[lx.loop]]$noun1$shuffle.975ci.lower,
                             .shuffle.dist.color = adjust.transparency('red', alpha = .1),
                             .polygons.x = unlist(stage.ranges[stage.loop,c('start.time','end.time')]),
                             .horizontal.line.at = 1/6,
                             .sig.windows = sig.windows,
                             .sig.color = c('red','blue'),
                             .center.sig.bars.vertically = FALSE,
                             .title = paste0(roi.loop,' - ',stage.loop,' - ',lx.loop))
            
          }; rm(lx.loop)
          
          dev.off()
          
        }; rm(comparison.loop)
        
        
        
        
        # Progress update
        message('Completed: ',
                model.loop,' - ',
                band.loop,' - ',
                roi.loop,' - ',
                stage.loop,'.',
                Sys.time())
        
        
    } # foreach() %:% foreach() %dopar% {
    # }#; rm(stage.loop)
  # }#; rm(roi.loop)
  
  ## Close parallel backend
  stopCluster(cl)
  unregister_dopar()
  
}; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


