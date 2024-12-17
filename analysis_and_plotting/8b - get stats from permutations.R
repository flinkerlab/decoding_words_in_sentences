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
  n.cores.to.use = 20
  n.cores.to.use.nmf = 32
}

### Lemmas
# Function to convert between samples, sample labels, and times
source(paste0(path,'/analysis/R/functions/time_convert.R'))
# Function for rolling average smoothing
source(paste0(path,'/analysis/R/functions/smoothing.R'))
# Close any old parallel backends
source(paste0(path,'/analysis/R/functions/unregister_dopar.R'))
# Plot time series
source(paste0(path,'/analysis/R/functions/plot_time_series.R'))
# Add colored text to plots
source(paste0(path,'/analysis/R/functions/add_text_line_multiple_colors.R'))
# Change color transparency
source(paste0(path,'/analysis/R/functions/adjust_transparency.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Elementwise matrix apply
source(paste0(path,'/analysis/R/functions/elementwise_matrix_apply.R'))
# Squish Bayes Factors
source(paste0(path,'/analysis/R/functions/squish_bayes_factors.R'))
# Weighted mean and standard error functions
source(paste0(path,'/analysis/R/functions/weighted_summary_stats.R'))
# Classification
source(paste0(path,'/analysis/R/functions/classify_parallel_chunks.R'))
# Inverse get significant windows
source(paste0(path,'/analysis/R/functions/get_significant_windows_inverse.R'))

# Downsample data
downsample.trial <- function(trial.data,
                             input.rate = 512,
                             output.rate = 256){
  output <- signal::resample(trial.data,
                             p = output.rate,
                             q = input.rate)
  return(output)
} # downsample.trial()

# Define function for getting sig windows at p<.05 for 100ms or p<.01 for 50ms
get.sig.05_100ms.or.01_50ms.vals <- function(..acc, 
                                             ..05.thresh,
                                             ..01.thresh,
                                             ..exclude.times.before.ms = NULL,
                                             ..sample.labels){
  # Get significant windows for each criterion
  sw.05_100 <- get.significant.windows(sig.vals = ..acc > ..05.thresh,
                                       .sample.labels = ..sample.labels,
                                       output.class = 'data.frame',
                                       .exclude.times.before.ms = ..exclude.times.before.ms,
                                       include.duration = TRUE,
                                       .exclude.sig.durations.under.ms = 100,
                                       .sampling.rate = 256)
  sw.01_50 <- get.significant.windows(sig.vals = ..acc > ..01.thresh,
                                      .sample.labels = ..sample.labels,
                                      output.class = 'data.frame',
                                      .exclude.times.before.ms = ..exclude.times.before.ms,
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
} # get.sig.05_100ms.or.01_50ms.vals()



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


### Get number of trials per patient (weights for combining accuracies within roi-windows)
print("Attaching trial info...")
attach(paste0(path,'analysis/R/downsample data/warped data/high_gamma/warped high_gamma data at 256 Hz.RData'))
n.trials.per.patient <- sapply(patient.trial.info, function(x){nrow(x[['pn']])})
detach()
print("...detached!")


### Load ROIs and stages
load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/stages and rois/stages and rois.RData'))
roi.elecs <- lapply(roi.elecs, function(.elec){
  .elec[.elec %in% use.these.elecs]
})

### Load RTs
load(paste0(path, 'analysis/R/downsample data/warped data/high_gamma/median.rt.samples 256 Hz.RData'))
median.rt.times <- time.convert(median.rt.samples, "samples", "times", 256)
rm(median.rt.samples) # risky to keep since dealing with 256 and 512 Hz datasets

### Loop thru models
model.types <- list.files(paste0(output.path, 'data/8a - cross-validated PN classification - 50ms/'))
model.keep <- c(ls(), 'model.keep', 'model.loop')
for(model.loop in model.types){ ## UNCOMMENT
  # model.loop = model.types[1] ## UNCOMMENT
  
  # Clean up
  rm(list = ls()[! ls() %in% model.keep])
  gc()
  
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
  
  foreach(roi.loop = rois) %:% 
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
    # stage.loop = rev(names(stage.sample.labels))[17] ## UNCOMMENT
    
    # Get patients with >=3 elecs this roi
    current.patients <- patients.per.roi[[roi.loop]]
    
    ### Save data dir
    save.data.dir <- paste0(output.path, 'data/8b - get stats from permutation - 50ms/',
                            model.loop,'/',
                            roi.loop,'/',
                            stage.loop,'/')
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
      patient.pn.stats <- list() # sentence list stats
      
      ### Loop thru patients
      for(patient in current.patients){
        # patient = current.patients[6]
        
        message('Starting: ',model.loop,' - ',band.loop,' - ',roi.loop,' - ',stage.loop,' - ',patient,'. ',Sys.time())
        
        real.accs.dir <- 
          paste0(output.path,
                 'data/8a - cross-validated PN classification - 50ms/',
                 model.loop,'/',
                 roi.loop,'/',
                 stage.loop,'/',
                 patient,'/',
                 'training samples averaged - unweighted/accuracies/real/')
        shuffled.accs.dir <- 
          paste0(output.path,
                 'data/8a - cross-validated PN classification - 50ms/',
                 model.loop,'/',
                 roi.loop,'/',
                 stage.loop,'/',
                 patient,'/',
                 'training samples averaged - unweighted/accuracies/shuffled/')
        
        ### Load real data
        print("Attaching real accuracies...")
        attach(paste0(real.accs.dir,
                      'classifier.accuracies_seed=18.RData'))
        real.accs[[patient]] <- classifier.accuracies
        detach()
        print("...detached!")
        
        
        ### Load shuffled accuracies
        shuff.acc.files <- 
          list.files(shuffled.accs.dir)
        
        shuff.accs[[patient]] <- list()
        print("Attaching shuffle accuracies...")
        for(shuff.acc.file in shuff.acc.files){
          attach(paste0(shuffled.accs.dir, shuff.acc.file))
          shuff.accs[[patient]][[shuff.acc.file]] <- classifier.accuracies
          detach()
        }; rm(shuff.acc.file)
        print(paste0("...detached! Read in ",length(shuff.acc.files)," shuffled data files.")); rm(shuff.acc.files)
        
        # Clean up
        shuff.accs[[patient]] <- data.frame(bind_rows(shuff.accs[[patient]]))
        
        ## Smoothing parameters
        half.n.smoothing.samples <- time.convert(75, "times", "samples", 256) # 150ms window (in half)
        
        
        # ### Load trial information (to get number of trials per role for binomial distribution/Bayesian analyses)
        # trial.info <- read.csv(paste0(path,'/data/',patient,'/data/epoched/',patient,'_trial_labels_without_bad_trials.csv'))
        # rt.quantile.range <- c(.025, .95)
        # just.good.pic.naming(x, .remove.outliers = TRUE, .rt.quantile.range = rt.quantile.range)})
        
        
        ### Summary stats
        temp <-
          data.frame('sample.label' = names(real.accs[[patient]]),
                     'accuracy' = real.accs[[patient]],
                     'shuffle.mean' = apply(shuff.accs[[patient]], 2, mean),
                     'shuffle.sd' = apply(shuff.accs[[patient]], 2, sd),
                     'shuffle.95ci.upper' = apply(shuff.accs[[patient]], 2, quantile, .95, na.rm = TRUE),
                     'shuffle.95ci.lower' = apply(shuff.accs[[patient]], 2, quantile, .05, na.rm = TRUE),
                     'shuffle.99ci.upper' = apply(shuff.accs[[patient]], 2, quantile, .99, na.rm = TRUE),
                     'shuffle.975ci.lower' = apply(shuff.accs[[patient]], 2, quantile, .025, na.rm = TRUE),
                     'shuffle.975ci.upper' = apply(shuff.accs[[patient]], 2, quantile, .975, na.rm = TRUE))
        
        ## Add Bayes Factor
        # Add column
        temp$bayes.factor.log10 <- NA
        if(! any(c(is.na(temp$shuffle.mean), is.na(temp$accuracy)))){
          # Null hypothesis (conservative: chance is defined as whichever's higher: theoretical (1/6) or empirical (mean of shuffle accuracies))
          p0 <- sapply(temp$shuffle.mean, function(x){max(c(x, 1/6))})
          n.correct.predictions <- round(temp$accuracy * n.trials.per.patient[patient])
          for(row.loop in 1:nrow(temp)){
            # row.loop = 1
            
            # Get BF using default priors (Jeffreys prior: Beta(.5, .5))
            temp$bayes.factor.log10[row.loop] <- 
              log10(extractBF(proportionBF(y = n.correct.predictions[row.loop], 
                                           N = n.trials.per.patient[patient], 
                                           p = p0[row.loop],
                                           nullInterval = c(p0, 1)))[1,"bf"]) # test: y greater than p?
          }; rm(row.loop)
          rm(p0, n.correct.predictions)
        } # If any non-NA data
        
        ## Smooth
        for(smooth.loop in c('accuracy','shuffle.mean','shuffle.sd','shuffle.95ci.upper','shuffle.95ci.lower','shuffle.99ci.upper','shuffle.975ci.lower','shuffle.975ci.upper','bayes.factor.log10')){
          temp[,smooth.loop] <- smoothing(temp[,smooth.loop], n.samples.pre = half.n.smoothing.samples)  
        }; rm(smooth.loop)
        
        # Remove NAs introduced by smoothing
        n.rows <- nrow(temp)
        temp <- temp[-c((1:half.n.smoothing.samples),
                        ((n.rows - half.n.smoothing.samples + 1):n.rows)),]
        rm(n.rows)
        
        # # Analytical stats
        # patient.pn.stats[[patient]]$z.score <- with(patient.pn.stats[[patient]], (accuracy - shuffle.mean) / shuffle.sd)
        # patient.pn.stats[[patient]]$p.value <- 1 - pnorm(patient.pn.stats[[patient]]$z.score) # one-tailed
        
        # Store
        patient.pn.stats[[patient]] <- temp
        
      }; rm(patient)
      
      
      ### Get sig windows
      patient.pn.sig.windows <- 
        lapply(patient.pn.stats, function(x){
          get.sig.05_100ms.or.01_50ms.vals(..acc = x$accuracy,
                                           ..05.thresh = x$shuffle.95ci.upper,
                                           ..01.thresh = x$shuffle.99ci.upper,
                                           ..exclude.times.before.ms = -median.rt.times['pn'],
                                           ..sample.labels = x$sample.label)})
      
      ### Get Bayesian "sig" windows
      # Thresholded at BF > 3
      bf3.threshold <- log10(3)
      patient.pn.BF.over.3.windows <-
        lapply(patient.pn.stats, function(x){
          get.significant.windows(x$bayes.factor.log10 > bf3.threshold,
                                  .sample.labels = x$sample.label,
                                  output.class = 'data.frame',
                                  .exclude.sig.durations.under.ms = 100,
                                  .exclude.times.before.ms = -median.rt.times['pn'],
                                  .sampling.rate = 256,
                                  include.duration = TRUE)})
      
      # Thresholded at BF > 3
      bf10.threshold <- log10(10)
      patient.pn.BF.over.10.windows <-
        lapply(patient.pn.stats, function(x){
          get.significant.windows(x$bayes.factor.log10 > bf10.threshold,
                                  .sample.labels = x$sample.label,
                                  output.class = 'data.frame',
                                  .exclude.sig.durations.under.ms = 100,
                                  .exclude.times.before.ms = -median.rt.times['pn'],
                                  .sampling.rate = 256,
                                  include.duration = TRUE)})
      
      
      ##
      ## Combine across patients
      ##
      
      ### Real data
      combined.accs <- data.frame(do.call(cbind, real.accs))
      combined.accs <- apply(combined.accs, 1, function(x){
        weighted.average(x, w = n.trials.per.patient[colnames(combined.accs)])
      })
      
      ### Shuffled data
      combined.shuffs <- list()
      n.shuffles.done <- min(sapply(shuff.accs, nrow))
      for(shuffle.loop in 1:n.shuffles.done){
        # shuffle.loop = 1
        combined.shuffs[[shuffle.loop]] <- data.frame(do.call(cbind, lapply(shuff.accs, function(x){unlist(x[shuffle.loop,])})))
        combined.shuffs[[shuffle.loop]] <- apply(combined.shuffs[[shuffle.loop]], 1, function(x){
          weighted.average(x, w = n.trials.per.patient[colnames(combined.shuffs[[shuffle.loop]])])
        })
      }; rm(shuffle.loop)
      
      ## Clean up
      combined.shuffs <- data.frame(bind_rows(combined.shuffs))
      
      # # Z-scores
      # combined.pn.stats <- data.frame('sample.label' = names(combined.accs),
      #                              'accuracy' = combined.accs,
      #                              'shuffle.mean' = apply(combined.shuffs, 2, mean),
      #                              'shuffle.sd' = apply(combined.shuffs, 2, sd))
      # combined.pn.stats$z.score <- with(combined.pn.stats, (accuracy - shuffle.mean) / shuffle.sd)
      # combined.pn.stats$p.value <- 1 - pnorm(combined.pn.stats$z.score) # one-tailed
      
      # Summary stats
      combined.pn.stats <- 
        data.frame('sample.label' = names(combined.accs),
                   'accuracy' = combined.accs,
                   'shuffle.mean' = apply(combined.shuffs, 2, mean, na.rm = TRUE),
                   'shuffle.sd' = apply(combined.shuffs, 2, sd, na.rm = TRUE),
                   'shuffle.95ci.upper' = apply(combined.shuffs, 2, quantile, .95, na.rm = TRUE),
                   'shuffle.95ci.lower' = apply(combined.shuffs, 2, quantile, .05, na.rm = TRUE),
                   'shuffle.99ci.upper' = apply(combined.shuffs, 2, quantile, .99, na.rm = TRUE),
                   'shuffle.975ci.lower' = apply(combined.shuffs, 2, quantile, .025, na.rm = TRUE),
                   'shuffle.975ci.upper' = apply(combined.shuffs, 2, quantile, .975, na.rm = TRUE))
      
      # # Analytical stats
      # combined.pn.stats$z.score <- with(combined.pn.stats, (accuracy - shuffle.mean) / shuffle.sd)
      # combined.pn.stats$p.value <- 1 - pnorm(combined.pn.stats$z.score) # one-tailed
      
      ## Smooth
      for(smooth.loop in c('accuracy','shuffle.mean','shuffle.sd','shuffle.95ci.upper','shuffle.95ci.lower','shuffle.99ci.upper','shuffle.975ci.lower','shuffle.975ci.upper')){
        combined.pn.stats[,smooth.loop] <- smoothing(combined.pn.stats[,smooth.loop], n.samples.pre = half.n.smoothing.samples)  
      }; rm(smooth.loop)
      
      # Remove NAs introduced by smoothing
      n.rows <- nrow(combined.pn.stats)
      combined.pn.stats <- combined.pn.stats[-c((1:half.n.smoothing.samples),
                                                ((n.rows - half.n.smoothing.samples + 1):n.rows)),]
      rm(n.rows)
      
      
      ### Get sig windows
      combined.pn.sig.windows <- 
        get.sig.05_100ms.or.01_50ms.vals(..acc = combined.pn.stats$accuracy,
                                         ..05.thresh = combined.pn.stats$shuffle.95ci.upper,
                                         ..01.thresh = combined.pn.stats$shuffle.99ci.upper,
                                         ..exclude.times.before.ms = -median.rt.times['pn'],
                                         ..sample.labels = combined.pn.stats$sample.label)
      
      
      ### Save
      save.these <- c('patient.pn.stats',
                      'patient.pn.sig.windows',
                      'combined.pn.stats',
                      'combined.pn.sig.windows',
                      'half.n.smoothing.samples')
      save.dir <- paste0(output.path, 'data/8b - get stats from permutation - 50ms/',
                         model.loop,'/',
                         roi.loop,'/',
                         stage.loop,'/')
      
      dir.create(save.data.dir, showWarnings = FALSE, recursive = TRUE)
      save(list = save.these,
           file = paste0(save.dir, 'patient and combined stats - ',n.shuffles.done,' shuffles.RData'))
      
      
      ### Quick and dirty plot
      # Save dir
      save.plot.dir <- paste0(output.path, 
                              'figures/8b - get stats from permutation - 50ms/',
                              model.loop,
                              '/quick and dirty time series/')
      dir.create(save.plot.dir, showWarnings = FALSE, recursive = TRUE)
      
      # Data
      ys <- lapply(patient.pn.stats, function(x){x$accuracy})
      ys[['combined']] <- combined.pn.stats$accuracy
      sig.windows <- patient.pn.sig.windows
      sig.windows[['combined']] <- combined.pn.sig.windows
      
      # Save!
      pdf(paste0(save.plot.dir, 'accuracy - ',roi.loop,' - ',stage.loop,'.pdf'),
          width = 7, height = 6)
      plot.time.series(.y.values = ys,
                       .x.values = combined.pn.stats$sample.label,
                       .sampling.rate = 256,
                       .colors = c(grey.colors(length(patient.pn.stats)), 'magenta'),
                       .y.lwd = c(rep(3, times = length(patient.pn.stats)), 6),
                       .y.limits.min.at.least = .08,
                       .y.limits.max.at.least = .4,
                       .y.ticks = c(.1,.2,.3,.4),
                       .y.label = 'accuracy',
                       .shuffle.dist.mean = combined.pn.stats$shuffle.mean,
                       .shuffle.dist.error.bars = combined.pn.stats$shuffle.sd * qnorm(.95),
                       .shuffle.dist.color = adjust.transparency('magenta', alpha = .3),
                       .polygons.x = unlist(stage.ranges[stage.loop,c('start.time','end.time')]),
                       .polygons.color = 'white',
                       .horizontal.line.at = 1/6,
                       .sig.windows = sig.windows,
                       .sig.color = c(grey.colors(length(patient.pn.stats)), 'magenta'),
                       .center.sig.bars.vertically = FALSE,
                       .title = paste0(band.loop,' - ',roi.loop,' - ',stage.loop))
      dev.off()
      
      message('Completed: ',band.loop,' - ',roi.loop,' - ',stage.loop,'. ',Sys.time())
      
    } # if file exists {} else {
  } # foreach
  
}#; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())






