### 8g - mean accuracy time series
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
source(paste0(path,'/analysis/R/functions/plot_axes.R'))
# Add colored text to plots
source(paste0(path,'/analysis/R/functions/add_text_line_multiple_colors.R'))
# Adjust color transparency
source(paste0(path,'/analysis/R/functions/adjust_transparency.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Get sig windows multiple thresholds
source(paste0(path,'/analysis/R/functions/get_significant_windows_multiple_thresholds.R'))
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))
# Plot time series stack with significance
source(paste0(path,'/analysis/R/functions/plot_time_series_stack.R'))
# Violin plot
source(paste0(path,'/analysis/R/functions/plot_violin.R'))



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

### Load ROIs and stages # OUT OF ORDER -- I guess switch 8 and 9 eventually?
load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/stages and rois/stages and rois.RData'))


### Load RTs
load(paste0(path, 'analysis/R/downsample data/warped data/high_gamma/median.rt.samples 256 Hz.RData'))
median.rt.times <- time.convert(median.rt.samples, "samples", "times", 256)
rm(median.rt.samples) # risky to keep since dealing with 256 and 512 Hz datasets

### Load colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))
load(paste0(output.path, 'data/0a - definitions/character colors/character colors.RData'))
load(paste0(output.path, 'data/0a - definitions/roi metadata/ROIs and temporal lobe splits.RData'))
load(paste0(output.path, 'data/0a - definitions/noun colors/noun colors.RData'))


time.series <- 
  time.series.by.stage <-
  sig.windows <-
  train.time.dfs <- list()


### Loop thru models
model.types <- list.files(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/'))
for(model.loop in model.types){ ## UNCOMMENT
  # model.loop = model.types[1] ## UNCOMMENT
  
  
  ###
  ### Get data: Peak accuracy vs. central train time
  ###
  
  # Storage
  time.series.by.stage[[model.loop]] <- list()
  time.series[[model.loop]] <- 
    sig.windows[[model.loop]] <-
    train.time.dfs[[model.loop]] <-
    list('collapsed.across.patients' = list(),
         'individually.by.patient' = list())
  train.time.df.template <- data.frame('start.time' = NA,
                                       'end.time' = NA)
  
  ### Loop thru rois
  for(roi.loop in rev(rois)){ ## UNCOMMENT
    # roi.loop = rois[1] ## UNCOMMENT
    
    ### Loop thru stages
    # For clean up
    stage.keep <- c(ls(), 'stage.keep', 'stage.loop')
    
    # Storage 
    time.series.by.stage[[model.loop]][[roi.loop]] <- list()
    
    # Loop!
    for(stage.loop in rev(names(stage.sample.labels))){ ## UNCOMMENT
      # stage.loop = names(stage.sample.labels)[11] ## UNCOMMENT
      
      # Clean up
      rm(list = ls()[! ls() %in% stage.keep])
      gc()
      
      # Storage 
      time.series.by.stage[[model.loop]][[roi.loop]][[stage.loop]] <- list()
      
      # Load stats ("combined.pn.sig.windows","combined.pn.stats","half.n.smoothing.samples","patient.pn.sig.windows","patient.pn.stats")
      current.data.path <- paste0(output.path, 
                                  'data/8b - get stats from permutation - 50ms/',
                                  model.loop,'/',
                                  roi.loop,'/',
                                  stage.loop,'/')
      # current.data.file <- current.data.file[grepl("patient and",current.data.file)]
      current.data.file <- 'patient and combined stats - 1000 shuffles.RData'
      # current.data.file <- current.data.file[grepl("patient and",current.data.file)]
      if(! grepl("1000",current.data.file)){message("WARNING! You don't have 1000 shuffles for ",model.loop," - ",roi.loop," - ",stage.loop,". Loading what you do have...")}
      load(paste0(current.data.path, current.data.file))
      
      
      ### Collapsed across patients
      ## If any sig windows
      if(nrow(combined.pn.sig.windows) > 0){
        
        ## Get data
        collapsed.label <- paste0(roi.loop,'.',stage.loop)
        
        # Stats
        time.series[[model.loop]][['collapsed.across.patients']][[collapsed.label]] <- combined.pn.stats
        sig.windows[[model.loop]][['collapsed.across.patients']][[collapsed.label]] <- combined.pn.sig.windows
        
        # Train times
        train.time.dfs[[model.loop]][['collapsed.across.patients']][[collapsed.label]] <- train.time.df.template
        train.time.dfs[[model.loop]][['collapsed.across.patients']][[collapsed.label]]$start.time <- 
          stage.ranges[stage.loop, 'start.time']
        train.time.dfs[[model.loop]][['collapsed.across.patients']][[collapsed.label]]$end.time <- 
          stage.ranges[stage.loop, 'end.time']
        
        # Clean up
        rm(collapsed.label)
        
      } # if(nrow(combined.pn.sig.windows) > 0){
      
      ### By patient
      for(patient in names(patient.pn.stats)){
        # patient = names(patient.pn.stats)[2]
        
        time.series.by.stage[[model.loop]][[roi.loop]][[stage.loop]][[patient]] <- max(patient.pn.stats[[patient]]$accuracy)
        
        ## If any sig windows
        if(nrow(patient.pn.sig.windows[[patient]]) > 0){
          
          ## Get data
          individual.label <- paste0(roi.loop,'.',stage.loop,'.',patient)
          
          # Stats
          time.series[[model.loop]][['individually.by.patient']][[individual.label]] <- patient.pn.stats[[patient]]
          sig.windows[[model.loop]][['individually.by.patient']][[individual.label]] <- patient.pn.sig.windows[[patient]]
          
          # Train times
          train.time.dfs[[model.loop]][['individually.by.patient']][[individual.label]] <- train.time.df.template
          train.time.dfs[[model.loop]][['individually.by.patient']][[individual.label]]$start.time <- 
            stage.ranges[stage.loop, 'start.time']
          train.time.dfs[[model.loop]][['individually.by.patient']][[individual.label]]$end.time <- 
            stage.ranges[stage.loop, 'end.time']
          
          # Clean up
          rm(individual.label)
          
          
          ### ROI plot data
          roi.label <- paste0(stage.loop,'.',patient)
          
          # Stats
          time.series[[model.loop]][[roi.loop]][[roi.label]] <- patient.pn.stats[[patient]]
          sig.windows[[model.loop]][[roi.loop]][[roi.label]] <- patient.pn.sig.windows[[patient]]
          
          
          # Train times
          train.time.dfs[[model.loop]][[roi.loop]][[roi.label]] <- train.time.df.template
          train.time.dfs[[model.loop]][[roi.loop]][[roi.label]]$start.time <- 
            stage.ranges[stage.loop, 'start.time']
          train.time.dfs[[model.loop]][[roi.loop]][[roi.label]]$end.time <- 
            stage.ranges[stage.loop, 'end.time']
          
          # Clean up
          rm(roi.label)
          
        } # if(nrow(patient.pn.sig.windows[[patient]]) > 0){
        
      }; rm(patient)
    }; rm(stage.loop)
  }; rm(roi.loop)
  
  
  ### Clean up
  train.time.dfs[[model.loop]] <- 
    lapply(train.time.dfs[[model.loop]], function(x){
      x <- bind_rows(x, .id = 'label')
      rownames(x) <- x$label
      x$label <- NULL
      return(x)
    })
  
}; rm(model.loop)
message('Done getting data! ',Sys.time()) ## < 1 min




### Plot parameters
save.plot.dir.main <- paste0(output.path, 'figures/8g - mean accuracy time series/')
zoom <- 1.8
text.size.big <- 1.8
text.size.med <- 1.6
text.size.small <- 1.4


### ROI Plots

for(model.loop in model.types){
  
  ##
  ## Plot: Accuracy - individual plots
  ##
  
  for(theme.loop in c('white','black')){ ## UNCOMMENT
    # theme.loop = c('white','black')[1] ## UNCOMMENT
    
    for(data.loop in names(time.series[[model.loop]])){
      # data.loop = names(time.series[[model.loop]])[3]
      
      # Save dir
      save.dir.all.rois.all.stages <- 
        paste0(save.plot.dir.main,
               'mean accuracy - indiviudal ROIs (bad shuffle error dists fix!)/',
               theme.loop,'/')
      dir.create(save.dir.all.rois.all.stages, showWarnings = FALSE, recursive = TRUE)
      
      # Get data
      current.ys <- apply(bind_cols(lapply(time.series[[model.loop]][[data.loop]], function(.classifier){.classifier$accuracy})), 1, mean, na.rm = TRUE)
      current.xs <- time.series[[model.loop]][[data.loop]][[1]]$sample.label
      current.shuffle.mean <- 
        apply(bind_cols(lapply(time.series[[model.loop]][[data.loop]], function(.classifier){.classifier$shuffle.mean})), 1, mean, na.rm = TRUE)
      current.shuffle.upper <- 
        apply(bind_cols(lapply(time.series[[model.loop]][[data.loop]], function(.classifier){.classifier$shuffle.95ci.upper})), 1, mean, na.rm = TRUE) - 
        current.shuffle.mean
      current.shuffle.lower <- 
        current.shuffle.mean - 
        apply(bind_cols(lapply(time.series[[model.loop]][[data.loop]], function(.classifier){.classifier$shuffle.95ci.lower})), 1, mean, na.rm = TRUE)
      
      # # Save
      pdf(paste0(save.dir.all.rois.all.stages, 
                 data.loop,
                 ' - picture naming cross validation accuracies - averaged - ',
                 model.loop,'.pdf'),
          width = 6, height = 4.5)
      par(oma = c(0,0,0,0))
      plot.time.series(.y.values = current.ys,
                       .x.values = current.xs,
                       .sampling.rate = 256,
                       .x.limits = c(-median.rt.times['pn'], 500),
                       .x.ticks = c(-median.rt.times['pn'], 0, 500),
                       .y.limits = c(.1, .25),
                       .y.ticks = c(.1, .25),
                       .shuffle.dist.mean = current.shuffle.mean,
                       .shuffle.dist.error.bars.upper = current.shuffle.upper,
                       .shuffle.dist.error.bars.lower = current.shuffle.lower,
                       .title = '',
                       .y.label = '',
                       .x.label = '',
                       .margin = c(3,3,1,1),
                       .background = rgb(1,1,1,0),
                       .theme = theme.loop,
                       .zoom = zoom)
      
      
      dev.off()
      
    }#; rm(data.loop)
  }#; rm(theme.loop)
  
  
  
  ##
  ## Plot: Mean accuracy - all ROIs one plot
  ##
  
  # Get data
  current.ys <- lapply(time.series[[model.loop]][rois], function(.roi){apply(bind_cols(lapply(.roi, function(.classifier){.classifier$accuracy})), 1, mean, na.rm = TRUE)})
  current.xs <- time.series[[model.loop]][[1]][[1]]$sample.label
  current.sig <- lapply(time.series[[model.loop]][rois], function(.roi){
    apply(bind_cols(lapply(.roi, function(.classifier){
      # Get a "z"score: vals over 1 significant (over 99% CI), under 1 non sig
      .z.ish <- (.classifier$accuracy - .classifier$shuffle.mean) / 
        (.classifier$shuffle.99ci.upper - .classifier$shuffle.mean)
      return(.z.ish)
    })), 1, max, na.rm = TRUE) > 1}) # get the max -- if any sig
  current.sig.windows <- lapply(current.sig, function(.sig){
    get.significant.windows(.sig,
                            .sample.labels = current.xs,
                            output.class = "data.frame",
                            include.duration = TRUE,
                            .exclude.sig.durations.under.ms = 100,
                            .sampling.rate = 256,
                            .exclude.times.before.ms = -median.rt.times['pn'] + 50)
  })
  
  
  # # Smooth
  # current.ys <- lapply(current.ys, smoothing, n.samples.pre = time.convert(25, 'times', 'samples', 256))
  
  # Order by order of when they come online
  roi.order <- names(sort(sapply(current.sig.windows, function(x){min(x$start.time)})))
  current.ys <- current.ys[roi.order]
  current.sig.windows <- current.sig.windows[roi.order]
  
  for(theme.loop in c('black','white')){ ## UNCOMMENT
    # theme.loop = c('white','black')[1] ## UNCOMMENT
    
    # Save dir
    save.dir.all.rois.all.stages <- 
      paste0(save.plot.dir.main,
             'mean accuracy - all ROIs/',
             theme.loop,'/')
    dir.create(save.dir.all.rois.all.stages, showWarnings = FALSE, recursive = TRUE)
    
    # # Save
    pdf(paste0(save.dir.all.rois.all.stages, 
               'picture naming cross validation accuracies - averaged - ',
               model.loop,'.pdf'),
        width = 7, height = 5.5)
    par(oma = c(0,0,1,0))
    plot.time.series(.y.values = current.ys,
                     .x.values = current.xs,
                     .sampling.rate = 256,
                     .x.limits = c(-median.rt.times['pn'], 500),
                     .x.ticks = c(-median.rt.times['pn'], 0, 500),
                     # .horizontal.line.at = 1/6,
                     .y.limits = c(.15, .25),
                     .y.ticks = c(.15, .25),
                     .colors = roi.colors$rostral.caudal[names(current.ys)],
                     .sig.windows = current.sig.windows,
                     .sig.color = roi.colors$rostral.caudal[names(current.ys)],
                     .title = '',
                     .y.label = '',
                     .x.label = '',
                     .margin = c(3,3,1,1),
                     .background = rgb(1,1,1,0),
                     .theme = theme.loop,
                     .zoom = zoom)
    add.text.line.multiple.colors(names(current.ys), 
                                  text.colors = roi.colors$rostral.caudal[names(current.ys)],
                                  .outer = TRUE)
    
    dev.off()
    
  }#; rm(theme.loop)
  rm(current.sig, current.sig.windows, current.ys, current.xs)
  
  
  ##
  ## Plot: MAX accuracy - all ROIs one plot
  ##
  
  # Get data
  current.ys <- lapply(time.series[[model.loop]][rois], function(.roi){apply(bind_cols(lapply(.roi, function(.classifier){.classifier$accuracy})), 1, max, na.rm = TRUE)})
  current.xs <- time.series[[model.loop]][[1]][[1]]$sample.label
  
  current.sig <- lapply(time.series[[model.loop]][rois], function(.roi){
    apply(bind_cols(lapply(.roi, function(.classifier){
      # Get a "z"score: vals over 1 significant (over 95% CI), under 1 non sig
      .z.ish <- (.classifier$accuracy - .classifier$shuffle.mean) / 
        (.classifier$shuffle.99ci.upper - .classifier$shuffle.mean)
      return(.z.ish)
    })), 1, max, na.rm = TRUE) > 1})
  current.sig.windows <- lapply(current.sig, function(.sig){
    get.significant.windows(.sig,
                            .sample.labels = current.xs,
                            output.class = "data.frame",
                            include.duration = TRUE,
                            .exclude.sig.durations.under.ms = 100,
                            .sampling.rate = 256,
                            .exclude.times.before.ms = -median.rt.times['pn'] + 50)
  })
  
  # Smooth
  current.ys <- lapply(current.ys, smoothing, n.samples.pre = time.convert(25, 'times', 'samples', 256))
  
  # Sort
  roi.order <- names(sort(sapply(current.sig.windows, function(x){min(x$start.time)})))
  current.ys <- current.ys[roi.order]
  current.sig.windows <- current.sig.windows[roi.order]
  
  for(theme.loop in c('white','black')){ ## UNCOMMENT
    # theme.loop = c('white','black')[1] ## UNCOMMENT
    
    # Save dir
    save.dir.all.rois.all.stages <- 
      paste0(save.plot.dir.main,
             'max accuracy - all ROIs/',
             theme.loop,'/')
    dir.create(save.dir.all.rois.all.stages, showWarnings = FALSE, recursive = TRUE)
    
    # # Save
    pdf(paste0(save.dir.all.rois.all.stages, 
               'picture naming cross validation accuracies - maxed - ',
               model.loop,'.pdf'),
        width = 7, height = 5.5)
    # par(oma = c(0,0,1,0))
    plot.time.series(.y.values = current.ys,
                     .x.values = current.xs,
                     .sampling.rate = 256,
                     .x.limits = c(-median.rt.times['pn'], 500),
                     .x.ticks = c(-median.rt.times['pn'], 0, 500),
                     .y.limits = c(.15, .35),
                     .y.ticks = c(.15, .35),
                     .colors = roi.colors$rostral.caudal[names(current.ys)],
                     .sig.windows = current.sig.windows,
                     .sig.color = roi.colors$rostral.caudal[names(current.ys)],
                     .title = '',
                     .y.label = '',
                     .x.label = '',
                     .margin = c(3,3,1,1),
                     .background = rgb(1,1,1,0),
                     .theme = theme.loop,
                     .zoom = zoom)
    # plot.axes(.y.limits = c(0, .4),
    #           .y.ticks = c(0, .4),
    #           show.x.axis = FALSE,
    #           .zoom = zoom)
    add.text.line.multiple.colors(names(current.ys), 
                                  text.colors = roi.colors$rostral.caudal[names(current.ys)],
                                  .outer = TRUE)
    
    dev.off()
    
  }#; rm(theme.loop)
  rm(current.sig, current.sig.windows, current.ys, current.xs)
  
}; rm(model.loop)


###
### Model comparisons, single ROIs
###

model.colors <- c('multinom' = rgb(0, 0, 0),
                  'nnet' = rgb(.5, .5, .5))

##
## Plot: Mean accuracy - one ROI, both models
##

for(roi.loop in rois){
  # roi.loop = rois[1]
  
  # Get data
  current.ys <- lapply(time.series, function(.model){
    apply(bind_cols(lapply(.model[[roi.loop]], function(.classifier){.classifier$accuracy})), 1, mean, na.rm = TRUE)})
  current.xs <- time.series[[1]][[1]][[1]]$sample.label
  current.sig <- lapply(time.series, function(.model){
    apply(bind_cols(lapply(.model[[roi.loop]], function(.classifier){
      # Get a "z"score: vals over 1 significant (over 95% CI), under 1 non sig
      .z.ish <- (.classifier$accuracy - .classifier$shuffle.mean) / 
        (.classifier$shuffle.99ci.upper - .classifier$shuffle.mean)
      return(.z.ish)
    })), 1, max, na.rm = TRUE) > 1}) # get the max -- if any sig
  current.sig.windows <- lapply(current.sig, function(.sig){
    get.significant.windows(.sig,
                            .sample.labels = current.xs,
                            output.class = "data.frame",
                            include.duration = TRUE,
                            .exclude.sig.durations.under.ms = 100,
                            .sampling.rate = 256,
                            .exclude.times.before.ms = -median.rt.times['pn'] + 50)
  })
  
  # # Smooth
  # current.ys <- lapply(current.ys, smoothing, n.samples.pre = time.convert(25, 'times', 'samples', 256))
  
  for(theme.loop in c('black','white')){ ## UNCOMMENT
    # theme.loop = c('white','black')[1] ## UNCOMMENT
    
    # Save dir
    save.dir.all.rois.all.stages <- 
      paste0(save.plot.dir.main,
             'mean accuracy - by ROI - both models/',
             theme.loop,'/')
    dir.create(save.dir.all.rois.all.stages, showWarnings = FALSE, recursive = TRUE)
    
    # # Save
    pdf(paste0(save.dir.all.rois.all.stages, 
               'picture naming cross validation accuracies - averaged - ',
               roi.loop,'.pdf'),
        width = 6, height = 5.5)
    par(oma = c(0,0,1,0))
    plot.time.series(.y.values = current.ys,
                     .x.values = current.xs,
                     .sampling.rate = 256,
                     .x.limits = c(-median.rt.times['pn'], 500),
                     .x.ticks = c(-median.rt.times['pn'], 0, 500),
                     # .horizontal.line.at = 1/6,
                     .y.limits = c(.15, .25),
                     .y.ticks = c(.15, .25),
                     .colors = model.colors[names(current.ys)],
                     .sig.windows = current.sig.windows,
                     .sig.color = model.colors[names(current.ys)],
                     .title = '',
                     .y.label = '',
                     .x.label = '',
                     .margin = c(3,3,1,1),
                     .background = rgb(1,1,1,0),
                     .theme = theme.loop,
                     .zoom = zoom)
    add.text.line.multiple.colors(names(current.ys), 
                                  text.colors = model.colors[names(current.ys)],
                                  .outer = TRUE)
    
    dev.off()
    
  }#; rm(theme.loop)
  rm(current.sig, current.sig.windows, current.ys, current.xs)
  
}#; rm(roi.loop)



##
## Plot: MAX accuracy - one ROI, both models
##

for(roi.loop in rois){
  # roi.loop = rois[1]
  
  # Get data
  current.ys <- lapply(time.series, function(.model){
    apply(bind_cols(lapply(.model[[roi.loop]], function(.classifier){.classifier$accuracy})), 1, max, na.rm = TRUE)})
  current.xs <- time.series[[1]][[1]][[1]]$sample.label
  current.sig <- lapply(time.series, function(.model){
    apply(bind_cols(lapply(.model[[roi.loop]], function(.classifier){
      # Get a "z"score: vals over 1 significant (over 95% CI), under 1 non sig
      .z.ish <- (.classifier$accuracy - .classifier$shuffle.mean) / 
        (.classifier$shuffle.99ci.upper - .classifier$shuffle.mean)
      return(.z.ish)
    })), 1, max, na.rm = TRUE) > 1}) # get the max -- if any sig
  current.sig.windows <- lapply(current.sig, function(.sig){
    get.significant.windows(.sig,
                            .sample.labels = current.xs,
                            output.class = "data.frame",
                            include.duration = TRUE,
                            .exclude.sig.durations.under.ms = 100,
                            .sampling.rate = 256,
                            .exclude.times.before.ms = -median.rt.times['pn'] + 50)
  })
  
  # # Smooth
  # current.ys <- lapply(current.ys, smoothing, n.samples.pre = time.convert(25, 'times', 'samples', 256))
  
  for(theme.loop in c('black','white')){ ## UNCOMMENT
    # theme.loop = c('white','black')[1] ## UNCOMMENT
    
    # Save dir
    save.dir.all.rois.all.stages <- 
      paste0(save.plot.dir.main,
             'max accuracy - by ROI - both models/',
             theme.loop,'/')
    dir.create(save.dir.all.rois.all.stages, showWarnings = FALSE, recursive = TRUE)
    
    # # Save
    pdf(paste0(save.dir.all.rois.all.stages, 
               'picture naming cross validation accuracies - maxed - ',
               roi.loop,'.pdf'),
        width = 6, height = 5.5)
    par(oma = c(0,0,1,0))
    plot.time.series(.y.values = current.ys,
                     .x.values = current.xs,
                     .sampling.rate = 256,
                     .x.limits = c(-median.rt.times['pn'], 500),
                     .x.ticks = c(-median.rt.times['pn'], 0, 500),
                     # .horizontal.line.at = 1/6,
                     .y.limits = c(.15, .25),
                     .y.ticks = c(.15, .25),
                     .colors = model.colors[names(current.ys)],
                     .sig.windows = current.sig.windows,
                     .sig.color = model.colors[names(current.ys)],
                     .title = '',
                     .y.label = '',
                     .x.label = '',
                     .margin = c(3,3,1,1),
                     .background = rgb(1,1,1,0),
                     .theme = theme.loop,
                     .zoom = zoom)
    add.text.line.multiple.colors(names(current.ys), 
                                  text.colors = model.colors[names(current.ys)],
                                  .outer = TRUE)
    
    dev.off()
    
  }#; rm(theme.loop)
  rm(current.sig, current.sig.windows, current.ys, current.xs)
  
}#; rm(roi.loop)



#time.series.by.stage$multinom$IPL$window_200_to_250ms$NY869
max.accs.by.stage <- list()
for(stage.loop in names(time.series.by.stage[[1]][[1]])){
  # stage.loop = names(time.series.by.stage[[1]][[1]])[1]
  
  max.accs.by.stage[[stage.loop]] <- list()
  
  for(model.loop in model.types){
    # model.loop = model.types[1]
    
    max.accs.by.stage[[stage.loop]][[model.loop]] <-
      unlist(sapply(time.series.by.stage[[model.loop]], function(.roi){
        unlist(.roi[[stage.loop]])
      }))
    
  }; rm(model.loop)
  
}; rm(stage.loop)

# Sort by train window time
stages.in.order <- stage.ranges[order(stage.ranges$start.time),]$label
max.accs.by.stage <- max.accs.by.stage[stages.in.order]
names(max.accs.by.stage) <- as.character(stage.ranges[stages.in.order, 'start.time'])

# Save dir
save.model.comparison.barplot.dir <- 
  paste0(save.plot.dir.main,
         'model comparison barplot/')
dir.create(save.model.comparison.barplot.dir, showWarnings = FALSE, recursive = TRUE)

for(theme.loop in c('black','white')){ ## UNCOMMENT
  # theme.loop = c('white','black')[1] ## UNCOMMENT
  
  # # Save
  pdf(paste0(save.model.comparison.barplot.dir, 
             'mean classifier accuracies for neural net vs logistic regression classifiers by training stage - ',
             theme.loop,'.pdf'),
      width = 16, height = 4.8)
  
  # Colors
  multinom.color <- 
    colors[[theme.loop]]$rainbow['pink','hex'] 
    #ifelse(theme.loop == 'white', rgb(0,0,0), rgb(1,1,1))
  nnet.color <- 
    ifelse(theme.loop == 'white', rgb(.7,.7,.7), rgb(.4,.4,.4))
  multinom.error.color <- 
    ifelse(theme.loop == 'white', rgb(.6,.6,.6), rgb(.7,.7,.7))
  nnet.error.color <- 
    ifelse(theme.loop == 'white', rgb(.6,.6,.6), rgb(.7,.7,.7))
    # ifelse(theme.loop == 'white', rgb(.7,.7,.7), rgb(.3,.3,.3))
  dot.colors <- 
    # ifelse(theme.loop == 'white', 'black', 'white')
    # ifelse(theme.loop == 'white', rgb(.65,.65,.65), rgb(.6, .6, .6))
    ifelse(theme.loop == 'white', rgb(.8,.8,.8), rgb(.6, .6, .6))
  
  plot.violin(.values = max.accs.by.stage,
              .show.bars = TRUE,
              .bar.spacing = 0,
              .show.violin = FALSE,
              .show.mean = FALSE,
              .dots.side = "both",
              .show.dots = FALSE,
              .point.size = .75,
              .jitter.width = .4,
              .colors = c(multinom.color, nnet.color),
              .error.bar.colors = c(multinom.error.color, nnet.error.color),
              .dot.colors = dot.colors,
              .dot.alpha = .4,
              # .y.limits = c(.15, .35),
              # .y.ticks = c(.15, .25, .35),
              .y.limits = c(1/6, .27),
              .y.ticks = c(1/6, .25),
              .show.x.axis = FALSE,
              .y.tick.labels = c('0.167','0.25'),
              .theme = theme.loop,
              .margin = c(1.5, 3.5, 1.5, 1.5),
              .background = rgb(1,1,1,0),
              .zoom = 1.6)
  
  dev.off()
  
}#; rm(theme.loop)



# }; rm(band.loop)










# Finish!
message('Script completed successfully. ',Sys.time())


