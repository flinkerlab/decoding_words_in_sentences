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
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))
# Elementwise matrix apply
source(paste0(path,'/analysis/R/functions/elementwise_matrix_apply.R'))
# Squish Bayes Factors
source(paste0(path,'/analysis/R/functions/squish_bayes_factors.R'))
# Weighted mean and standard error functions
source(paste0(path,'/analysis/R/functions/weighted_summary_stats.R'))
# Weighted mean and standard error functions
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

### Load ROIs and stages
load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/stages and rois/stages and rois.RData'))

### Load RTs
load(paste0(path, 'analysis/R/downsample data/warped data/high_gamma/median.rt.samples 256 Hz.RData'))
median.rt.times <- time.convert(median.rt.samples, "samples", "times", 256)

### Load colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))
load(paste0(output.path, 'data/0a - definitions/character colors/character colors.RData'))
load(paste0(output.path, 'data/0a - definitions/roi metadata/ROIs and temporal lobe splits.RData'))
load(paste0(output.path, 'data/0a - definitions/noun colors/noun colors.RData'))

### Load anova BFs
load(paste0(output.path,'data/2a - bayesian anovas/anova Bayes Factors.RData'))


### Define time windows
window.bin.sizes <- c(250, 150, 100)
onset.time <- -time.convert(median.rt.samples['pn'], 'samples', 'times', 256)
windows.list <- list()
for(bin.size.loop in window.bin.sizes){
  # bin.size.loop = window.bin.sizes[1]
  
  # Label
  bin.size.label <- paste0('bin_size_',bin.size.loop,'ms')
  
  # Get start and end times
  starts <- rev(seq(0, onset.time, by = -bin.size.loop))
  starts[1] <- onset.time # make the first window longer to include everything from stim onset
  ends <- rev(seq(bin.size.loop, starts[2], by = -bin.size.loop))
  
  # Put in dataframe
  windows.list[[bin.size.label]] <-
    data.frame('start.time' = starts,
               'end.time' = ends,
               'start.sample' = time.convert(starts, 'times', 'samples', 256),
               'end.sample' = time.convert(ends, 'times', 'samples', 256))
  
  # Remove overlap in adjacent windows
  windows.list[[bin.size.label]]$end.sample <-
    windows.list[[bin.size.label]]$end.sample - 1
  
  # Add label
  windows.list[[bin.size.label]]$window.label <-
    with(windows.list[[bin.size.label]],
         gsub("-","neg",paste0('bin_',start.time,'ms_to_',end.time), fixed = TRUE))
  
}; rm(bin.size.loop)



### Loop thru bin.sizes and get means and maxes for each elec
for(bin.size.loop in window.bin.sizes){
  # bin.size.loop = window.bin.sizes[1]
  
  # Label
  bin.size.label <- paste0('bin_size_',bin.size.loop,'ms')
  
  # Storage
  elec.means <- list()
  elec.maxes <- list()
  
  # Loop thru windows
  for(window.loop in 1:nrow(windows.list[[bin.size.label]])){
    # window.loop = c(1:nrow(windows.list[[bin.size.label]]))[1]
    
    # Label
    current.window <- windows.list[[bin.size.label]][window.loop,'window.label']
    
    current.sample.labels <- 
      time.convert(windows.list[[bin.size.label]][window.loop,'start.sample']:
                     windows.list[[bin.size.label]][window.loop,'end.sample'],
                   'samples','sample.labels')
    
    # Store summary stats
    elec.means[[current.window]] <-
      apply(word.aov.bfs[current.sample.labels,], 2, function(x){
        # log(mean(x), base = 10)})
        mean(log(x, base = 10))})
    elec.maxes[[current.window]] <-
      apply(word.aov.bfs[current.sample.labels,], 2, function(x){
        # log(max(x), base = 10)})
        max(log(x, base = 10))})
    
  }; rm(window.loop)
  
  # Make into dataframes
  elec.means <- data.frame(bind_rows(elec.means, .id = 'window'))
  elec.maxes <- data.frame(bind_rows(elec.maxes, .id = 'window'))
  
  # Rownames
  rownames(elec.means) <- elec.means$window
  rownames(elec.maxes) <- elec.maxes$window
  elec.means$window <- elec.maxes$window <- NULL
  
  # Rearrange
  elec.means <- data.frame(t(elec.means))
  elec.maxes <- data.frame(t(elec.maxes))
  
  # Add elec column
  elec.means$elec <- rownames(elec.means)
  elec.maxes$elec <- rownames(elec.maxes)
  
  # Add region column
  elec.means$region_clinical <- elec.info[elec.means$elec,'region_clinical']
  elec.maxes$region_clinical <- elec.info[elec.maxes$elec,'region_clinical']
  
  # Add MNI
  elec.means$MNI_x <- elec.info[elec.means$elec,'MNI_x']
  elec.means$MNI_y <- elec.info[elec.means$elec,'MNI_y']
  elec.means$MNI_z <- elec.info[elec.means$elec,'MNI_z']
  elec.maxes$MNI_x <- elec.info[elec.maxes$elec,'MNI_x']
  elec.maxes$MNI_y <- elec.info[elec.maxes$elec,'MNI_y']
  elec.maxes$MNI_z <- elec.info[elec.maxes$elec,'MNI_z']
  
  ### Save
  save.data.dir <- paste0(output.path, 'data/2c - anova BF in windows for brain plots/',
                          bin.size.label,'/')
  dir.create(save.data.dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(elec.means,
            paste0(save.data.dir, 'elec_mean_log_BF_in_bins.csv'),
            row.names = FALSE, quote = FALSE)
  write.csv(elec.maxes,
            paste0(save.data.dir, 'elec_max_log_BF_in_bins.csv'),
            row.names = FALSE, quote = FALSE)
    
  
  ###  
  ### Get ROI means in bins for barplots and do barplots
  ###
  
  ## Storage of summary stats
  # Mean BFs in bin
  roi.mean.means <- list()
  roi.sd.means <- list()
  roi.se.means <- list()
  
  # Max BFs in bin
  roi.mean.maxes <- list()
  roi.sd.maxes <- list()
  roi.se.maxes <- list()
  
  # Loop thru ROIs and fill in
  for(roi.loop in rois){
    # roi.loop = rois[1]
    
    current.elecs <- roi.elecs[[roi.loop]]
    current.elecs <- current.elecs[current.elecs %in% rownames(elec.means)]
    
    # Stats on mean BFs
    current.data <- elec.means[current.elecs,windows.list[[bin.size.label]]$window.label]
    roi.mean.means[[roi.loop]] <- apply(current.data, 2, mean)
    roi.sd.means[[roi.loop]] <- apply(current.data, 2, sd)
    roi.se.means[[roi.loop]] <- apply(current.data, 2, function(x){
      sd(x) / sqrt(length(nrow(x)))
    })
    
    # Stats on max BFs
    current.data <- elec.maxes[current.elecs,windows.list[[bin.size.label]]$window.label]
    roi.mean.maxes[[roi.loop]] <- apply(current.data, 2, mean)
    roi.sd.maxes[[roi.loop]] <- apply(current.data, 2, sd)
    roi.se.maxes[[roi.loop]] <- apply(current.data, 2, function(x){
      sd(x) / sqrt(length(nrow(x)))
    })
    
  }; rm(roi.loop)
  
  ## Reorganize data
  roi.mean.means <- as.matrix(data.frame(bind_rows(roi.mean.means, .id = 'roi'), row.names = 'roi'))
  
  
  ###
  ### Barplots!
  ###
  
  for(tl.split.loop in names(temporal.lobe.splits)){
    # tl.split.loop = names(temporal.lobe.splits)[2]
    
    for(theme.loop in c('black','white')){
      # theme.loop = 'white'
    
    current.rois <- temporal.lobe.splits[[tl.split.loop]]
    
    for(data.type.loop in c('window.means','window.maxes','window.count.bf.over.3')){
      # data.type.loop = c('window.means','window.maxes','window.count.bf.over.3')[1]
  
      if(data.type.loop == 'window.means'){
        current.data <- elec.means[,windows.list[[bin.size.label]]$window.label]
        }
      if(data.type.loop == 'window.maxes'){
        current.data <- elec.maxes[,windows.list[[bin.size.label]]$window.label]
      }
      if(data.type.loop == 'window.count.bf.over.3'){
        current.data <- elec.maxes[,windows.list[[bin.size.label]]$window.label]
        for(col.loop in 1:ncol(current.data)){
          current.data[,col.loop] <- as.numeric(current.data[,col.loop] > log(3, base = 10))
        }; rm(col.loop)
      }
      
    ##
    ## Plot grid
    ##
    
    save.barplot.dir.grid <- paste0(output.path, 
                             'figures/2c - anova BF in windows for brain plots/barplots/',
                             'grids/',
                             tl.split.loop,'/',
                             data.type.loop,'/',
                             theme.loop,'/')
  dir.create(save.barplot.dir.grid, showWarnings = FALSE, recursive = TRUE)
  
  # pdf(paste0(save.barplot.dir.grid, 
  #            window.loop,' - mean and SE by ',bin.size.loop,'ms window.pdf'),
  #     width = 2 * length(current.rois), height = 4)
    
  # dev.off()
    
  
    ##
    ## Plot individual windows
    ##
    
    for(window.loop in windows.list[[bin.size.label]]$window.label){
      # window.loop = windows.list[[bin.size.label]]$window.label[1]
      
      save.barplot.dir.individual <- paste0(output.path, 
                                      'figures/2c - anova BF in windows for brain plots/barplots/',
                                      'individual plots for publication/',
                                      tl.split.loop,'/',
                                      data.type.loop,'/',
                                      theme.loop,'/',
                                      bin.size.label,'/')
      dir.create(save.barplot.dir.individual, showWarnings = FALSE, recursive = TRUE)
      
      pdf(paste0(save.barplot.dir.individual, 
                 window.loop,' - mean and SE by ',bin.size.loop,'ms window.pdf'),
          width = 1 * length(current.rois), height = 5)
      par(oma = c(0,0,0,0))
      
      ### Plot
      violin.data <- list()
      for(roi.loop in current.rois){
        # roi.loop = rois[1]
        current.elecs <- roi.elecs[[roi.loop]]
        current.elecs <- current.elecs[current.elecs %in% rownames(current.data)]
        violin.data[[roi.loop]] <- current.data[current.elecs,window.loop]
      }; rm(roi.loop)
      current.y.limits <- list(
        'window.means' = c(-2, 8),
        'window.maxes' = c(-2, 8),
        'window.count.bf.over.3' = c(0, .6)
      )[[data.type.loop]]
      current.y.ticks <- list(
        'window.means' = c(-2, 0, 8),
        'window.maxes' = c(-2, 0, 8),
        'window.count.bf.over.3' = c(0, .5)
      )[[data.type.loop]]
      
      # Tick labels - make percents for the proprotion loop
      current.y.tick.labels <- current.y.ticks * ifelse(data.type.loop == 'window.count.bf.over.3', 100, 1)
      if(window.loop != windows.list[[bin.size.label]]$window.label[1]){
        current.y.tick.labels <- rep('', times = length(current.y.tick.labels))
      }
      
        # Colors
      if(data.type.loop == 'window.count.bf.over.3'){
        # current.border.colors <- 
        #   current.error.bar.colors <-
        #   cubicl(length(current.rois))
        current.border.colors <- NA
        current.error.bar.colors <- 'black'
        current.bar.colors <- cubicl(length(current.rois))
      }else{
        current.bar.colors <- 'white'
        current.border.colors <- 
          current.error.bar.colors <-
          'black'
      }
      
      par(mfrow = c(1,1))
      set.seed(18)
      plot.violin(.values = violin.data,
                  .colors = cubicl(length(current.rois)),
                  # .colors = rgb(.8, .8, .8),
                  # .horizontal.line.at = 0,
                  .show.bars = TRUE,
                  .show.dots = ifelse(data.type.loop == 'window.count.bf.over.3', FALSE, TRUE),
                  .bar.colors = current.bar.colors,
                  .bar.border.colors = current.border.colors,
                  .bar.spacing = .15,
                  .error.bar.type = 'se',
                  .error.bar.caps = .25,
                  .error.bar.colors = current.error.bar.colors,
                  .y.limits = current.y.limits,
                  .y.ticks = current.y.ticks,
                  .y.tick.labels = current.y.tick.labels,
                  .show.violin = FALSE,
                  .dots.side = 'both',
                  .show.mean = FALSE,
                  .scale.jitter.by.density = TRUE,
                  .max.density.y = .2,
                  .theme = theme.loop,
                  .background = rgb(1,1,1,0),
                  .x.tick.labels.vertical = TRUE,
                  .margin = c(5,4,1,1),
                  .zoom = 1.8)
      # abline(h = 0, lwd = 1.5 * 1.2)
      
      dev.off()
      
    }#; rm(window.loop)

  }#; rm(data.type.loop)
    
    }#; rm(theme.loop)
  
  }#; rm(tl.split.loop)
  
  
}; rm(bin.size.loop)




###
### Plots
###



# } # band.loop








# Finish!
message('Script completed successfully. ',Sys.time())







