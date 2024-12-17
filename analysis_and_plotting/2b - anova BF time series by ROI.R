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

### Load ECoG means ('median.rt.samples','region.means','region.ses','sample.labels','times')
load(paste0(output.path,'data/0b - plot ECoG - ROIs - warped data/roi stats/ECoG ROI means.RData'))
# Convert from 512Hz to 256
median.rt.times <- time.convert(median.rt.samples, 'samples', 'times', 512)
median.rt.samples <- time.convert(median.rt.samples, 'samples', 'samples', input.sampling.rate = 512, output.sampling.rate = 256)

### Load anova BFs
load(paste0(output.path,'data/2a - bayesian anovas/anova Bayes Factors.RData'))

### Define ROIs
load(paste0(output.path, 'data/0a - definitions/roi metadata/ROIs and temporal lobe splits.RData'))
rois <- names(roi.categories)
roi.elecs <- lapply(roi.elecs, function(x){
  x[x %in% colnames(word.aov.bfs)]
})

### Load colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))

### Loop thru ROIs and store means and maxes of AOV BFs
roi.means <- list()
roi.sds <- list()
roi.ses <- list()
roi.proportion.over.bf3 <- list()
for(roi.loop in rois){
  # roi.loop = rois[1]
  
  # Get anova BFs for elecs in this ROI
  current.bfs <- word.aov.bfs[, roi.elecs[[roi.loop]]]
  
  # What proportion over 3 over time?
  roi.proportion.over.bf3[[roi.loop]] <- 
    setNames(apply(apply(current.bfs, 2, function(x){as.numeric(x >= 3)}), 1, mean),
             rownames(current.bfs))
  
  # Limit to just those with any BF > 3
  current.elecs.high.bf <- names(which(apply(current.bfs, 2, function(x){any(x > 3)})))
  current.bfs <- current.bfs[, current.elecs.high.bf]
  
  # Convert to linear scale for averaging etc.
  # current.bfs <- log(current.bfs, base = 10)
  
  # Store summary stats if more than one elec
  if(ncol(current.bfs) > 1){
    
    # Get summary stats
    roi.means[[roi.loop]] <- apply(current.bfs, 1, mean)
    roi.sds[[roi.loop]] <- apply(current.bfs, 1, sd)
    roi.ses[[roi.loop]] <- roi.sds[[roi.loop]] / sqrt(ncol(current.bfs))
    
    message('n elecs in ',roi.loop,': ',ncol(current.bfs))
    
    # # Unlog
    # roi.means[[roi.loop]] <- 10^roi.means[[roi.loop]]
    # roi.sds[[roi.loop]] <- 10^roi.sds[[roi.loop]]
    # roi.ses[[roi.loop]] <- 10^roi.ses[[roi.loop]]

  } # if(ncol(current.bfs) > 1)
  
}; rm(roi.loop)


###
### Plots
###

##
## Regions 
##

xs.sample.labels <- rownames(word.aov.bfs)
xs.samples <- time.convert(xs.sample.labels, "sample.labels", "samples", 256)
x.stim.locked <- time.convert(xs.samples + median.rt.samples['pn'], "samples", "times", 256)


for(theme.loop in c('black','white')){
  # theme.loop = 'white'
  
  ##
  ## Plot mean BF per region
  ##
  
  save.mean.fig.dir <- paste0(output.path, 'figures/2b - plot anova BF time series by ROI/mean BF/',theme.loop,'/')
  dir.create(save.mean.fig.dir, showWarnings = FALSE, recursive = TRUE)
  
  save.prop.fig.dir <- paste0(output.path, 'figures/2b - plot anova BF time series by ROI/proprotion BF over 3/',theme.loop,'/')
  dir.create(save.prop.fig.dir, showWarnings = FALSE, recursive = TRUE)
  
  save.mean.with.ecog.fig.dir <- paste0(output.path, 'figures/2b - plot anova BF time series by ROI/individual regions - mean BF and ECoG/',theme.loop,'/')
  dir.create(save.mean.with.ecog.fig.dir, showWarnings = FALSE, recursive = TRUE)
  
  corrplot.dir <- paste0(output.path, 'figures/2b - plot anova BF time series by ROI/ROI correlations between mean BF and ECoG/',theme.loop,'/')
  dir.create(corrplot.dir, showWarnings = FALSE, recursive = TRUE)
  
  for(tl.split.loop in names(temporal.lobe.splits)){
    # tl.split.loop = names(temporal.lobe.splits)[1]
    
    pdf(paste0(save.mean.fig.dir, 'anova BFs by ROI - temporal lobe split along ',tl.split.loop,' axis - mean BF.pdf'),
        width = 9, height = 4.5)
    # par(oma = c(0,0,2,0))
    
    current.rois <- temporal.lobe.splits[[tl.split.loop]]
    current.ys <- roi.means[current.rois]
    current.ys <- lapply(current.ys, log10)
    current.ys <- lapply(current.ys,
                         smoothing,
                         n.samples.pre = time.convert(50, 'times', 'samples', 256))
    current.xs <- lapply(current.ys, function(x){time.convert(names(x), 'sample.labels', 'times', 256)})
    keep.indices <- lapply(current.xs, function(x){which(x > (-median.rt.times['pn'] - 100))})
    for(i in 1:length(current.ys)){
      current.ys[[i]] <- current.ys[[i]][keep.indices[[i]]]
      current.xs[[i]] <- current.xs[[i]][keep.indices[[i]]]
    }; rm(i)
    
    plot.time.series(.y.values = current.ys,
                     .x.values = current.xs,
                     .sampling.rate = 256,
                     .colors = cubicl(length(current.rois)),
                     .y.ticks = c(0, 40),
                     .y.limits = c(0, 40),
                     .y.label = '',
                     .x.label = '',
                     # .x.limits = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256) -100, 500),
                     # .x.limits = c(-median.rt.times['sp'], 1000),
                     .x.limits = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256) -100, 500),
                     # .x.limits = c(-1000, 500),
                     .x.ticks = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256), 0, 1000),
                     .theme = theme.loop,
                     .background = rgb(1,1,1,0),
                     .zoom = 1.8,
                     .margin = c(3.5,3.5,1.5,1.5))
    add.text.line.multiple.colors(text.segments = current.rois,
                                  text.colors = cubicl(length(current.rois)),
                                  .outer = TRUE)
    
    dev.off()
  
    
    ##
    ## Plot proportion
    ##
    
      pdf(paste0(save.prop.fig.dir, 'anova BFs by ROI - temporal lobe split along ',tl.split.loop,' axis - proportion elecs over BF=3.pdf'),
          width = 7, height = 4.5)
      # par(oma = c(0,0,2,0))
      
      current.rois <- temporal.lobe.splits[[tl.split.loop]]
      current.ys <- roi.proportion.over.bf3[current.rois]
      current.ys <- lapply(current.ys,
                           smoothing,
                           n.samples.pre = time.convert(50, 'times', 'samples', 256))
      
      plot.time.series(.y.values = current.ys,
                       .x.values = lapply(current.ys, names),
                       .sampling.rate = 256,
                       .colors = cubicl(length(current.rois)),
                       .y.ticks = c(0, .3),
                       .y.limits = c(0, .3),
                       .y.label = '',
                       .x.label = '',
                       # .x.limits = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256) -100, 500),
                       .x.limits = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256) -100, 500),
                       .x.ticks = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256), 0, 500),
                       # .vertical.line.at = c(time.convert(median.rt.samples['pn'], "samples", "times", 256)),
                       # show.t0 = FALSE,
                       .theme = theme.loop,
                       .background = rgb(1,1,1,0),
                       .zoom = 1.8,
                       .margin = c(3.5,3.5,1.5,1.5))
      add.text.line.multiple.colors(text.segments = current.rois,
                                    text.colors = cubicl(length(current.rois)),
                                    .outer = TRUE)
      
      dev.off()
      
      
      ##
      ## Plot individual regions: mean BF and mean ECoG
      ##
      
      # Set up storage for correlations
      ys.for.corr <- list()
      
      # Loop thru ROIs and plot time series
      for(roi.loop in temporal.lobe.splits[[tl.split.loop]]){
        # roi.loop = rois[1]
        
        # Get data
        half.n.smoothing.samples.512 <- 
          time.convert(50, 'times', 'samples', 512)
        half.n.smoothing.samples.256 <- 
          time.convert(half.n.smoothing.samples.512, 'samples', 'samples', 
                       input.sampling.rate = 512, output.sampling.rate = 256)
        current.ys <- list('bf' = smoothing(log10(roi.means[[roi.loop]]), n.samples.pre = half.n.smoothing.samples.256, na.pad = FALSE),
                           'ecog' = smoothing(region.means[['pn']][[roi.loop]], n.samples.pre = half.n.smoothing.samples.512, na.pad = FALSE))
        current.xs <- list('bf' = time.convert(names(current.ys[['bf']]), 'sample.labels', 'times', 256),
                           'ecog' = time.convert(names(current.ys[['ecog']]), 'sample.labels', 'times', 512))
        current.error.bars <- list('bf' = rep(0, times = length(current.ys[['bf']])),
                                   'ecog' = region.ses[['pn']][[roi.loop]][(half.n.smoothing.samples.512 + 1):(length(region.ses[['pn']][[roi.loop]]) - half.n.smoothing.samples.512)])
        
        # Make starting time sample the same
        earliest.time <- max(sapply(current.xs, min))
        latest.time <- 540
        keep.samples <- lapply(current.xs, function(x){which((x >= earliest.time) & (x < latest.time))})
        for(i in 1:length(current.ys)){
          # i = 2
          current.ys[[i]] <- current.ys[[i]][keep.samples[[i]]]
          current.xs[[i]] <- current.xs[[i]][keep.samples[[i]]]
          current.error.bars[[i]] <- current.error.bars[[i]][keep.samples[[i]]]
        }; rm(i)
        
        
        # Plot
        pdf(paste0(save.mean.with.ecog.fig.dir, 'mean anova BF and ECoG - ',roi.loop,'.pdf'),
            width = 6, height = 4)
        
        plot.time.series(.y.values = current.ys,
                         .x.values = current.xs,
                         .which.y.axis = c(2,1),
                         .sampling.rate = 256,
                         .colors = c(colors[[theme.loop]]$rainbow_bright['pink','hex'],
                                   # roi.colors[[tl.split.loop]][roi.loop]),
                           ifelse(theme.loop == 'white', rgb(.4, .4, .4), rgb(.6, .6, .6))),
                         .error.bars = current.error.bars,
                         .y.limits = c(0, 1),
                         .y.ticks = c(0, 1),
                         # .y.tick.labels = c(ifelse(roi.loop == 'IFG', 0, ''), ifelse(roi.loop == 'IFG', 1, '')),
                         # show.y2.axis = roi.loop == 'ATL', # only show for last in row
                         # show.y.axis = roi.loop == 'IFG', # only show for first in row
                         # .y.axis.color = roi.colors[[tl.split.loop]][roi.loop],#rgb(.5,.5,.5),
                         .y2.limits = c(0, 40),
                         .y2.ticks = c(0, 40),
                         # .y2.tick.labels = c(ifelse(roi.loop == 'ATL', 0, ''), ifelse(roi.loop == 'ATL', 40, '')),
                         .y2.axis.color = colors[[theme.loop]]$rainbow_bright['pink','hex'],
                         .y.label = '',
                         .x.label = '',
                         .x.limits = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256) -100, 500),
                         .x.ticks = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256), 0, 500),
                         # .vertical.line.at = c(time.convert(median.rt.samples['pn'], "samples", "times", 256)),
                         # show.t0 = FALSE,
                         .theme = theme.loop,
                         .background = rgb(1,1,1,0),
                         .zoom = 1.8,
                         .margin = c(3.5,3.5,1.5,3.5))
        
        dev.off()
        
      
        ### Store time series in same sampling rate for correlations
        current.ys$ecog <- downsample.trial(current.ys$ecog,
                                            input.rate = 512, 
                                            output.rate = 256)
        ys.for.corr[[roi.loop]] <- current.ys
        
      }#; rm(roi.loop)
      
      
      ##
      ## Plot correlation matrices
      ##
      
      # Storage matrix
      anova.ecog.cor <- matrix(nrow = length(temporal.lobe.splits[[tl.split.loop]]),
                               ncol = length(temporal.lobe.splits[[tl.split.loop]]),
                               dimnames = list(temporal.lobe.splits[[tl.split.loop]],
                                               temporal.lobe.splits[[tl.split.loop]]))
      
      for(ecog.roi.loop in temporal.lobe.splits[[tl.split.loop]]){
        # ecog.roi.loop = temporal.lobe.splits[[tl.split.loop]][1]
        for(bf.roi.loop in temporal.lobe.splits[[tl.split.loop]]){
          # bf.roi.loop = temporal.lobe.splits[[tl.split.loop]][2]
          
          anova.ecog.cor[ecog.roi.loop, bf.roi.loop] <-
            cor(ys.for.corr[[ecog.roi.loop]]$ecog,
                ys.for.corr[[bf.roi.loop]]$bf)
          
        }#; rm(bf.roi.loop)
      }#; rm(ecog.roi.loop)
      
      
      pdf(paste0(corrplot.dir,'x=ecog correlation with y=anova bf - temporal lobe split=',tl.split.loop,'.pdf'),
          height = 7, width = 8)
      
      corrplot(anova.ecog.cor, 
               method = 'color',
               is.corr = FALSE,
               bg = rgb(1,1,1,0),
               # tl.col = roi.colors[[tl.split.loop]][rownames(anova.ecog.cor)],
               tl.col = 'black',
               tl.cex = 1.6 * 1.8,
               tl.offset = .5,
               cl.cex = 1.4 * 1.8,
               cl.length = 2,
               cl.ratio = .2,
               col.lim = c(0, 1))
      dev.off()
  
  }#; rm(tl.split.loop)
  
}#; rm(theme.loop)


# } # band.loop








# Finish!
message('Script completed successfully. ',Sys.time())







