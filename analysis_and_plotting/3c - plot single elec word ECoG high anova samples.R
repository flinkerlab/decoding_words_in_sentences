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
library('pheatmap')
# library('shape') # for Arrows()
# library('lme4')

# Colors
library('scico')
library('pals')
library('cmocean')
library('shades')

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


### Loop thru beta/high gamma data# For cleanup
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


###
### Load data
###

##
## ECoG (for mean word activity time series plots)
##

# Load data
message("Loading ECoG data (<1 min)... ",Sys.time())
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
                          'end.time' = time.convert(605, "times", "samples", 256))
# Get times
median.rt.times <- time.convert(median.rt.samples, 'samples', 'times', 256)

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

# Subset data to just use.these.elecs
for(patient in names(pn.data)){
  pn.data[[patient]] <- pn.data[[patient]][patient.elecs[[patient]]]
}; rm(patient)
gc()

# Words
words <- c('chicken','dog','dracula','frankenstein','ninja','nurse')
n.words <- length(words)

### Load AOV BFs (for significance)
load(paste0(output.path, 'data/2a - bayesian anovas/anova Bayes Factors.RData'))

# # Limit to just those between 50ms after stim onset and 500ms post production
# rows.to.keep <- time.convert(rownames(word.aov.bfs), 'sample.labels', 'times', 256)
# rows.to.keep <- rows.to.keep[(rows.to.keep >= (-median.rt.times['pn'] + 50)) & (rows.to.keep <= 500)]
# rows.to.keep <- time.convert(rows.to.keep, "times", "sample.labels", 256)
# word.aov.bfs <- word.aov.bfs[rows.to.keep,]

# Get sig windows
elecs.with.bfs.over.3 <- list()
for(elec.loop in colnames(word.aov.bfs)){
  # elec.loop = colnames(word.aov.bfs)[1]
  
  elecs.with.bfs.over.3[[elec.loop]] <- 
    get.significant.windows(
      sig.vals = word.aov.bfs[,elec.loop] >= 3,
      .sample.labels = rownames(word.aov.bfs),
      output.class = 'data.frame',
      include.duration = TRUE,
      .exclude.sig.durations.under.ms = 100,
      .exclude.times.before.ms = -median.rt.times['pn'],
      .exclude.times.after.ms = 540,
      .sampling.rate = 256)
  
}; rm(elec.loop)

### Load ROIs
load(paste0(output.path, 'data/0a - definitions/roi metadata/ROIs and temporal lobe splits.RData'))
elec.roi.df <- stack(roi.elecs[temporal.lobe.splits$rostral.caudal])
colnames(elec.roi.df) <- c('elec','roi')
rownames(elec.roi.df) <- elec.roi.df$elec

###
### Plots
###

### Define colors
## Character colors (for word ECoG activity time series plots)
load(paste0(output.path, 'data/0a - definitions/character colors/character colors.RData'))
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))

### Loop thru NMFS and create data for brain plots
# for(theme.loop in c('white','black)){ ## UNCOMMENT
theme.loop = 'white' ## UNCOMMENT

###
### Plot elec word ECoG time series
###


for(elec.loop in colnames(word.aov.bfs)){
  # elec.loop = colnames(word.aov.bfs)[1]
  
  ### Plot if significant
  if(nrow(elecs.with.bfs.over.3[[elec.loop]]) > 0){
    
    ##
    ## Plot ECoG
    ##
    
    # Patient
    patient <- substr(elec.loop, 1, 5)
    
    # Get time series data
    elec.word.data <- split(pn.data[[patient]][[elec.loop]], pn.info[[patient]]$word)
    elec.word.means <- lapply(elec.word.data, colMeans)
    elec.word.ses <- lapply(elec.word.data, function(x){apply(x, 2, function(y){sd(y)/sqrt(length(y))})})
    
    # Smooth
    half.n.smoothing.samples <- time.convert(50, 'times', 'samples', 256)
    elec.word.means <- lapply(elec.word.means, smoothing, n.samples.pre = half.n.smoothing.samples, na.pad = FALSE)
    elec.word.ses <- lapply(elec.word.ses, smoothing, n.samples.pre = half.n.smoothing.samples, na.pad = FALSE)
    
    # Save directory
    save.word.time.series.dir <- 
      paste0(output.path,
             'figures/3c - plot single elec word time series for high anova elecs/',
             theme.loop,
             '/individual plots for publication - word time series',
             'BF over 3 for 100ms/',
             elec.roi.df[elec.loop,'roi'],'/')
    dir.create(save.word.time.series.dir, showWarnings = FALSE, recursive = TRUE)
    
    # Plot time series
    pdf(paste0(save.word.time.series.dir, elec.loop, ' - ',elec.info[elec.loop,'region_clinical'],' - ECoG - max 4 - h5.pdf'),
        width = 7.5, height = 5)
    plot.time.series(.y.values = elec.word.means,
                     .x.values = lapply(elec.word.means, names),
                     .sampling.rate = 256,
                     .error.bars = elec.word.ses,
                     .colors = character.colors[names(elec.word.means), 'hex'],
                     .x.limits = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256) -100, 500),
                     .x.ticks = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256), 0, 500),
                     .sig.windows = elecs.with.bfs.over.3[[elec.loop]],
                     .sig.color = colors[[theme.loop]]$rainbow_bright['pink','hex'],
                     .y.limits.min.at.least = -.4,
                     .y.limits.max.at.least = 4,
                     .y.ticks = c(0, 4),
                     .margin = c(3.5,3.5,1.5,3.5),
                     .theme = theme.loop,
                     .horizontal.line.at = 0,
                     .y.label = '',
                     .x.label = '',
                     .background = rgb(1,1,1,0),
                     .zoom = 1.8)
    
    dev.off()
    
    
    ##
    ## Plot ANOVA BF
    ##
    
    
    # Plot time series
    pdf(paste0(save.word.time.series.dir, elec.loop, ' - ',elec.info[elec.loop,'region_clinical'],' - ANoVA BF - max 25 - h4.5.pdf'),
        width = 7.5, height = 4.5)
    plot.time.series(.y.values = smoothing(log10(word.aov.bfs[,elec.loop]), n.samples.pre = half.n.smoothing.samples),
                     .x.values = rownames(word.aov.bfs),
                     .sampling.rate = 256,
                     .colors = colors[[theme.loop]]$rainbow_bright['pink','hex'],
                     .x.limits = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256) -100, 500),
                     .x.ticks = c(-time.convert(median.rt.samples['pn'], 'samples', 'times', 256), 0, 500),
                     .sig.windows = elecs.with.bfs.over.3[[elec.loop]],
                     .sig.color = colors[[theme.loop]]$rainbow_bright['pink','hex'],
                     .y.limits.min.at.least = -4,
                     .y.limits.max.at.least = 25,
                     .y.ticks = c(0, 25),
                     .margin = c(3.5,3.5,1.5,3.5),
                     .theme = theme.loop,
                     .horizontal.line.at = 0,
                     .y.label = '',
                     .x.label = '',
                     .background = rgb(1,1,1,0),
                     .zoom = 1.8)
    
    dev.off()
    
    # Clean up
    rm(elec.word.data, elec.word.means, elec.word.ses)
    
  } # if any sig windows
} #; rm(elec.loop)


# }#; rm(theme.loop)

# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())







