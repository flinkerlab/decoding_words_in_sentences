### Plot stages (from multivariate change point detection) per cluster
### July 2024
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
library('grid')
library('cmocean') # for cmocean() colors
library('effsize') # for hedges.g()


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
# Touch up heatmap
source(paste0(path,'/analysis/R/functions/touchup_heatmap.R'))

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
message("Loading ECoG data (<1 min)... ",Sys.time())
load(paste0(path,
            'analysis/R/downsample data/warped data/',band.loop,'/warped ',band.loop,' data at 256 Hz.RData')) # loads "warped.data", "median.rt.samples", "patient.trial.info", and "sampling.rate"
message('...done! ', Sys.time())

### Clean up and get metadata
pn.data <- lapply(warped.data, function(x){x[['pn']]})
rm(warped.data)
pn.info <- lapply(patient.trial.info, function(x){x[['pn']]})
rm(patient.trial.info)

# Words
words <- c('chicken','dog','dracula','frankenstein','ninja','nurse')
n.words <- length(words)


### Load AOV BFs for weighting PN data in windows (word.aov.bfs)
load(paste0(path,'analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/data/2a - bayesian anovas/anova Bayes Factors.RData'))

### Load stages (roi.elecs, rois, stage.ranges, stage.sample.labels)
load(paste0(path,'analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/data/9a - train classifiers - get variable importances - 50ms/stages and rois/stages and rois.RData'))

###
### Get plot data -- Hedges g per elec per stage
###

word.hedges.gs <- list()
for(word.loop in words){
  # word.loop = words[1]
  
  # Storage
  word.hedges.gs[[word.loop]] <-
    data.frame(matrix(nrow = length(use.these.elecs),
                      ncol = length(stage.sample.labels),
                      dimnames = list(use.these.elecs,
                                      names(stage.sample.labels))))
  
  for(elec.loop in use.these.elecs){
    # elec.loop = use.these.elecs[1]
    
    patient <- substr(elec.loop, 1, 5)
    
    for(stage.loop in names(stage.sample.labels)){
      # stage.loop = names(stage.sample.labels)[1]
      
      current.word.data <- 
        pn.data[[patient]][[elec.loop]][pn.info[[patient]]$word == word.loop, 
                                        stage.sample.labels[[stage.loop]]]
      current.other.word.data <- 
        pn.data[[patient]][[elec.loop]][pn.info[[patient]]$word %in% words[words != word.loop], 
                                        stage.sample.labels[[stage.loop]]]
      current.weights <-  word.aov.bfs[stage.sample.labels[[stage.loop]], elec.loop]
      
      # Average across time
      current.word.data <- 
        apply(current.word.data, 1, function(x){
          weighted.average(x, w = current.weights)})
      
      current.other.word.data <- 
        apply(current.other.word.data, 1, function(x){
          weighted.average(x, w = current.weights)})
      
      # Get normalized distance between current word and others
      word.hedges.gs[[word.loop]][elec.loop, stage.loop] <- 
        cohen.d(
          d = c(current.word.data, current.other.word.data), # data
          f = as.factor(c(rep(1, times = length(current.word.data)),
                          rep(0, times = length(current.other.word.data)))), # labels
          hedges.correction = TRUE # correct for imbalanced sample sizes
        )$estimate
      
    }#; rm(stage.loop)      
    
  }#; rm(elec.loop)
  
  # Get max hedges g over time
  word.hedges.gs[[word.loop]] <-
    data.frame(hedges.g = apply(word.hedges.gs[[word.loop]], 1, max))
  colnames(word.hedges.gs[[word.loop]]) <- word.loop
  
}#; rm(word.loop)

# Combine into dataframe
word.hedges.gs <- 
  data.frame(do.call(cbind, word.hedges.gs))

# Define threshold for high hedges g
high.g <- .3

# Get number of words each elec has a high g for
n.words.per.elec <- apply(word.hedges.gs > high.g, 1, sum)
one.word.elecs <- names(n.words.per.elec[n.words.per.elec == 1])

# Get elecs that are only selective for one word
one.word.elecs <- word.hedges.gs[one.word.elecs, ]
one.word.elecs <- data.frame(apply(one.word.elecs, 2, function(x){ifelse(x < high.g, NA, x)}))

# Split by word
selective.elec.gs <- list()
for(word.loop in words){
  # word.loop = words[1]
  
  selective.elec.gs[[word.loop]] <- one.word.elecs[, word.loop, drop = FALSE]
  colnames(selective.elec.gs[[word.loop]]) <- 'hedges.g'
  
  # Add elec details
  selective.elec.gs[[word.loop]]$elec <- rownames(selective.elec.gs[[word.loop]])
  selective.elec.gs[[word.loop]] <- cbind(selective.elec.gs[[word.loop]],
                                          elec.info[selective.elec.gs[[word.loop]]$elec,
                                                    c('MNI_x','MNI_y','MNI_z','region_clinical_for_plots')])
  colnames(selective.elec.gs[[word.loop]])[
    colnames(selective.elec.gs[[word.loop]]) == 'region_clinical_for_plots'] <- 'region_clinical'
  
  # Get rid of NAs
  selective.elec.gs[[word.loop]] <- 
    droplevels(subset(selective.elec.gs[[word.loop]], ! is.na(hedges.g)))
  

  ### Save
save.data.dir <- paste0(output.path,'/data/7c - word-specific brain plots/hedges g values/')
dir.create(save.data.dir, showWarnings = FALSE, recursive = TRUE)
write.csv(selective.elec.gs[[word.loop]],
          file = paste0(save.data.dir, 'elecs selective for just ',word.loop,' - hedges g values.csv'),
          row.names = FALSE, quote = FALSE)

}; rm(word.loop)



# }; rm(band.loop)















# Finish!
message('Script completed successfully. ',Sys.time())







