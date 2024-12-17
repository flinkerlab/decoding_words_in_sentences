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
# Reverse hierarchy of list of dataframes
source(paste0(path,'/analysis/R/functions/reverse.data.frame.list.hierarchy.R'))


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
load(paste0(output.path, 'data/0a - definitions/roi metadata/ROIs and temporal lobe splits.RData'))
load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/stages and rois/stages and rois.RData'))


### Load RTs
load(paste0(path, 'analysis/R/downsample data/warped data/high_gamma/median.rt.samples 256 Hz.RData'))
median.rt.times <- time.convert(median.rt.samples, "samples", "times", 256)
rm(median.rt.samples) # risky to keep since dealing with 256 and 512 Hz datasets

### Load colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))
load(paste0(output.path, 'data/0a - definitions/character colors/character colors.RData'))
load(paste0(output.path, 'data/0a - definitions/noun colors/noun colors.RData'))


### Loop thru models
model.types <- c('multinom','nnet')
for(model.loop in model.types){ ## UNCOMMENT
# model.loop = model.types[1] ## UNCOMMENT

### Loop thru rois
for(roi.loop in rev(rois)){ ## UNCOMMENT
  # roi.loop = rois[1] ## UNCOMMENT
  
  ### Loop thru stages
  # For clean up
  stage.keep <- c(ls(), 'stage.keep', 'stage.loop')
  
  # Loop!
  for(stage.loop in rev(names(stage.sample.labels))){ ## UNCOMMENT
    # stage.loop = names(stage.sample.labels)[11] ## UNCOMMENT
    
    # Clean up
    rm(list = ls()[! ls() %in% stage.keep])
    gc()
    
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
    
    
    ##
    ## Metadata
    ##
    
    save.plot.dir.main <- paste0(output.path, 'figures/8c - plot picture naming cross-validated prediction time series - 50ms/')
    zoom <- 1.8
    
    
    ##
    ## Plot: Means of patients
    ##
    
    # Get plot data
    ys <- combined.pn.stats$accuracy
    xs <- combined.pn.stats$sample.label
    sig.windows <- combined.pn.sig.windows
    shuffle.mean <- combined.pn.stats$shuffle.mean
    shuffle.error.upper <- combined.pn.stats$shuffle.95ci.upper
    shuffle.error.lower <- combined.pn.stats$shuffle.95ci.lower
    
    # Y limits
    all.y.vals <- c(unlist(ys),
                    unlist(shuffle.error.upper),
                    unlist(shuffle.error.lower))
    y.limits <- range(all.y.vals)
    y.limits[1] <- floor(100 * y.limits[1] / 2) / 100 * 2
    y.limits[2] <- ceiling(100 * y.limits[2] / 2) / 100 * 2
    
    # Standardize a bit
    y.limits[1] <- min(y.limits[1], .14)
    y.limits[2] <- max(y.limits[2], .22)
    
    for(theme.loop in c('black','white')){ ## UNCOMMENT
      # theme.loop = c('white','black')[1] ## UNCOMMENT
      
      # Save dir
      save.means.dir <- paste0(save.plot.dir.main, 
                               'mean accuracy/',
                               theme.loop,'/',
                               model.loop,'/')
      dir.create(save.means.dir, showWarnings = FALSE, recursive = TRUE)
      
      # Save individual plot
      pdf(paste0(save.means.dir, 
                 'PN CV prediction accuracy - ',
                 roi.loop,' - ',
                 stage.loop,'.pdf'),
          width = 7, 
          height = 4)
      
      # Plot
      plot.time.series(.y.values = ys,
                       .x.values = xs,
                       .sampling.rate = 256,
                       # .colors = colors[[theme.loop]]$rainbow['pink','hex'],
                       .y.lwd = 3,
                       .y.limits = y.limits,
                       .y.ticks = c(y.limits[1], 
                                    mean(y.limits),
                                    y.limits[2]),
                       .y.tick.labels = c(y.limits[1], 
                                          "",
                                          y.limits[2]),
                       .y.label = '',
                       .x.label = '',
                       .x.limits = c(-median.rt.times['pn'] - 50, 250),
                       .x.ticks = c(-median.rt.times['pn'], 0, 250),
                       .shuffle.dist.mean = shuffle.mean,
                       .shuffle.dist.error.bars.upper = shuffle.error.upper - shuffle.mean,
                       .shuffle.dist.error.bars.lower = shuffle.mean - shuffle.error.lower,
                       .shuffle.dist.color = 'grey',
                       .polygons.x = unlist(stage.ranges[stage.loop, c('start.time','end.time')]),
                       .polygons.color = ifelse(theme.loop == 'white', 'black', 'white'),
                       .sig.windows = sig.windows,
                       .sig.color = colors[[theme.loop]]$rainbow_bright['pink','hex'],
                       # .sig.color = roi.colors[['rostral.caudal']][roi.loop],
                       show.sig.bars = FALSE,
                       show.sig.highlights = TRUE,
                       .center.sig.bars.vertically = FALSE,
                       .title = '',
                       .theme = theme.loop,
                       .background = rgb(1,1,1,0),
                       .margin = c(3,3,1.5,1.5),
                       .zoom = zoom)
      
      # Save
      dev.off()
      
    }#; rm(theme.loop)
    
    ##
    ## Individual patients
    ##
    
    for(patient in names(patient.pn.stats)){ ## UNCOMMENT
      # patient = names(patient.pn.stats)[1] ## UNCOMMENT
      
      ys <- patient.pn.stats[[patient]]$accuracy
      xs <- patient.pn.stats[[patient]]$sample.label
      sig.windows <- patient.pn.sig.windows[[patient]]
      shuffle.mean <- patient.pn.stats[[patient]]$shuffle.mean
      shuffle.error.upper <- patient.pn.stats[[patient]]$shuffle.95ci.upper
      shuffle.error.lower <- patient.pn.stats[[patient]]$shuffle.95ci.lower
      
      # Y limits
      all.y.vals <- c(unlist(ys),
                      unlist(shuffle.error.upper),
                      unlist(shuffle.error.lower))
      y.limits <- range(all.y.vals)
      y.limits[1] <- floor(100 * y.limits[1] / 2) / 100 * 2
      y.limits[2] <- ceiling(100 * y.limits[2] / 2) / 100 * 2
      
      # Standardize a bit
      y.limits[1] <- min(y.limits[1], .12)
      y.limits[2] <- max(y.limits[2], .28)
      
      for(theme.loop in c('black','white')){ ## UNCOMMENT
        # theme.loop = c('white','black')[1] ## UNCOMMENT
        
        # Save dir
        save.patients.dir <-
          paste0(save.plot.dir.main,
                 'individual patients/',
                 theme.loop,'/',
                 model.loop,'/')
        dir.create(save.patients.dir, showWarnings = FALSE, recursive = TRUE)
        
        # Sort grids into folders by whether any significant results this patient/roi/stage
        current.sig <- nrow(sig.windows) > 0
        
        ## Just don't plot if not sig
        if(current.sig){
          
        current.save.patients.sig.dir <- 
          paste0(save.patients.dir,
                 ifelse(current.sig,'significant/', 'not significant/'))
        dir.create(current.save.patients.sig.dir, showWarnings = FALSE, recursive = TRUE)
        
        # Begin figure save
        pdf(paste0(current.save.patients.sig.dir,
                   'PN CV prediction accuracy - ',
                   roi.loop,' - ',
                   stage.loop,' - ',
                   patient,'.pdf'),
            width = 7, height = 4)
        
        ### Plot
        plot.time.series(.y.values = ys,
                         .x.values = xs,
                         .sampling.rate = 256,
                         # .colors = colors[[theme.loop]]$rainbow['pink','hex'],
                         .y.lwd = 3,
                         # .y.limits = y.limits,
                         .y.limits = c(.12, .3),
                         .y.ticks = c(.15, .3),
                         # .y.ticks = c(y.limits[1] y.limits[2]),
                         .y.label = '',
                         .x.label = '',
                         .x.limits = c(-median.rt.times['pn'] - 50, 500),
                         .x.ticks = c(-median.rt.times['pn'], 0, 500),
                         .shuffle.dist.mean = shuffle.mean,
                         .shuffle.dist.error.bars.upper = shuffle.error.upper - shuffle.mean,
                         .shuffle.dist.error.bars.lower = shuffle.mean - shuffle.error.lower,
                         .shuffle.dist.color = 'grey',
                         .polygons.x = unlist(stage.ranges[stage.loop, c('start.time','end.time')]),
                         .polygons.color = ifelse(theme.loop == 'white', 'black', 'white'),
                         .sig.windows = sig.windows,
                         .sig.color = colors[[theme.loop]]$rainbow_bright['pink','hex'],
                         # .sig.color = roi.colors[['rostral.caudal']][roi.loop],
                         show.sig.bars = FALSE,
                         show.sig.highlights = TRUE,
                         .center.sig.bars.vertically = FALSE,
                         .title = '',
                         .theme = theme.loop,
                         .background = rgb(1,1,1,0),
                         .margin = c(3,3,1.5,1.5),
                         .zoom = zoom)
        
        # Save
        dev.off()
        
        } # if(current.sig){}
        
      }#; rm(theme.loop)
      
    }#; rm(patient.loop)
    
    
    
    # Progress update
    message('Completed: ',
            model.loop,' - ',
            band.loop,' - ',
            roi.loop,' - ',
            stage.loop,'.',
            Sys.time())
    
    
  }; rm(stage.loop)
}; rm(roi.loop)
}; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


