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
# Adjust color transparency
source(paste0(path,'/analysis/R/functions/adjust_transparency.R'))
# Get sig windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Get Fisher z-transformation functions
source(paste0(path,'/analysis/R/functions/fisher_z_transform.R'))
# Plot time series stack with significance
source(paste0(path,'/analysis/R/functions/plot_time_series_stack.R'))



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


### Loop thru models
model.types <- list.files(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/'))
for(model.loop in model.types){ ## UNCOMMENT
# model.loop = model.types[2] ## UNCOMMENT


###
### Get data: Peak accuracy vs. central train time
###

# Storage
time.series <- 
  sig.windows <-
  train.time.dfs <-
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
    
    
    ### Collapsed across patients
    ## If any sig windows
    if(nrow(combined.pn.sig.windows) > 0){
      
      ## Get data
      collapsed.label <- paste0(roi.loop,'.',stage.loop)
      
      # Stats
      time.series[['collapsed.across.patients']][[collapsed.label]] <- combined.pn.stats
      sig.windows[['collapsed.across.patients']][[collapsed.label]] <- combined.pn.sig.windows
      
      # Train times
      train.time.dfs[['collapsed.across.patients']][[collapsed.label]] <- train.time.df.template
      train.time.dfs[['collapsed.across.patients']][[collapsed.label]]$start.time <- 
        stage.ranges[stage.loop, 'start.time']
      train.time.dfs[['collapsed.across.patients']][[collapsed.label]]$end.time <- 
        stage.ranges[stage.loop, 'end.time']
      
      # Clean up
      rm(collapsed.label)
      
    } # if(nrow(combined.pn.sig.windows) > 0){
    
    ### By patient
    for(patient in names(patient.pn.stats)){
      # patient = names(patient.pn.stats)[2]
      
      ## If any sig windows
      if(nrow(patient.pn.sig.windows[[patient]]) > 0){
        
        ## Get data
        individual.label <- paste0(roi.loop,'.',stage.loop,'.',patient)
        
        # Stats
        time.series[['individually.by.patient']][[individual.label]] <- patient.pn.stats[[patient]]
        sig.windows[['individually.by.patient']][[individual.label]] <- patient.pn.sig.windows[[patient]]
        
        # Train times
        train.time.dfs[['individually.by.patient']][[individual.label]] <- train.time.df.template
        train.time.dfs[['individually.by.patient']][[individual.label]]$start.time <- 
          stage.ranges[stage.loop, 'start.time']
        train.time.dfs[['individually.by.patient']][[individual.label]]$end.time <- 
          stage.ranges[stage.loop, 'end.time']
        
        # Clean up
        rm(individual.label)
        
        
        ### ROI plot data
        roi.label <- paste0(stage.loop,'.',patient)
        
        # Stats
        time.series[[roi.loop]][[roi.label]] <- patient.pn.stats[[patient]]
        sig.windows[[roi.loop]][[roi.label]] <- patient.pn.sig.windows[[patient]]
        
        # Train times
        train.time.dfs[[roi.loop]][[roi.label]] <- train.time.df.template
        train.time.dfs[[roi.loop]][[roi.label]]$start.time <- 
          stage.ranges[stage.loop, 'start.time']
        train.time.dfs[[roi.loop]][[roi.label]]$end.time <- 
          stage.ranges[stage.loop, 'end.time']
        
        # Clean up
        rm(roi.label)
        
      } # if(nrow(patient.pn.sig.windows[[patient]]) > 0){
      
    }; rm(patient)
  }; rm(stage.loop)
}; rm(roi.loop)
message('Done getting data! ',Sys.time()) ## < 1 min


### Clean up
train.time.dfs <- 
  lapply(train.time.dfs, function(x){
    x <- bind_rows(x, .id = 'label')
    rownames(x) <- x$label
    x$label <- NULL
    return(x)
  })


### Plot parameters
save.plot.dir.main <- paste0(output.path, 'figures/8f - stack prediction time series/')
zoom <- 1.8
text.size.big <- 1.8
text.size.med <- 1.6
text.size.small <- 1.4


##
## Plot: Accuracy
##

for(theme.loop in c('white','black')){ ## UNCOMMENT
  # theme.loop = c('white','black')[1] ## UNCOMMENT
  
  for(data.loop in names(train.time.dfs)){
    # data.loop = names(train.time.dfs)[2]
  
    # Save dir
    save.dir.all.rois.all.stages <- 
      paste0(save.plot.dir.main,
             'all rois all stages/',
             theme.loop,'/')
    dir.create(save.dir.all.rois.all.stages, showWarnings = FALSE, recursive = TRUE)
    
    # # Save
    plot.height <- ifelse(data.loop %in% rois, 5, 20)
    plot.height <- ifelse(data.loop == 'SMC', 10, plot.height)
    pdf(paste0(save.dir.all.rois.all.stages, 
               data.loop,
               ' - picture naming cross validation accuracies stacked - ',
               model.loop,' - n_sig=',nrow(train.time.dfs[[data.loop]]),'.pdf'),
        width = 6, height = plot.height)
    plot.time.series.stack(.cluster.times.df = train.time.dfs[[data.loop]],
                           .sig.windows = sig.windows[[data.loop]],
                           .time.series.y = lapply(time.series[[data.loop]], function(x){x$accuracy - 1/6}),
                           .time.series.x = lapply(time.series[[data.loop]], function(x){x$sample.label}),
                           .sampling.rate = 256,
                           .title = '',
                           .y.label = '',
                           .x.label = '',
                           .x.limits = c(-median.rt.times['pn'],500),
                           .x.ticks = c(-median.rt.times['pn'],0,500),
                           # .sig.colors = shades::saturation(character.colors['nurse','hex'], values = delta(.8)),
                           # .time.series.color = colorspace::lighten(character.colors['nurse','hex'], amount = .5, method = 'relative'),
                           .sig.colors = ifelse(data.loop %in% rois,
                                                roi.colors$rostral.caudal[data.loop],
                                                colors[[theme.loop]]$rainbow['pink','hex']),
                           .time.series.color = rgb(.5,.5,.5),
                           .train.time.color = 'black',
                           .omit.clusters.with.no.sig = FALSE,
                           show.sig.bars = FALSE,
                           show.sig.highlights = TRUE,
                           .time.series.scale.factor = 3,
                           .time.series.lwd = 1,
                           .cluster.line.width.zoom = .5, # factor for width of sig bar/highlight
                           show.density = TRUE,
                           .density.line.color = ifelse(data.loop %in% rois,
                                                        roi.colors$rostral.caudal[data.loop],
                                                        colors[[theme.loop]]$rainbow['pink','hex']),
                           .density.height.proportion = 20 / plot.height * .05,#ifelse(data.loop %in% rois, .15, .05),
                           .density.variable = 'significance',
                           show.row.labels = FALSE,#plot.grid.loop,
                           plot.guidelines = FALSE,
                           .background = rgb(1,1,1,0),
                           .theme = theme.loop,
                           .margin = c(3,1.5,1,1.5),
                           .zoom = 1.8)
    
  dev.off()
  
  }#; rm(data.loop)
}#; rm(theme.loop)



}; rm(model.loop)
# }; rm(band.loop)










# Finish!
message('Script completed successfully. ',Sys.time())


