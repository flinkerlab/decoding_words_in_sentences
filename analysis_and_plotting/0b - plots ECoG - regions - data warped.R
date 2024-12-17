### Plot electrode high gamma and encoding importance over time to complement brain plots
### Feb 2023
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
library('viridis') # for viridis color schemes
library('pals') # for ocean.balance, ocean.delta, ocean.curl colors
#pal.bands(ocean.balance, ocean.delta, ocean.curl)
library('scico') # for palettes cork, vik, lisbon, tofino, berlin, roma
library('reticulate') # for Python
library('stringr') # for str_split()
library('lme4') # for lmer()
library('lmerTest') # for lmer() that gives p-value
library('NMF') # for non-negative matrix factorization

# Clean up
rm(list=ls())
cat("\014")
message("Begin plotting ECoG for each patient, electrode, task, and time-lock. ",Sys.time())

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 2
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 6
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
    n.cores.to.use = 18
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
  n.cores.to.use = 28
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
# Adjust color transparency
source(paste0(path,'/analysis/R/functions/adjust_transparency.R'))
# Plot time series
source(paste0(path,'/analysis/R/functions/plot_time_series.R'))
# Add a line of text to plot with varying font color
source(paste0(path,'/analysis/R/functions/add_text_line_multiple_colors.R'))
# Get significant windows
source(paste0(path,'/analysis/R/functions/get_significant_windows.R'))
# Good trial subset functions
source(paste0(path,'/analysis/R/functions/just_good_trial_functions.R'))


### Set up elec info
# Read in elecs
elec.info <- 
  read.csv(paste0(path,
                  '/analysis/R/brain plots/ecog/output/data/elec info/patients - combined/row_labels_and_localizations.csv'))
# Subset to just good elecs:
elec.info <- droplevels(subset(elec.info,
                               (bad_elec == 0) &
                                 (active.stim.locked == 1) &
                                 (active.prod.locked == 1) &
                                 (bad_localization == 0) &
                                 (! region_clinical %in% c('NaN','Unknown',''))))
use.these.elecs <- elec.info$patient_elec


### Loop thru frequency bands
keep <- c(ls(), 'keep', 'band.loop')
# for(band.loop in c('high_gamma','beta')){
band.loop = c('high_gamma','beta')[1]

rm(list = ls()[which(! ls() %in% keep)])

# Save directory
output.path <- paste0(path, 'analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/',band.loop,'/')


### Load ROIs
load(paste0(output.path, 'data/0a - definitions/roi metadata/ROIs and temporal lobe splits.RData'))

# Turn elecs into a dataframe for easy lookup
roi.elecs.df <- data.frame('elec' = unlist(roi.elecs))
roi.elecs.df$region.category <- rownames(roi.elecs.df)
for(i in 0:9){roi.elecs.df$region.category <- gsub(i, '', roi.elecs.df$region.category)}


### Metadata
patients = c('NY765',
             'NY794',
             'NY798',
             'NY799',
             'NY829',
             # 'NY834',
             'NY837',
             'NY844',
             'NY857',
             'NY863',
             'NY869')
tasks <- c('pn','sp','lp') # picture naming, sentence production, list production

# Colors
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))

### Begin reading in stats and stuff
## Takes ~2 hours with stats so if possible load previous results
save.data.dir <- paste0(output.path,'data/0b - plot ECoG - ROIs - warped data/roi stats/')
save.data.file <- 'ECoG ROI means.RData'
if(file.exists(paste0(save.data.dir, save.data.file))){
  load(paste0(save.data.dir, save.data.file))
}else{
  
  # Initialize storage of stats etc
  region.means <- list()
  region.ses <- list()
  # sentence.vs.list.sig.windows <- list()
  # sentence.vs.naming.sig.windows <- list()
  # list.vs.naming.sig.windows <- list()
  
  message("Begin data set up for loops: ",band.loop,'. ', Sys.time())
  
  ## Read in electrode data
  elec.data <- list()
  trial.info <- list()
  message('Attaching warped data...')
  attach(paste0(path,
                'analysis/R/warp time series to standard RT/simple linear stretch to median RT/output - ',
                band.loop,
                '/data/elec data with RTs stretched to global median by task - bad trials and .025 to .95 RT outliers excluded.RData')
  ) # Loads: median.rt.samples, patient.trial.info, stretch.data, stretch.sample.labels
  for(task.loop in tasks){
    elec.data[[task.loop]] <- unlist(lapply(stretch.data, function(x){x[[task.loop]]}), recursive = FALSE)
    names(elec.data[[task.loop]]) <- gsub("^[^.]*\\.", "", names(elec.data[[task.loop]]))
    trial.info[[task.loop]] <- lapply(patient.trial.info, function(x){x[[task.loop]]})
  }; rm(task.loop)
  median.rt.samples <- median.rt.samples
  stretch.sample.labels <- stretch.sample.labels
  detach()
  message('...detached!')
  
  # Subset to just time samples you care about
  sample.labels <- list('locked_to_production_onset' = list(),
                        'locked_to_stimulus_onset' = list())
  for(task.loop in tasks){
    # task.loop = tasks[1]
    
    prod.locked.start.sample <- (-median.rt.samples[task.loop] - time.convert(100, "times", "samples"))
    prod.locked.end.sample <- time.convert(1100, "times", "samples")
    sample.labels$locked_to_production_onset[[task.loop]] <- 
      time.convert(prod.locked.start.sample:prod.locked.end.sample, "samples", "sample.labels")
    sample.labels$locked_to_stimulus_onset[[task.loop]] <- 
      stretch.sample.labels[[task.loop]]$locked_to_stimulus_onset[
        which(stretch.sample.labels[[task.loop]]$locked_to_production_onset %in% 
                sample.labels$locked_to_production_onset[[task.loop]])]
    
    elec.data[[task.loop]] <- lapply(elec.data[[task.loop]], function(x){x[, sample.labels[['locked_to_production_onset']][[task.loop]]]})
  }; rm(task.loop)
  
  # Get metadata
  times <- lapply(sample.labels, function(x){lapply(x, function(y){time.convert(y, 'sample.labels', 'times')})})
  
  
  ###
  ### Get ECoG means and SEs by region
  ###
  
  # Initialize storage for ECoG means
  region.elec.data <- list()
  region.elec.means <- list()
  region.means <- list()
  region.ses <- list()
  
  # Loop thru tasks
  for(task.loop in tasks){
    # task.loop = tasks[3]
    region.elec.data[[task.loop]] <- list()
    region.elec.means[[task.loop]] <- list()
    region.means[[task.loop]] <- list()
    region.ses[[task.loop]] <- list()
    
    # Loop thru region categories
    for(region.loop in 1:length(roi.categories)){
      # region.loop = 7
      current.region <- names(roi.categories)[region.loop]
      current.elecs <- roi.elecs.df$elec[roi.elecs.df$region.category == current.region]
      # Get rid of elecs with no data for this task
      current.elecs <- current.elecs[current.elecs %in% names(elec.data[[task.loop]])]
      
      # Store data for stats below
      region.elec.data[[task.loop]][[current.region]] <- elec.data[[task.loop]][current.elecs]
      current.n.rows <- sapply(region.elec.data[[task.loop]][[current.region]], nrow)
      # Add electrode and task to each elec dataframe
      for(elec.loop in current.elecs){
        # elec.loop <- current.elecs[1]
        region.elec.data[[task.loop]][[current.region]][[elec.loop]] <- 
          cbind(data.frame('elec' = rep(elec.loop, times = current.n.rows[elec.loop]),
                           'task' = rep(task.loop, times = current.n.rows[elec.loop])),
                region.elec.data[[task.loop]][[current.region]][[elec.loop]])
      }; rm(elec.loop)
      # Collapse elec dataframes for this region into one big one
      region.elec.data[[task.loop]][[current.region]] <-
        do.call(rbind, region.elec.data[[task.loop]][[current.region]])
      
      # Get means for each electrode
      region.elec.means[[task.loop]][[current.region]] <- 
        data.frame(bind_rows(lapply(elec.data[[task.loop]][current.elecs],
                                    function(x){apply(x, 2, mean)})))
      # Get rid of electrodes where there's no data
      region.elec.means[[task.loop]][[current.region]] <- 
        na.omit(region.elec.means[[task.loop]][[current.region]])
      # SE across elecs
      region.ses[[task.loop]][[current.region]] <- 
        apply(region.elec.means[[task.loop]][[current.region]], 2, sd) /
        sqrt(nrow(region.elec.means[[task.loop]][[current.region]]))
      # Average of elecs
      region.means[[task.loop]][[current.region]] <- 
        apply(region.elec.means[[task.loop]][[current.region]], 2, mean)
    }; rm(region.loop)
  }; rm(task.loop)
  
  
  # ### Get stats on differences between tasks
  # for(lock.loop in names(sample.labels)){
  #   # lock.loop = names(sample.labels)[1]
  #   
  #   # Storage
  #   sentence.vs.list <- list()
  #   sentence.vs.naming <- list()
  #   list.vs.naming <- list()
  #   sentence.vs.list.sig.windows[[lock.loop]] <- list()
  #   sentence.vs.naming.sig.windows[[lock.loop]] <- list()
  #   list.vs.naming.sig.windows[[lock.loop]] <- list()
  #   
  #   # Loop thru regions
  #   for(region.loop in names(roi.categories)){
  #     # region.loop = 'MTG'#names(roi.categories)[1]
  #     message("Beginning EGoG activity models for ",region.loop," (",which(names(roi.categories) == region.loop)," of ",length(roi.categories)," regions), ",lock.loop,". ",Sys.time())
  #     
  #     # Storage
  #     sentence.vs.list[[region.loop]] <- list()
  #     sentence.vs.naming[[region.loop]] <- list()
  #     list.vs.naming[[region.loop]] <- list()
  #     sentence.vs.list.sig.windows[[lock.loop]][[region.loop]] <- list()
  #     sentence.vs.naming.sig.windows[[lock.loop]][[region.loop]] <- list()
  #     list.vs.naming.sig.windows[[lock.loop]][[region.loop]] <- list()
  #     
  #     # Loop thru time samples
  #     for(sample.loop in sample.labels[[lock.loop]][['sp']]){ # sp has the most sample labels
  #       # sample.loop = sample.labels[[lock.loop]][['sp']][1]
  #       
  #       if(lock.loop == 'locked_to_production_onset'){
  #         # If locked to production onset
  #         current.prod.locked.sample.sp <- sample.loop
  #         current.prod.locked.sample.lp <- sample.labels[['locked_to_production_onset']][['lp']][
  #           which(sample.labels[['locked_to_production_onset']][['lp']] == sample.loop)]
  #         current.prod.locked.sample.pn <- sample.labels[['locked_to_production_onset']][['pn']][
  #           which(sample.labels[['locked_to_production_onset']][['pn']] == sample.loop)]
  #       }else{ # if(lock.loop == 'locked_to_production_onset'){
  #         # If locked to stimulus onset
  #         current.prod.locked.sample.sp <- 
  #           sample.labels[['locked_to_production_onset']][['sp']][
  #             which(sample.labels[['locked_to_stimulus_onset']][['sp']] == sample.loop)]
  #         current.prod.locked.sample.lp <- 
  #           sample.labels[['locked_to_production_onset']][['lp']][
  #             which(sample.labels[['locked_to_stimulus_onset']][['lp']] == sample.loop)]
  #         current.prod.locked.sample.pn <- 
  #           sample.labels[['locked_to_production_onset']][['pn']][
  #             which(sample.labels[['locked_to_stimulus_onset']][['pn']] == sample.loop)]
  #       } # if(lock.loop == 'locked_to_production_onset'){
  #       
  #       ### Sentence vs. list
  #       if(length(c(current.prod.locked.sample.sp, current.prod.locked.sample.lp)) == 2){
  #         current.model <-
  #           wilcox.test(region.elec.data$sp[[region.loop]][, current.prod.locked.sample.sp],
  #                       region.elec.data$lp[[region.loop]][, current.prod.locked.sample.lp],
  #                       alternative = 'greater')
  #         
  #         # Store p-value
  #         sentence.vs.list[[region.loop]][[sample.loop]] <-
  #           current.model$p.value
  #         rm(current.model)
  #       } # if(length(c(current.prod.locked.sample.sp, current.prod.locked.sample.lp)) == 2){
  #       
  #       
  #       ### Sentence vs. naming
  #       if(length(c(current.prod.locked.sample.sp, current.prod.locked.sample.pn)) == 2){
  #         current.model <- 
  #           wilcox.test(region.elec.data$sp[[region.loop]][, current.prod.locked.sample.sp],
  #                       region.elec.data$pn[[region.loop]][, current.prod.locked.sample.pn],
  #                       alternative = 'greater')
  #         
  #         # Store p-value
  #         sentence.vs.naming[[region.loop]][[sample.loop]] <- 
  #           current.model$p.value
  #         rm(current.model)
  #       } # if(sample.loop %in% sample.labels[[lock.loop]][['pn']]){
  #       
  #       ### List vs. naming
  #       if(length(c(current.prod.locked.sample.lp, current.prod.locked.sample.pn)) == 2){
  #         current.model <- 
  #           wilcox.test(region.elec.data$lp[[region.loop]][, current.prod.locked.sample.lp],
  #                       region.elec.data$pn[[region.loop]][, current.prod.locked.sample.pn],
  #                       alternative = 'greater')
  #         
  #         # Store p-value
  #         list.vs.naming[[region.loop]][[sample.loop]] <-
  #           current.model$p.value
  #         rm(current.model)
  #       } # if(sample.loop %in% sample.labels[[lock.loop]][['pn']]){
  #       
  #     }; rm(sample.loop)
  #     
  #     
  #     ### Fix p-values etc.
  #     # Set up sig windows lists
  #     sentence.vs.list.sig.windows[[lock.loop]][[region.loop]] <- list()
  #     sentence.vs.naming.sig.windows[[lock.loop]][[region.loop]] <- list()
  #     list.vs.naming.sig.windows[[lock.loop]][[region.loop]] <- list()
  #     
  #     # Reorganize p values
  #     sentence.vs.list[[region.loop]] <- 
  #       data.frame('p' = t(bind_rows(sentence.vs.list[[region.loop]])))
  #     sentence.vs.naming[[region.loop]] <- 
  #       data.frame('p' = t(bind_rows(sentence.vs.naming[[region.loop]])))
  #     list.vs.naming[[region.loop]] <- 
  #       data.frame('p' = t(bind_rows(list.vs.naming[[region.loop]])))
  #     
  #     ## Fill in sig windows
  #     for(alpha.loop in c(0.01, 0.05)){
  #       # alpha.loop = c(.01, .05)[1]
  #       
  #       alpha.label = paste0('alpha=',alpha.loop)
  #       
  #       ## Sentence vs. list
  #       sentence.vs.list[[region.loop]][,alpha.label] <- 
  #         as.numeric(sentence.vs.list[[region.loop]]$p < alpha.loop)
  #       sentence.vs.list.sig.windows[[lock.loop]][[region.loop]][[alpha.label]] <- 
  #         get.significant.windows(sentence.vs.list[[region.loop]][,alpha.label],
  #                                 .sample.labels = row.names(sentence.vs.list[[region.loop]]),
  #                                 include.duration = TRUE,
  #                                 # .exclude.sig.durations.under.ms = 100,
  #                                 output.class = 'data.frame')
  #       
  #       ## Sentence vs. naming
  #       sentence.vs.naming[[region.loop]][,alpha.label] <- 
  #         as.numeric(sentence.vs.naming[[region.loop]]$p < alpha.loop)
  #       sentence.vs.naming.sig.windows[[lock.loop]][[region.loop]][[alpha.label]] <- 
  #         get.significant.windows(sentence.vs.naming[[region.loop]][,alpha.label],
  #                                 .sample.labels = row.names(sentence.vs.naming[[region.loop]]),
  #                                 include.duration = TRUE,
  #                                 # .exclude.sig.durations.under.ms = 100,
  #                                 output.class = 'data.frame')
  #       
  #       ## List vs. naming
  #       list.vs.naming[[region.loop]][,alpha.label] <- 
  #         as.numeric(list.vs.naming[[region.loop]]$p < alpha.loop)
  #       list.vs.naming.sig.windows[[lock.loop]][[region.loop]][[alpha.label]] <- 
  #         get.significant.windows(list.vs.naming[[region.loop]][,alpha.label],
  #                                 .sample.labels = row.names(list.vs.naming[[region.loop]]),
  #                                 include.duration = TRUE,
  #                                 # .exclude.sig.durations.under.ms = 100,
  #                                 output.class = 'data.frame')
  #       
  #     }; rm(alpha.loop)
  #     
  #   }; rm(region.loop) # duration
  #   
  #   beep(2)
  # }; rm(lock.loop)
  
  ## Save
  # Make a list of things you want to save:
  save.these <- c('region.means',
                  'region.ses',
                  'sample.labels',
                  'times',
                  'median.rt.samples')
  # 'sentence.vs.list.sig.windows',
  # 'sentence.vs.naming.sig.windows',
  # 'list.vs.naming.sig.windows')
  
  # Save
  dir.create(save.data.dir, showWarnings = FALSE, recursive = TRUE)
  save(list = save.these,
       file = paste0(save.data.dir, save.data.file))
  
} # if(file.exists(paste0(save.data.dir, 'ecog wilcoxon task tests.RData'))){



###
### Plots
###

### For region means and region ses, duplicate with sample labels
region.means <- list('locked_to_production_onset' = region.means,
                     'locked_to_stimulus_onset' = region.means)
region.ses <- list('locked_to_production_onset' = region.ses,
                   'locked_to_stimulus_onset' = region.ses)
# Update stim-locked labels
for(task.loop in tasks){
  # task.loop = tasks[1]
  for(region.loop in names(region.means[[1]][[task.loop]])){
    names(region.means$locked_to_stimulus_onset[[task.loop]][[region.loop]]) <-
      sample.labels$locked_to_stimulus_onset[[task.loop]]
    names(region.ses$locked_to_stimulus_onset[[task.loop]][[region.loop]]) <-
      sample.labels$locked_to_stimulus_onset[[task.loop]]
  }; rm(region.loop)
}; rm(task.loop)

### Loop thru TL splits
for(tl.split.loop in names(temporal.lobe.splits)){
  # tl.split.loop = names(temporal.lobe.splits)[1]
  
  ## Get y limits across prod-locked and stim-locked data
  y.limits <- list()
  regions.with.data <- names(which(sapply(region.means$locked_to_production_onset$pn, length) != 0))
  plot.regions <- regions.with.data[regions.with.data %in% temporal.lobe.splits[[tl.split.loop]]]
  for(region.loop in plot.regions){
    # region.loop = plot.regions[1]
    y.limits[[region.loop]] <- c()
    
    for(lock.loop in c('locked_to_production_onset','locked_to_stimulus_onset')){
      # lock.loop = c('locked_to_production_onset','locked_to_stimulus_onset')[2]
      
      y.limits[[region.loop]] <-
        c(y.limits[[region.loop]],
          region.means[[lock.loop]]$pn[[region.loop]] + region.ses[[lock.loop]]$pn[[region.loop]],
          region.means[[lock.loop]]$pn[[region.loop]] - region.ses[[lock.loop]]$pn[[region.loop]],
          region.means[[lock.loop]]$lp[[region.loop]] + region.ses[[lock.loop]]$lp[[region.loop]],
          region.means[[lock.loop]]$lp[[region.loop]] - region.ses[[lock.loop]]$lp[[region.loop]],
          region.means[[lock.loop]]$sp[[region.loop]] + region.ses[[lock.loop]]$sp[[region.loop]],
          region.means[[lock.loop]]$sp[[region.loop]] - region.ses[[lock.loop]]$sp[[region.loop]])
    }; rm(lock.loop)
    
    y.limits[[region.loop]] <-
      c('min' = min(y.limits[[region.loop]]),
        'max' = max(y.limits[[region.loop]]))
  }; rm(region.loop)
  
  # Y limits
  pub.regions.y.limits <- 
    c('min' = min(unlist(y.limits[plot.regions])),
      'max' = max(unlist(y.limits[plot.regions])))
  
  # Times
  median.rt.times <- time.convert(median.rt.samples, 'samples', 'times', 512)
  
  # Loop thru black and white backgrounds
  for(theme.loop in c('black','white')){
    # theme.loop = 'white'
    
    for(lock.loop in c('locked_to_production_onset','locked_to_stimulus_onset')){
      # lock.loop = c('locked_to_production_onset','locked_to_stimulus_onset')[2]
      
      ##
      ## No statistical comparisons plotted - PDF FOR PUBLICATION (no axis labels etc.)
      ##
      
      for(smooth.loop in c(TRUE, FALSE)){
        
        # Save directory
        save.dir <- 
          paste0(output.path,
                 'figures/0b - plot ECoG - ROIs - warped data/by region/no comparisons - for manuscript/',
                 theme.loop,'/',
                 'temporal lobe split along ',tl.split.loop,' axis/')
        dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
        
        for(region.loop in 1:length(plot.regions)){
          # region.loop = 1
          current.region <- plot.regions[region.loop]
          
          # Any electrodes in this region?
          if((length(region.means[[lock.loop]]$pn[[current.region]]) > 0) &
             ((! any(is.na(region.means[[lock.loop]]$pn[[current.region]])))) &
             (length(region.ses[[lock.loop]]$pn[[current.region]]) > 0) &
             ((! any(is.na(region.ses[[lock.loop]]$pn[[current.region]]))))){
            
            current.y <- list(region.means[[lock.loop]]$pn[[current.region]],
                              region.means[[lock.loop]]$lp[[current.region]],
                              region.means[[lock.loop]]$sp[[current.region]])
            current.x <- times[[lock.loop]][c('pn','lp','sp')]
            current.error.bar <- list(region.ses[[lock.loop]]$pn[[current.region]],
                                      region.ses[[lock.loop]]$lp[[current.region]],
                                      region.ses[[lock.loop]]$sp[[current.region]])
            
            
            # Smooth
            if(smooth.loop){
              half.n.smoothing.samples <- time.convert(50, 'times', 'samples', 256)
              current.y <- lapply(current.y, function(x){
                smoothing(x, 
                          n.samples.pre = half.n.smoothing.samples, 
                          na.pad = FALSE)})
              current.error.bar <- lapply(current.error.bar, function(x){
                smoothing(x, 
                          n.samples.pre = half.n.smoothing.samples, 
                          na.pad = FALSE)})
              current.x <- lapply(current.x, function(x){
                x[(half.n.smoothing.samples + 1):(length(x) - half.n.smoothing.samples)]
              })
            } # if(smooth.loop)
            
            
            pdf(paste0(save.dir,current.region,'_',lock.loop,' - ', ifelse(smooth.loop, 'smoothed', 'unsmoothed'),'.pdf'),
                width = 9, height = 4.5)
            # par(mfrow = c(1,1),
                # oma = c(0,0,1,0))
            
            plot.time.series(.y.values = current.y,
                             .x.values = current.x,
                             #.y.limits = pub.regions.y.limits,
                             .y.limits.min.at.least = pub.regions.y.limits['min'],
                             .y.limits.max.at.least = pub.regions.y.limits['max'],
                             # .title = current.region,
                             .y.label = '',#band.loop,
                             .y.ticks = list('locked_to_production_onset' = list('beta' = c(0,.5,1),
                                                                                 'high_gamma' = c(0,1))[[band.loop]],
                                             'locked_to_stimulus_onset' = list('beta' = c(-.25, 0, .25, .5),
                                                                               'high_gamma' = c(0,1))[[band.loop]])[[lock.loop]],
                             .y.tick.labels = list('locked_to_production_onset' = list('beta' = c("0",'.5','1'),
                                                                                       'high_gamma' = c("0","1"))[[band.loop]],
                                                   'locked_to_stimulus_onset' = list('beta' = c('', '0', ' ', '.5'),
                                                                                     'high_gamma' = c("0","1"))[[band.loop]])[[lock.loop]],
                             .x.ticks = list('locked_to_production_onset' = seq(-1000, 1000, by = 500),
                                             'locked_to_stimulus_onset' = seq(0, 2000, by = 500))[[lock.loop]],
                             .x.tick.labels = list('locked_to_production_onset' = c('-1000', '', '0', '', '1000'),
                                                   'locked_to_stimulus_onset' = c('0', '', '1000', '', '2000'))[[lock.loop]],
                             .x.label = '',#time (ms)',
                             .x.limits = list('locked_to_production_onset' = c(-time.convert(median.rt.samples['sp'], 'samples', 'times'), 1000),
                                              'locked_to_stimulus_onset' = c(0, 2000))[[lock.loop]],
                             show.t0 = (lock.loop == 'locked_to_production_onset'),
                             .horizontal.line.at = list('high_gamma' = NA,
                                                        'beta' = 0)[[band.loop]],
                             # show.y.axis = (lock.loop == 'locked_to_stimulus_onset'),
                             .colors = colors[[theme.loop]]$task_line_plots[c('pn','lp','sp'),'hex'],
                             .error.bars = current.error.bar,
                             .zoom = 1.8,
                             .margin = c(3.5,3.5,1.5,1.5),
                             .theme = theme.loop,
                             .background = ifelse(theme.loop == 'white', rgb(1,1,1,0), theme.loop))
            dev.off()
            
          } # if any elecs
        }; rm(region.loop)
        
      } # smooth.loop
      
      
      ##
      ## Plot ECoG for each region by task
      ## 
      
      
      # Save directory
      save.dir <- paste0(output.path,
                         'figures/0b - plot ECoG - ROIs - warped data/by task/no comparisons - for manuscript/',
                         theme.loop,'/',
                         'temporal lobe split along ',tl.split.loop,' axis/')
      dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
      
      x.limits <- list(
        'locked_to_stimulus_onset' = 
          lapply(as.list(time.convert(median.rt.samples, 'samples', 'times')), function(x){
            c(-50, round((x + 500) / 500) * 500)}),
        'locked_to_production_onset' =
          lapply(time.convert(median.rt.samples, 'samples', 'times'), function(x){
            c(-round((x + 500) / 500) * 500, 500)})
      )
      x.ticks <- lapply(x.limits, function(x){lapply(x, function(y){
        .samples <- time.convert(y, "times", "samples")
        .samples <- .samples[1]:.samples[2]
        .tick.times.superset <- seq(-2000, 2000, by = 500)
        .times <- time.convert(.samples, "samples", "times")
        .ticks <- .times[.times %in% .tick.times.superset]
        return(.ticks)
      })})
      
      for(smooth.loop in c(TRUE, FALSE)){
        for(task.loop in tasks){
          # task.loop = tasks[3]
          
          
          current.y <- region.means[[lock.loop]][[task.loop]][plot.regions]
          current.x <- times[[lock.loop]][task.loop]
          current.error.bar <- region.ses[[lock.loop]][[task.loop]][plot.regions]
          
          # Smooth
          if(smooth.loop){
            half.n.smoothing.samples <- time.convert(50, 'times', 'samples', 256)
            current.y <- lapply(current.y, function(x){
              smoothing(x, 
                        n.samples.pre = half.n.smoothing.samples, 
                        na.pad = FALSE)})
            current.error.bar <- lapply(current.error.bar, function(x){
              smoothing(x, 
                        n.samples.pre = half.n.smoothing.samples, 
                        na.pad = FALSE)})
            current.x <- lapply(current.x, function(x){
              x[(half.n.smoothing.samples + 1):(length(x) - half.n.smoothing.samples)]
            })
          } # if(smooth.loop)
          
          pdf(paste0(save.dir,lock.loop,' - ',task.loop,ifelse(smooth.loop, ' - smoothed', ' - unsmoothed'),'.pdf'),
              width = 9, height = 4.5)
          # par(mfcol = c(1,1),
              # oma = c(0,0,1,0))
          plot.time.series(.y.values = current.y,
                           .x.values = current.x,
                           .title = '',
                           .colors = roi.colors[[tl.split.loop]][plot.regions],
                           #.y.limits = c(0, 1.6),
                           .error.bars = current.error.bar,
                           # .x.limits = x.limits[[lock.loop]][[task.loop]],
                           # .x.ticks = x.ticks[[lock.loop]][[task.loop]],
                           ###
                           ###
                           .y.limits.min.at.least = pub.regions.y.limits['min'],
                           .y.limits.max.at.least = pub.regions.y.limits['max'],
                           .y.label = '',
                           .y.ticks = list('locked_to_production_onset' = list('beta' = c(0,.5,1),
                                                                               'high_gamma' = c(0,1))[[band.loop]],
                                           'locked_to_stimulus_onset' = list('beta' = c(-.25, 0, .25, .5),
                                                                             'high_gamma' = c(0,1))[[band.loop]])[[lock.loop]],
                           .y.tick.labels = list('locked_to_production_onset' = list('beta' = c("0",'.5','1'),
                                                                                     'high_gamma' = c("0","1"))[[band.loop]],
                                                 'locked_to_stimulus_onset' = list('beta' = c('', '0', ' ', '.5'),
                                                                                   'high_gamma' = c("0","1"))[[band.loop]])[[lock.loop]],
                           .x.ticks = list(#'locked_to_production_onset' = seq(-1000, 1000, by = 500),
                             'locked_to_production_onset' = c(-median.rt.times[task.loop], 0, 1000),
                                           'locked_to_stimulus_onset' = seq(0, 2000, by = 500))[[lock.loop]],
                           .x.tick.labels = list(#'locked_to_production_onset' = c('-1000', '', '0', '', '1000'),
                             'locked_to_production_onset' = c(-median.rt.times[task.loop], 0, 1000),
                                                 'locked_to_stimulus_onset' = c('0', '', '1000', '', '2000'))[[lock.loop]],
                           .x.label = '',#time (ms)',
                           .x.limits = list('locked_to_production_onset' = c(-time.convert(median.rt.samples[task.loop], 'samples', 'times'), 
                                                                             time.convert(median.rt.samples['sp'], 'samples', 'times') + 1000),
                                            'locked_to_stimulus_onset' = c(0, 2000))[[lock.loop]],
                           show.t0 = (lock.loop == 'locked_to_production_onset'),
                           .horizontal.line.at = list('high_gamma' = NA,
                                                      'beta' = 0)[[band.loop]],
                           .vertical.line.at = ifelse(lock.loop == 'locked_to_stimulus_onset',time.convert(median.rt.samples[task.loop], 'samples', 'times'), NA),
                           .zoom = 1.8,
                           .margin = c(3.5,3.5,1.5,1.5),
                           .theme = theme.loop,
                           .background = rgb(1,1,1,0))
          
          # add.text.line.multiple.colors(text.segments = plot.regions,
          #                               text.colors = roi.colors[[tl.split.loop]][plot.regions],
          #                               .outer = TRUE)
          dev.off()  
        }#; rm(task.loop)
        
      } # smooth.loop
      
    }#; rm(lock.loop)
  }#; rm(theme.loop)
}#; rm(tl.split.loop)

# } # alpha.label and alpha.loop








# } # band.loop

#stopCluster(cl)

message('Script completed successfully. ',Sys.time())
