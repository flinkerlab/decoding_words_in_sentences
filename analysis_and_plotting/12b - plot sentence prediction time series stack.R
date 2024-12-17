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
library('colorspace')
library('shades')


# Clean up
rm(list=ls())
cat("\014")
message("Begin plotting stacked sentence noun prediction accuracies. ",Sys.time())

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
    n.cores.to.use = 4
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

### Load ROIs and stages
load(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/stages and rois/stages and rois.RData'))

## Load RTs
load(paste0(path, 'analysis/R/downsample data/warped data/high_gamma/median.rt.samples 256 Hz.RData'))
median.rt.times <- time.convert(median.rt.samples, "samples", "times", 256)

## Load colors
load(paste0(output.path,'data/0a - definitions/noun colors/noun colors.RData'))


### Loop thru models
model.types <- list.files(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/'))
for(model.loop in model.types){ ## UNCOMMENT
  # model.loop = model.types[1] ## UNCOMMENT
  
  
  ### Get all significant windows and z-scored accuracies
  combined.stage.windows <- list()
  patient.stage.windows <- list()
  combined.stage.stats <- list()
  patient.stage.stats <- list()
  message(band.loop,' - ',model.loop,' - Attaching significant windows... ',Sys.time())
  for(roi.loop in rois){ ## UNCOMMENT
    # roi.loop = rois[1] ## UNCOMMENT
    
    ### Loop thru stages
    # Storage
    combined.stage.windows[[roi.loop]] <- list()
    patient.stage.windows[[roi.loop]] <- list()
    combined.stage.stats[[roi.loop]] <- list()
    patient.stage.stats[[roi.loop]] <- list()
    
    # For clean up
    stage.keep <- c(ls(), 'stage.keep', 'stage.loop')
    
    # Loop!
    for(stage.loop in rev(names(stage.sample.labels))){ ## UNCOMMENT
      # stage.loop = names(stage.sample.labels)[6] ## UNCOMMENT
      
      # Clean up
      rm(list = ls()[! ls() %in% stage.keep])
      gc()
      
      ### Load data
      read.data.dir <- paste0(output.path, 'data/10b_vE - get sentence stats from permutations - 50ms/',
                              model.loop,'/',
                              roi.loop,'/',
                              stage.loop,'/')
      
      # Load stats: 'patient.sl.stats','patient.sl.sig.windows','combined.sl.stats','combined.sl.sig.windows','half.n.smoothing.samples.sl'
      # read.filename <- list.files(read.data.dir)[grepl("patient and combined sentence stats", list.files(read.data.dir))]
      read.filename <- 'patient and combined sentence stats - 1000 shuffles.RData'
      # if(read.filename != "patient and combined sentence stats - 1000 shuffles.RData"){message("Warning: ",model.loop,' - ',roi.loop,' - ',stage.loop,": Not 1000 shuffles!!!!")}
      attach(paste0(read.data.dir, read.filename)) 
      combined.stage.windows[[roi.loop]][[stage.loop]] <- combined.sl.sig.windows
      patient.stage.windows[[roi.loop]][[stage.loop]] <- patient.sl.sig.windows
      combined.stage.stats[[roi.loop]][[stage.loop]] <- combined.sl.stats
      patient.stage.stats[[roi.loop]][[stage.loop]] <- patient.sl.stats
      detach()
      
      # Clean up
      rm(read.filename)
      
    }; rm(stage.loop)
  }; rm(roi.loop)
  message('...detached! ', Sys.time()) # < 1min
  
  
  ###
  ### Plots
  ###
  
  # Update
  message(band.loop, ' - ', model.loop,' - Beginning plotting! (May take ~5 minutes - inefficient plot code). ',Sys.time())
  
  # Plot grids
  task.grids <- list(
    'active' = c('subject.noun.active','verb.active','object.noun.active'),
    'passive' = c('subject.noun.passive','verb.passive','object.noun.passive'),
    'sentence' = c('subject.noun','verb','object.noun'),
    'list' = c('list1','list2'))
  task.grid.length <- lapply(task.grids, length)
  
  for(task.loop in names(task.grids)){ ## UNCOMMENT
    # task.loop = names(task.grids)[1] ## UNCOMMENT
    
    for(roi.loop in names(combined.stage.windows)){ ## UNCOMMENT
      # roi.loop = names(combined.stage.windows)[1] ## UNCOMMENT
      
      ### Get rois/stages with sig findings
      current.sig.windows <- list()
      for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
        # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
        
        current.sig.windows[[role.loop]] <- list()
        for(noun.loop in c('noun1','noun2')){
          # noun.loop = 'noun1'
          current.sig.windows[[role.loop]][[noun.loop]] <- 
            lapply(combined.stage.windows[[roi.loop]], 
                   function(x){x[[role.loop]][[noun.loop]]})
        }; rm(noun.loop)
      }; rm(role.loop)
      
      ### Get stages with any sig windows
      stages.with.any.sig <- c()
      for(..role in names(current.sig.windows)){
        # ..role = names(current.sig.windows)[1]
        for(..noun in names(current.sig.windows[[..role]])){
          # ..noun = names(current.sig.windows[[..role]])[1]
          stages.with.any.sig <- 
            c(stages.with.any.sig, 
              names(which(sapply(current.sig.windows[[..role]][[..noun]], 
                                 function(x){nrow(x) > 0}))))
        }; rm(..noun)
      }; rm(..role)
      stages.with.any.sig <- unique(stages.with.any.sig)
      stages.with.any.sig.starts <- sapply(strsplit(stages.with.any.sig, "_", fixed = TRUE), function(x){
        as.numeric(gsub("neg","-",x[2],fixed = TRUE))
      })
      stages.with.any.sig <- stages.with.any.sig[order(stages.with.any.sig.starts)]
      current.df <- stage.ranges[stages.with.any.sig,]
      
      ## Throw out all rest
      current.sig.windows <- lapply(current.sig.windows, function(x){
        lapply(x, function(y){y[stages.with.any.sig]})})
      
      ### Get time series
      # Noun1
      ts1.y <- list()
      ts1.x <- list()
      for(role.loop in names(current.sig.windows)){
        ts1.y[[role.loop]] <- 
          lapply(combined.stage.stats[[roi.loop]][rownames(current.df)], function(x){
            x[[role.loop]]$noun1$accuracy - (1/6)})
        ts1.x[[role.loop]] <- 
          lapply(combined.stage.stats[[roi.loop]][rownames(current.df)], function(x){
            x[[role.loop]]$noun1$sample.label})
      }; rm(role.loop)
      
      # Noun2
      ts2.y <- list()
      ts2.x <- list()
      for(role.loop in names(current.sig.windows)){
        ts2.y[[role.loop]] <- 
          lapply(combined.stage.stats[[roi.loop]][rownames(current.df)], function(x){
            x[[role.loop]]$noun2$accuracy - (1/6)})
        ts2.x[[role.loop]] <- 
          lapply(combined.stage.stats[[roi.loop]][rownames(current.df)], function(x){
            x[[role.loop]]$noun2$sample.label})
      }; rm(role.loop)
      
      
      
      
      
      ##
      ## Plot: mean of patients, one cluster, all stages with any sig windows
      ##
      
      # If any data
      if(any(sapply(ts1.y, length) > 0)){
        
        for(theme.loop in c('black','white')){
          # theme.loop = c('black','white')[2]
          
          print(paste0(model.loop,' - ',task.loop,' - ',theme.loop,' - Plotting: /average across patients - one cluster - all stages/. ',Sys.time()))  
          
          for(plot.grid.loop in c(TRUE,FALSE)){
            # plot.grid.loop = TRUE
            
            if(plot.grid.loop){
              save.fig.dir <- 
                paste0(output.path, 
                       'figures/12b_vA - plot sentence prediction time series stack/',
                       'average across patients - one cluster - all stages/',
                       theme.loop,'/',
                       'grids/',
                       task.loop,'/')
              dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
              pdf(paste0(save.fig.dir, roi.loop,' - ',model.loop,' - ',task.loop,' - n_sig=',nrow(current.df),' - 150ms smoothing.pdf'),
                  width = 7 * task.grid.length[[task.loop]], height = 2 + (.18 * length(ts1.y[[1]])))
              par(mfrow = c(1, task.grid.length[[task.loop]]))
            } # if(plot.grid.loop)
            
            for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
              # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
              
              if(! plot.grid.loop){
                save.fig.dir <- 
                  paste0(output.path, 
                         'figures/12b_vA - plot sentence prediction time series stack/',
                         'average across patients - one cluster - all stages/',
                         theme.loop,'/',
                         'individual plots/',
                         task.loop,'/')
                dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
                pdf(paste0(save.fig.dir, roi.loop,' - ',model.loop,' - ',task.loop,' - ',role.loop,' - n_sig=',nrow(current.df),'.pdf'),
                    width = 7, height = 7)
              } # if(! plot.grid.loop)
              
              plot.time.series.stack(.cluster.times.df = current.df,
                                     .sig.windows = current.sig.windows[[role.loop]]$noun1,
                                     .sig.windows2 = current.sig.windows[[role.loop]]$noun2,
                                     .time.series.y = ts1.y[[role.loop]],
                                     .time.series.x = ts1.x[[role.loop]],
                                     .time.series2.y = ts2.y[[role.loop]],
                                     .time.series2.x = ts2.x[[role.loop]],
                                     .sampling.rate = 256,
                                     .x.limits = c(-1000, 500),
                                     .x.ticks = seq(-1000, 500, by = 500),
                                     .title = ifelse(plot.grid.loop, paste0(roi.loop,' - patients combined'), ''),
                                     .y.label = ifelse(plot.grid.loop, 'stages', ''),
                                     .x.label = ifelse(plot.grid.loop, paste0('time from ',role.loop,' onset (ms)'),''),
                                     .sig.colors = noun.colors$active['noun1','hex'],
                                     .sig2.colors = noun.colors$active['noun2','hex'],
                                     # .time.series.color = colorspace::lighten(as.character(shades::saturation(character.colors['nurse','hex'], values = delta(-.23))), 
                                     #                                          amount = .04, method = 'relative'),
                                     # .time.series2.color = colorspace::darken(as.character(shades::saturation(character.colors['dog','hex'], values = delta(-.33))), 
                                     #                                          amount = .04, method = 'relative'),
                                     .time.series.color = noun.colors$passive['noun1','hex'],
                                     .time.series2.color = noun.colors$passive['noun2','hex'],
                                     .train.time.color = 'black',
                                     .omit.clusters.with.no.sig = FALSE,
                                     offset.sig.windows = TRUE,
                                     .row.labels = rownames(current.df),
                                     show.row.labels = FALSE,#plot.grid.loop,
                                     show.sig.bars = FALSE,
                                     show.sig.highlights = TRUE,
                                     show.train.times = (! grepl('verb', role.loop)), # doesn't make sense for verbs
                                     .theme = theme.loop,
                                     .background = rgb(1,1,1,0),
                                     .margin = c(3,3,1.5,1.5),
                                     .cluster.line.width.zoom = .7,
                                     .train.line.width.zoom = .5,
                                     .zoom = 1.8)
              
              # Save
              if(! plot.grid.loop){dev.off()}
              
            }; rm(role.loop)
            
            # Save
            if(plot.grid.loop){dev.off()}
            
          }; rm(plot.grid.loop)
        }#; rm(theme.loop)
      } # if(any(sapply(ts1.y, length) > 0)){
      
      # Clean up
      rm(current.sig.windows, current.df)
      
      
      ##
      ## Plot: individual patients, one cluster, all stages with any sig windows
      ##
      
      print(paste0(model.loop,' - ',task.loop,' - ',theme.loop,' - Plotting: /each patient - one cluster - all stages/. ',Sys.time()))  
      
      ### Get rois/stages with sig findings
      current.sig.windows <- list()
      for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
        # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
        
        current.sig.windows[[role.loop]] <- list()
        for(noun.loop in c('noun1','noun2')){
          # noun.loop = 'noun1'
          current.sig.windows[[role.loop]][[noun.loop]] <- 
            lapply(patient.stage.windows[[roi.loop]], 
                   function(y){lapply(y, function(x){
                     x[[role.loop]][[noun.loop]]})})
          
          current.sig.windows[[role.loop]][[noun.loop]] <- 
            unlist(current.sig.windows[[role.loop]][[noun.loop]], recursive = FALSE)
        }; rm(noun.loop)
      }; rm(role.loop)
      
      ## Get stages with any sig windows
      patient.stages.with.any.sig <- c()
      for(..role in names(current.sig.windows)){
        for(..noun in names(current.sig.windows[[..role]])){
          patient.stages.with.any.sig <- 
            c(patient.stages.with.any.sig, 
              names(which(sapply(current.sig.windows[[..role]][[..noun]], 
                                 function(x){nrow(x) > 0}))))
        }; rm(..noun)
      }; rm(..role)
      patient.stages.with.any.sig <- unique(patient.stages.with.any.sig)
      current.stages <- sapply(patient.stages.with.any.sig, function(x){
        strsplit(x, '.', fixed = TRUE)[[1]][1]})
      if(length(current.stages) == 0){current.stages <- 0}
      current.df <- stage.ranges[current.stages,]
      rownames(current.df) <- patient.stages.with.any.sig
      
      ## Throw out all rest
      current.sig.windows <- lapply(current.sig.windows, function(x){
        lapply(x, function(y){y[patient.stages.with.any.sig]})})
      
      ### Get time series
      # Noun1
      ts1.y <- list()
      ts2.y <- list()
      ts1.x <- list()
      ts2.x <- list()
      for(role.loop in names(current.sig.windows)){
        # role.loop = names(current.sig.windows)[1]
        
        ts1.y[[role.loop]] <- list()
        ts2.y[[role.loop]] <- list()
        ts1.x[[role.loop]] <- list()
        ts2.x[[role.loop]] <- list()
        
        for(deet.loop in rownames(current.df)){
          # deet.loop = rownames(current.df)[1]
          ..deets <- strsplit(deet.loop, '.', fixed = TRUE)[[1]]
          ..stage <- ..deets[1]
          ..patient <- ..deets[2]
          ts1.y[[role.loop]][[deet.loop]] <- 
            patient.stage.stats[[roi.loop]][[..stage]][[..patient]][[role.loop]]$noun1$accuracy - (1/6)
          ts1.x[[role.loop]][[deet.loop]] <- 
            patient.stage.stats[[roi.loop]][[..stage]][[..patient]][[role.loop]]$noun1$sample.label
          ts2.y[[role.loop]][[deet.loop]] <- 
            patient.stage.stats[[roi.loop]][[..stage]][[..patient]][[role.loop]]$noun2$accuracy - (1/6)
          ts2.x[[role.loop]][[deet.loop]] <- 
            patient.stage.stats[[roi.loop]][[..stage]][[..patient]][[role.loop]]$noun2$sample.label
          rm(..deets, ..stage, ..patient)
        }; rm(deet.loop)
        
      }; rm(role.loop)
      
      
      ### Plot!
      for(theme.loop in c('black','white')){
        # theme.loop = c('black','white')[2]
        
        for(plot.grid.loop in c(TRUE,FALSE)){
          # plot.grid.loop = c(TRUE, FALSE)[1]
          
          if(plot.grid.loop){
            save.fig.dir <- 
              paste0(output.path, 
                     'figures/12b_vA - plot sentence prediction time series stack/',
                     'each patient - one cluster - all stages/',
                     theme.loop,'/',
                     'grids/',
                     task.loop,'/')
            dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
            pdf(paste0(save.fig.dir, roi.loop,' - ',model.loop,' - ',task.loop,' - n_sig=',nrow(current.df),' - 150ms smoothing.pdf'),
                width = 7 * task.grid.length[[task.loop]], height = 7)
            par(mfrow = c(1, task.grid.length[[task.loop]]))
          } # if(plot.grid.loop)
          
          for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
            # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
            
            if(! plot.grid.loop){
              save.fig.dir <- 
                paste0(output.path, 
                       'figures/12b_vA - plot sentence prediction time series stack/',
                       'each patient - one cluster - all stages/',
                       theme.loop,'/',
                       'individual plots/',
                       task.loop,'/')
              dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
              pdf(paste0(save.fig.dir, roi.loop,' - ',model.loop,' - ',task.loop,' - ',role.loop,' - n_sig=',nrow(current.df),'.pdf'),
                  width = 7, height = 7)
            } # if(! plot.grid.loop)
            
            # patient.stage.times <- stage.ranges[gsub("\\..*","",names(current.sig.windows$noun1)),]
            # rownames(patient.stage.times) <- names(current.sig.windows$noun1)
            
            plot.time.series.stack(.cluster.times.df = current.df,
                                   .sig.windows = current.sig.windows[[role.loop]]$noun1,
                                   .sig.windows2 = current.sig.windows[[role.loop]]$noun2,
                                   .time.series.y = ts1.y[[role.loop]],
                                   .time.series.x = ts1.x[[role.loop]],
                                   .time.series2.y = ts2.y[[role.loop]],
                                   .time.series2.x = ts2.x[[role.loop]],
                                   .sampling.rate = 256,
                                   .title = ifelse(plot.grid.loop, paste0(roi.loop,' - patients'), ''),
                                   .y.label = ifelse(plot.grid.loop, 'patient-stages', ''),
                                   .x.label = ifelse(plot.grid.loop, paste0('time from ',role.loop,' onset (ms)'),''),
                                   .sig.colors = noun.colors$active['noun1','hex'],
                                   .sig2.colors = noun.colors$active['noun2','hex'],
                                   # .time.series.color = colorspace::lighten(as.character(shades::saturation(character.colors['nurse','hex'], values = delta(-.23))), 
                                   #                                          amount = .04, method = 'relative'),
                                   # .time.series2.color = colorspace::darken(as.character(shades::saturation(character.colors['dog','hex'], values = delta(-.33))), 
                                   #                                          amount = .04, method = 'relative'),
                                   .time.series.color = noun.colors$passive['noun1','hex'],
                                   .time.series2.color = noun.colors$passive['noun2','hex'],
                                   .train.time.color = 'black',
                                   .omit.clusters.with.no.sig = FALSE,
                                   offset.sig.windows = TRUE,
                                   .row.labels = rownames(current.df),
                                   show.row.labels = TRUE,#plot.grid.loop,
                                   show.sig.bars = FALSE,
                                   show.sig.highlights = TRUE,
                                   show.train.times = (! grepl('verb', role.loop)), # doesn't make sense for verbs
                                   .theme = theme.loop,
                                   .background = rgb(1,1,1,0),
                                   .margin = c(3,3,1.5,1.5),
                                   .cluster.line.width.zoom = .7,
                                   .train.line.width.zoom = .5,
                                   .zoom = 1.8)
            
            # Save
            if(! plot.grid.loop){dev.off()}
            
          }; rm(role.loop)
          
          # Save
          if(plot.grid.loop){dev.off()}
          
        }; rm(plot.grid.loop)
      }#; rm(theme.loop)
      
      # Clean up
      rm(current.sig.windows, current.df)
      
      
    }#; rm(roi.loop)  
    
    
    ###
    ### "All rois" plots
    ###
    
    ##
    ## Plot: combined patients, ALL rois, all stages with any sig windows
    ##
    
    
    print(paste0(model.loop,' - ',task.loop,' - Plotting: /average across patients - all rois - all stages/. ',Sys.time()))  
    
    ### Get rois/stages with sig findings
    current.sig.windows <- list()
    for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
      # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
      
      current.sig.windows[[role.loop]] <- list()
      for(noun.loop in c('noun1','noun2')){ # noun.loop = 'noun1'
        current.sig.windows[[role.loop]][[noun.loop]] <- 
          unlist(lapply(combined.stage.windows, function(y){
            lapply(y, function(x){
              x[[role.loop]][[noun.loop]]})}),
            recursive = FALSE)
      }; rm(noun.loop)
    }; rm(role.loop)
    
    ## Get stages with any sig windows
    cluster.stages.with.any.sig <- c()
    for(..role in names(current.sig.windows)){
      for(..noun in names(current.sig.windows[[..role]])){
        cluster.stages.with.any.sig <- 
          c(cluster.stages.with.any.sig, 
            names(which(sapply(current.sig.windows[[..role]][[..noun]], 
                               function(x){nrow(x) > 0}))))
      }; rm(..noun)
    }; rm(..role)
    cluster.stages.with.any.sig <- unique(cluster.stages.with.any.sig)
    current.df <- list()
    for(i in cluster.stages.with.any.sig){
      ..cluster <- strsplit(i, '.', fixed = TRUE)[[1]][1]
      ..stage <- strsplit(i, '.', fixed = TRUE)[[1]][2]
      ..df <- stage.ranges[..stage,]
      rownames(..df) <- paste0(..cluster,'.',rownames(..df))
      current.df[[length(current.df) + 1]] <- ..df
      rm(..cluster, ..stage, ..df)
    }; rm(i)
    current.df <- bind_rows(current.df)
    
    ## Throw out all rest
    current.sig.windows <- lapply(current.sig.windows, function(x){
      lapply(x, function(y){y[cluster.stages.with.any.sig]})})
    
    ### Get time series
    # Noun1
    ts1.y <- list()
    ts1.x <- list()
    ts2.y <- list()
    ts2.x <- list()
    for(role.loop in names(current.sig.windows)){ # role.loop = names(current.sig.windows)[1]
      
      ts1.y[[role.loop]] <- list()
      ts1.x[[role.loop]] <- list()
      ts2.y[[role.loop]] <- list()
      ts2.x[[role.loop]] <- list()
      for(deet.loop in rownames(current.df)){ # deet.loop = rownames(current.df)[1]
        ..deets <- strsplit(deet.loop, '.', fixed = TRUE)[[1]]
        ..cluster <- ..deets[1]
        ..stage <- ..deets[2]
        
        ts1.y[[role.loop]][[deet.loop]] <- combined.stage.stats[[..cluster]][[..stage]][[role.loop]]$noun1$accuracy - (1/6)
        ts1.x[[role.loop]][[deet.loop]] <- combined.stage.stats[[..cluster]][[..stage]][[role.loop]]$noun1$sample.label
        ts2.y[[role.loop]][[deet.loop]] <- combined.stage.stats[[..cluster]][[..stage]][[role.loop]]$noun2$accuracy - (1/6)
        ts2.x[[role.loop]][[deet.loop]] <- combined.stage.stats[[..cluster]][[..stage]][[role.loop]]$noun2$sample.label
        
        rm(..deets, ..cluster, ..stage)
      }; rm(deet.loop)
    }; rm(role.loop)
    
    
    ### Plot!
    for(theme.loop in c('black','white')){
      # theme.loop = c('black','white')[2]
      
      for(plot.grid.loop in c(TRUE,FALSE)){
        # plot.grid.loop = c(TRUE,FALSE)[1]
        
        if(plot.grid.loop){
          save.fig.dir <- 
            paste0(output.path, 
                   'figures/12b_vA - plot sentence prediction time series stack/',
                   'average across patients - all rois - all stages/',
                   theme.loop,'/',
                   'grids/',
                   task.loop,'/')
          dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
          pdf(paste0(save.fig.dir, model.loop,' - ',task.loop,' - n_sig=',nrow(current.df),' - 150ms smoothing.pdf'),
              width = 7 * task.grid.length[[task.loop]], height = 7)
          par(mfrow = c(1, task.grid.length[[task.loop]]))
        } # if(plot.grid.loop)
        
        for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
          # role.loop = task.grids[[task.loop]][2] ## UNCOMMENT
          
          if(! plot.grid.loop){
            save.fig.dir <- 
              paste0(output.path, 
                     'figures/12b_vA - plot sentence prediction time series stack/',
                     'average across patients - all rois - all stages/',
                     theme.loop,'/',
                     'individual plots/',
                     task.loop,'/')
            dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
            pdf(paste0(save.fig.dir, model.loop,' - ',task.loop,' - ',role.loop,' - n_sig=',nrow(current.df),'.pdf'),
                width = 7, height = 7)
          } # if(! plot.grid.loop)
          
          plot.time.series.stack(.cluster.times.df = current.df,
                                 .sig.windows = current.sig.windows[[role.loop]]$noun1,
                                 .sig.windows2 = current.sig.windows[[role.loop]]$noun2,
                                 .time.series.y = ts1.y[[role.loop]],
                                 .time.series.x = ts1.x[[role.loop]],
                                 .time.series2.y = ts2.y[[role.loop]],
                                 .time.series2.x = ts2.x[[role.loop]],
                                 .sampling.rate = 256,
                                 .title = ifelse(plot.grid.loop, 'patient means', ''),
                                 .y.label = ifelse(plot.grid.loop, 'cluster-stages', ''),
                                 .x.label = ifelse(plot.grid.loop, paste0('time from ',role.loop,' onset (ms)'),''),
                                 .sig.colors = noun.colors$active['noun1','hex'],
                                 .sig2.colors = noun.colors$active['noun2','hex'],
                                 # .time.series.color = colorspace::lighten(as.character(shades::saturation(character.colors['nurse','hex'], values = delta(-.23))), 
                                 #                                          amount = .04, method = 'relative'),
                                 # .time.series2.color = colorspace::darken(as.character(shades::saturation(character.colors['dog','hex'], values = delta(-.33))), 
                                 #                                          amount = .04, method = 'relative'),
                                 .time.series.color = noun.colors$passive['noun1','hex'],
                                 .time.series2.color = noun.colors$passive['noun2','hex'],
                                 .train.time.color = 'black',
                                 .omit.clusters.with.no.sig = FALSE,
                                 offset.sig.windows = TRUE,
                                 .row.labels = rownames(current.df),
                                 show.row.labels = FALSE,#plot.grid.loop,
                                 show.sig.bars = FALSE,
                                 show.sig.highlights = TRUE,
                                 show.train.times = (! grepl('verb', role.loop)), # doesn't make sense for verbs
                                 .theme = theme.loop,
                                 .background = rgb(1,1,1,0),
                                 .margin = c(3,3,1.5,1.5),
                                 .cluster.line.width.zoom = .7,
                                 .train.line.width.zoom = .5,
                                 .zoom = 1.8)
          
          # Save
          if(! plot.grid.loop){dev.off()}
          
        }; rm(role.loop)
        
        # Save
        if(plot.grid.loop){dev.off()}
        
      }; rm(plot.grid.loop)
    }#; rm(theme.loop)
    
    # Clean up
    rm(current.sig.windows, current.df)
    
    
    ##
    ## Plot: individual patients, ALL rois, all stages with any sig windows
    ##
    
    print(paste0(model.loop,' - ',task.loop,' - ',theme.loop,' - Plotting: /each patient - all rois - all stages/. ',Sys.time()))  
    
    ### Get rois/stages with sig findings
    current.sig.windows <- list()
    for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
      # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
      
      current.sig.windows[[role.loop]] <- list()
      for(noun.loop in c('noun1','noun2')){ # noun.loop = 'noun1'
        current.sig.windows[[role.loop]][[noun.loop]] <- 
          unlist(lapply(patient.stage.windows, function(z){
            unlist(lapply(z, function(y){
              lapply(y, function(x){
                x[[role.loop]][[noun.loop]]})}), 
              recursive = FALSE)}), 
            recursive = FALSE)
      }; rm(noun.loop)
    }; rm(role.loop)
    
    ## Get stages with any sig windows
    cluster.stage.patients.with.any.sig <- c()
    for(..role in names(current.sig.windows)){
      for(..noun in names(current.sig.windows[[..role]])){
        cluster.stage.patients.with.any.sig <- 
          c(cluster.stage.patients.with.any.sig, 
            names(which(sapply(current.sig.windows[[..role]][[..noun]], 
                               function(x){nrow(x) > 0}))))
      }; rm(..noun)
    }; rm(..role)
    cluster.stage.patients.with.any.sig <- unique(cluster.stage.patients.with.any.sig)
    
    ## Get dataframe for plot function to link times to stages
    current.df <- list()
    for(i in cluster.stage.patients.with.any.sig){
      ..deets <- strsplit(i, '.', fixed = TRUE)[[1]]
      ..cluster <- ..deets[1]
      ..stage <- ..deets[2]
      ..patient <- ..deets[3]
      ..df <- stage.ranges[..stage,]
      rownames(..df) <- i
      current.df[[length(current.df) + 1]] <- ..df
      rm(..deets, ..cluster, ..stage, ..df)
    }; rm(i)
    current.df <- bind_rows(current.df)
    
    ## Throw out all rest
    current.sig.windows <- lapply(current.sig.windows, function(x){
      lapply(x, function(y){y[cluster.stage.patients.with.any.sig]})})
    
    ### Get time series
    # Noun1
    ts1.y <- list()
    ts2.y <- list()
    ts1.x <- list()
    ts2.x <- list()
    for(role.loop in names(current.sig.windows)){
      # role.loop = names(current.sig.windows)[1]
      
      ts1.y[[role.loop]] <- list()
      ts2.y[[role.loop]] <- list()
      ts1.x[[role.loop]] <- list()
      ts2.x[[role.loop]] <- list()
      
      for(deet.loop in rownames(current.df)){
        # deet.loop = rownames(current.df)[1]
        ..deets <- strsplit(deet.loop, '.', fixed = TRUE)[[1]]
        ..cluster <- ..deets[1]
        ..stage <- ..deets[2]
        ..patient <- ..deets[3]
        ts1.y[[role.loop]][[deet.loop]] <- 
          patient.stage.stats[[..cluster]][[..stage]][[..patient]][[role.loop]]$noun1$accuracy - (1/6)
        ts1.x[[role.loop]][[deet.loop]] <- 
          patient.stage.stats[[..cluster]][[..stage]][[..patient]][[role.loop]]$noun1$sample.label
        ts2.y[[role.loop]][[deet.loop]] <- 
          patient.stage.stats[[..cluster]][[..stage]][[..patient]][[role.loop]]$noun2$accuracy - (1/6)
        ts2.x[[role.loop]][[deet.loop]] <- 
          patient.stage.stats[[..cluster]][[..stage]][[..patient]][[role.loop]]$noun2$sample.label
        rm(..deets, ..stage, ..patient)
      }; rm(deet.loop)
      
    }; rm(role.loop)
    
    ### Plot!
    for(theme.loop in c('black','white')){
      # theme.loop = c('black','white')[2]
      
      for(plot.grid.loop in c(TRUE,FALSE)){
        # plot.grid.loop = c(TRUE, FALSE)[1]
        
        if(plot.grid.loop){
          save.fig.dir <- 
            paste0(output.path, 
                   'figures/12b_vA - plot sentence prediction time series stack/',
                   'each patient - all rois - all stages/',
                   theme.loop,'/',
                   'grids/',
                   task.loop,'/')
          dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
          pdf(paste0(save.fig.dir, model.loop,' - ',task.loop,' - n_sig=',nrow(current.df),' - 150ms smoothing.pdf'),
              width = 7 * task.grid.length[[task.loop]], height = 10)
          par(mfrow = c(1, task.grid.length[[task.loop]]))
        } # if(plot.grid.loop)
        
        for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
          # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
          
          if(! plot.grid.loop){
            save.fig.dir <- 
              paste0(output.path, 
                     'figures/12b_vA - plot sentence prediction time series stack/',
                     'each patient - all rois - all stages/',
                     theme.loop,'/',
                     'individual plots/',
                     task.loop,'/')
            dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
            pdf(paste0(save.fig.dir, model.loop,' - ',task.loop,' - ',role.loop,' - n_sig=',nrow(current.df),'.pdf'),
                width = 7, height = 8)
          } # if(! plot.grid.loop)
          
          plot.time.series.stack(.cluster.times.df = current.df,
                                 .sig.windows = current.sig.windows[[role.loop]]$noun1,
                                 .sig.windows2 = current.sig.windows[[role.loop]]$noun2,
                                 .time.series.y = ts1.y[[role.loop]],
                                 .time.series.x = ts1.x[[role.loop]],
                                 .time.series2.y = ts2.y[[role.loop]],
                                 .time.series2.x = ts2.x[[role.loop]],
                                 .time.series.scale.factor = 3,
                                 .sampling.rate = 256,
                                 .title = ifelse(plot.grid.loop, 'patients', ''),
                                 .y.label = ifelse(plot.grid.loop, 'patient-cluster-stages', ''),
                                 .x.label = ifelse(plot.grid.loop, paste0('time from ',role.loop,' onset (ms)'),''),
                                 .sig.colors = noun.colors$active['noun1','hex'],
                                 .sig2.colors = noun.colors$active['noun2','hex'],
                                 # .time.series.color = noun.colors$passive['noun1','hex'],
                                 # .time.series2.color = noun.colors$passive['noun2','hex'],
                                 # .time.series.color = colorspace::lighten(as.character(shades::saturation(character.colors['nurse','hex'], values = delta(-.23))), 
                                 #                                          amount = .04, method = 'relative'),
                                 # .time.series2.color = colorspace::darken(as.character(shades::saturation(character.colors['dog','hex'], values = delta(-.33))), 
                                 #                                          amount = .04, method = 'relative'),
                                 .time.series.color = rgb(.6,.6,.6),
                                 .time.series2.color = rgb(.6,.6,.6),
                                 .train.time.color = 'black',
                                 .omit.clusters.with.no.sig = FALSE,
                                 offset.sig.windows = TRUE,
                                 .row.labels = rownames(current.df),
                                 show.row.labels = FALSE,#plot.grid.loop,
                                 show.sig.bars = FALSE,
                                 show.sig.highlights = TRUE,
                                 show.train.times = (! grepl('verb', role.loop)), # doesn't make sense for verbs
                                 .theme = theme.loop,
                                 .density.height.proportion = .1,
                                 .background = rgb(1,1,1,0),
                                 .margin = c(3,3,1.5,1.5),
                                 .cluster.line.width.zoom = .9,
                                 .train.line.width.zoom = .7,
                                 plot.guidelines = FALSE,
                                 .x.ticks = c(-1000, 0, 500),
                                 .zoom = 1.8)
          
          # Save
          if(! plot.grid.loop){dev.off()}
          
        }; rm(role.loop)
        
        # Save
        if(plot.grid.loop){dev.off()}
        
      }#; rm(plot.grid.loop)
    }#; rm(theme.loop)
    
  }#; rm(task.loop)
  
  
  
  
  
}#; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


