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
# Get sig windows inverse
source(paste0(path,'/analysis/R/functions/get_significant_windows_inverse.R'))
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
  
  
  ###
  ### Combined data (collapsed across patients)
  ###
  
  for(task.loop in names(task.grids)){ ## UNCOMMENT
    # task.loop = names(task.grids)[1] ## UNCOMMENT
    
    # Storage for combined data
    x.vals <-
      y.vals <-
      sig.vals.samples <-
      max.val.for.scaling <- list()
    
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
      
      
      ##
      ## Calculate sum and density of sig windows
      ##
      
      x.vals[[roi.loop]] <-
        y.vals[[roi.loop]] <-
        sig.vals.samples[[roi.loop]] <-
        max.val.for.scaling[[roi.loop]] <- list()
      x.vals.samples <- time.convert(-1000, 'times', 'samples', 256):time.convert(500, 'times', 'samples', 256)
      x.vals.times <- time.convert(x.vals.samples, 'samples', 'times', 256)
      
      for(metric.loop in c('density','sum')){
        # metric.loop = c('density','sum')[1]
        
        x.vals[[roi.loop]][[metric.loop]] <-
          y.vals[[roi.loop]][[metric.loop]] <-
          sig.vals.samples[[roi.loop]][[metric.loop]] <-
          max.val.for.scaling[[roi.loop]][[metric.loop]] <- list()
        
        for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
          # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
          
          # Storage
          x.vals[[roi.loop]][[metric.loop]][[role.loop]] <-
            y.vals[[roi.loop]][[metric.loop]][[role.loop]] <-
            sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]] <-
            max.val.for.scaling[[roi.loop]][[metric.loop]][[role.loop]] <- list()
          
          for(noun.loop in c('noun1','noun2')){
            # noun.loop = c('noun1','noun2')[1]
            
            sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <-
              lapply(current.sig.windows[[role.loop]][[noun.loop]], function(.window){
                .samples.sig <- get.significant.windows.inverse(
                  .window,
                  .samples = x.vals.samples,
                  .sampling.rate = 256, 
                  output.class = "vector")
                .samples.sig <- names(which(.samples.sig == 1))
                .samples.sig <- time.convert(.samples.sig, "sample.labels", "samples")
                return(.samples.sig)
              })
            if(length(sapply(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]], 
                             function(x){length(x) > 1})) == 0){
              current.n.sig.classifiers <- 0
            }else{
              current.n.sig.classifiers <- sum(sapply(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]], 
                                                      function(x){length(x) > 1}))
            }
            
            if(metric.loop == 'density'){
              if(current.n.sig.classifiers < 2){
                
                message("Fewer than 2 significant samples in .sig.windows (TS1), not enough to do density plot; setting density vals to 0.")
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- rep(0, times = length(x.vals.times) + 2)
                x.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(-1200, x.vals.times, 700)
                max.val.for.scaling[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- 1
                
              }else{
                
                density <- density(time.convert(unlist(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]]),
                                                'samples', 'times', 256), bw = 'nrd0')
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(0, density$y, 0)
                x.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(-1200, density$x, 700)
                max.val.for.scaling[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- 
                  max(table(unlist(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]])))
                
              }
            } # if(metric.loop == 'density')
            
            if(metric.loop == 'sum'){
              
              x.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(-1200, x.vals.times, 700)
              
              if(current.n.sig.classifiers == 0){
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- rep(0, times = length(x.vals.times) + 2)
              }
              if(current.n.sig.classifiers == 1){
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- rep(0, times = length(x.vals.times))
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]][x.vals.samples %in% sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]]] <- 1
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(0, y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]], 0)
              }
              if(current.n.sig.classifiers > 1){
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- 
                  lapply(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]], function(.window){
                    .output <- rep(0, times = length(x.vals.times))
                    .output[x.vals.samples %in% .window] <- 1
                    .output <- c(0, .output, 0)
                    return(.output)
                  })
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- 
                  rowSums(do.call(cbind, y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]]))
              }
              
            } # if(metric.loop == 'sum')
            
            # Clean up
            rm(current.n.sig.classifiers)
            
          }; rm(noun.loop)
        }#; rm(role.loop)
      }#; rm(metric.loop)
    }#; rm(roi.loop)
    
    ### Calculate y limits per ROI
    .y.max <- list()
    for(metric.loop in c('density','sum')){
      
      .y.max[[metric.loop]] <- sapply(y.vals, function(.roi){
        max(unlist(.roi[[metric.loop]]))})
      
    }; rm(metric.loop)
    
    
    
    ##
    ## Plot: mean of patients, one cluster, all stages with any sig windows
    ##
    
    # for(theme.loop in c('black','white')){
    theme.loop = c('black','white')[2]
    
    print(paste0(model.loop,' - ',task.loop,' - ',theme.loop,' - Plotting: /average across patients - one cluster - all stages/. ',Sys.time()))  
    
    for(roi.loop in names(combined.stage.windows)){
      # roi.loop = names(combined.stage.windows)[1]
      
      for(metric.loop in c('density','sum')){
        # metric.loop = c('density','sum')[2]
        
        for(plot.grid.loop in c(TRUE,FALSE)){
          # plot.grid.loop = TRUE
          
          if(plot.grid.loop){
            save.fig.dir <- 
              paste0(output.path, 
                     'figures/21a - plot significant prediction densities and sums by ROI/',
                     'average across patients - one cluster - all stages/',
                     theme.loop,'/',
                     'grids/',
                     task.loop,'/',
                     metric.loop,'/')
            dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
            pdf(paste0(save.fig.dir, model.loop,' - ',roi.loop,' - ',task.loop,'.pdf'),
                width = 7 * task.grid.length[[task.loop]], height = 5)
            par(mfrow = c(1, task.grid.length[[task.loop]]))
          } # if(plot.grid.loop)
          
          for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
            # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
            
            if(! plot.grid.loop){
              save.fig.dir <- 
                paste0(output.path, 
                       'figures/21a - plot significant prediction densities and sums by ROI/',
                       'average across patients - one cluster - all stages/',
                       theme.loop,'/',
                       'individual plots/',
                       task.loop,'/',
                       metric.loop,'/')
              dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
              pdf(paste0(save.fig.dir, model.loop,' - ',roi.loop,' - ',task.loop,' - ',role.loop,'.pdf'),
                  width = 7, height = 5)
            } # if(! plot.grid.loop)
            
            plot.time.series(.y.values = y.vals[[roi.loop]][[metric.loop]][[role.loop]],
                             .x.values = x.vals[[roi.loop]][[metric.loop]][[role.loop]],
                             .sampling.rate = 256,
                             .x.limits = c(-1000, 500),
                             .y.limits = c(0, max(.y.max[[metric.loop]])),
                             .y.ticks = c(0, max(unlist(y.vals))),
                             .x.ticks = seq(-1000, 500, by = 500),
                             .colors = noun.colors$active[c('noun1','noun2'),'hex'],
                             .y.label = '',
                             .x.label = '',
                             .title = ifelse(plot.grid.loop, role.loop, ''),
                             .zoom = 1.8,
                             .theme = theme.loop)
            
            # Save
            if(! plot.grid.loop){dev.off()}
            
          }; rm(role.loop)
          
          # Save
          if(plot.grid.loop){dev.off()}
          
        }; rm(plot.grid.loop)
        
      }#; rm(metric.loop)
      
    }; rm(roi.loop)
      
  # }; rm(theme.loop)
  }; rm(task.loop)
    
    
  ###
  ### Patient data
  ###
  
  for(task.loop in names(task.grids)){ ## UNCOMMENT
    # task.loop = names(task.grids)[1] ## UNCOMMENT
    
    # Storage for combined data
    x.vals <-
      y.vals <-
      sig.vals.samples <-
      max.val.for.scaling <- list()
    
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
      
      
      ##
      ## Calculate sum and density of sig windows
      ##
      
      x.vals[[roi.loop]] <-
        y.vals[[roi.loop]] <-
        sig.vals.samples[[roi.loop]] <-
        max.val.for.scaling[[roi.loop]] <- list()
      x.vals.samples <- time.convert(-1000, 'times', 'samples', 256):time.convert(500, 'times', 'samples', 256)
      x.vals.times <- time.convert(x.vals.samples, 'samples', 'times', 256)
      
      for(metric.loop in c('density','sum')){
        # metric.loop = c('density','sum')[1]
        
        x.vals[[roi.loop]][[metric.loop]] <-
          y.vals[[roi.loop]][[metric.loop]] <-
          sig.vals.samples[[roi.loop]][[metric.loop]] <-
          max.val.for.scaling[[roi.loop]][[metric.loop]] <- list()
        
        for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
          # role.loop = task.grids[[task.loop]][1] ## UNCOMMENT
          
          # Storage
          x.vals[[roi.loop]][[metric.loop]][[role.loop]] <-
            y.vals[[roi.loop]][[metric.loop]][[role.loop]] <-
            sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]] <-
            max.val.for.scaling[[roi.loop]][[metric.loop]][[role.loop]] <- list()
          
          for(noun.loop in c('noun1','noun2')){
            # noun.loop = c('noun1','noun2')[1]
            
            sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <-
              lapply(current.sig.windows[[role.loop]][[noun.loop]], function(.window){
                .samples.sig <- get.significant.windows.inverse(
                  .window,
                  .samples = x.vals.samples,
                  .sampling.rate = 256, 
                  output.class = "vector")
                .samples.sig <- names(which(.samples.sig == 1))
                .samples.sig <- time.convert(.samples.sig, "sample.labels", "samples")
                return(.samples.sig)
              })
            if(length(sapply(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]], 
                             function(x){length(x) > 1})) == 0){
              current.n.sig.classifiers <- 0
            }else{
              current.n.sig.classifiers <- sum(sapply(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]], 
                                                      function(x){length(x) > 1}))
            }
            
            if(metric.loop == 'density'){
              if(current.n.sig.classifiers < 2){
                
                message("Fewer than 2 significant samples in .sig.windows (TS1), not enough to do density plot; setting density vals to 0.")
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- rep(0, times = length(x.vals.times) + 2)
                x.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(-1200, x.vals.times, 700)
                max.val.for.scaling[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- 1
                
              }else{
                
                density <- density(time.convert(unlist(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]]),
                                                'samples', 'times', 256), bw = 'nrd0')
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(0, density$y, 0)
                x.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(-1200, density$x, 700)
                max.val.for.scaling[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- 
                  max(table(unlist(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]])))
                
              }
            } # if(metric.loop == 'density')
            
            if(metric.loop == 'sum'){
              
              x.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(-1200, x.vals.times, 700)
              
              if(current.n.sig.classifiers == 0){
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- rep(0, times = length(x.vals.times) + 2)
              }
              if(current.n.sig.classifiers == 1){
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- rep(0, times = length(x.vals.times))
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]][x.vals.samples %in% sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]]] <- 1
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- c(0, y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]], 0)
              }
              if(current.n.sig.classifiers > 1){
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- 
                  lapply(sig.vals.samples[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]], function(.window){
                    .output <- rep(0, times = length(x.vals.times))
                    .output[x.vals.samples %in% .window] <- 1
                    .output <- c(0, .output, 0)
                    return(.output)
                  })
                y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] <- 
                  rowSums(do.call(cbind, y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]]))
              }
              
            } # if(metric.loop == 'sum')
            
            # Clean up
            rm(current.n.sig.classifiers)
            
          }; rm(noun.loop)
        }#; rm(role.loop)
      }#; rm(metric.loop)
    }#; rm(roi.loop)
    
    ### Calculate y limits per ROI
    .y.max <- list()
    for(metric.loop in c('density','sum')){
      
      .y.max[[metric.loop]] <- sapply(y.vals, function(.roi){
        max(unlist(.roi[[metric.loop]]))})
      
    }; rm(metric.loop)
    
    ### Calculate y limits per ROI per noun per role
    .y.max.per.role <- list()
    for(metric.loop in c('density','sum')){
      
      .y.max.per.role[[metric.loop]] <- lapply(y.vals, function(.roi){
        sapply(.roi[[metric.loop]], function(.role){
          sapply(.role, max)})
        })
      
    }; rm(metric.loop)
    
    
    
    ##
    ## Plot: individual patients, one cluster, all stages with any sig windows
    ##
    
    # for(theme.loop in c('black','white')){
    theme.loop = c('black','white')[2]
    
    print(paste0(model.loop,' - ',roi.loop,' - ',task.loop,' - ',theme.loop,' - Plotting: /each patient - one cluster - all stages/. ',Sys.time()))  
    
    for(roi.loop in names(combined.stage.windows)){
      # roi.loop = names(combined.stage.windows)[3]
      
      for(metric.loop in c('density','sum')){
        # metric.loop = c('density','sum')[1]
        
        for(plot.grid.loop in c(TRUE,FALSE)){
          # plot.grid.loop = TRUE
          
          if(plot.grid.loop){
            save.fig.dir <- 
              paste0(output.path, 
                     'figures/21a - plot significant prediction densities and sums by ROI/',
                     'each patient - one cluster - all stages/',
                     theme.loop,'/',
                     'grids/',
                     task.loop,'/',
                     metric.loop,'/')
            dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
            pdf(paste0(save.fig.dir, model.loop,' - ',roi.loop,' - ',task.loop,'.pdf'),
                width = 7 * task.grid.length[[task.loop]], height = 5)
            par(mfrow = c(1, task.grid.length[[task.loop]]))
          } # if(plot.grid.loop)
          
          for(role.loop in task.grids[[task.loop]]){ ## UNCOMMENT
            # role.loop = task.grids[[task.loop]][3] ## UNCOMMENT
            
            if(! plot.grid.loop){
              save.fig.dir <- 
                paste0(output.path, 
                       'figures/21a - plot significant prediction densities and sums by ROI/',
                       'each patient - one cluster - all stages/',
                       theme.loop,'/',
                       'individual plots/',
                       task.loop,'/',
                       metric.loop,'/')
              dir.create(save.fig.dir, showWarnings = FALSE, recursive = TRUE)
              pdf(paste0(save.fig.dir, model.loop,' - ',roi.loop,' - ',task.loop,' - ',role.loop,'.pdf'),
                  width = 7, height = 5)
            } # if(! plot.grid.loop)
            
            # Scale y-value densities
            min.nonzero.density.val <- 
              min(unlist(.y.max.per.role$density)[unlist(.y.max.per.role$density) != 0])
            if(metric.loop == 'density'){
              current.y <- list()
              for(noun.loop in c('noun1','noun2')){
                current.y[[noun.loop]] <-
                  y.vals[[roi.loop]][[metric.loop]][[role.loop]][[noun.loop]] / 
                  max(min.nonzero.density.val, .y.max.per.role$density[[roi.loop]][noun.loop, role.loop]) *
                  .y.max.per.role$sum[[roi.loop]][noun.loop, role.loop]
              }; rm(noun.loop)
            }else{ # if(metric.loop == 'density')
              current.y <- y.vals[[roi.loop]][[metric.loop]][[role.loop]]
            }
            
            plot.time.series(.y.values = current.y,
                             .x.values = x.vals[[roi.loop]][[metric.loop]][[role.loop]],
                             .sampling.rate = 256,
                             .x.limits = c(-1000, 500),
                             .y.limits = c(0, ifelse(roi.loop == 'SMC', 30, 10)), # ifelse(metric.loop == 'sum', 30, max(.y.max[[metric.loop]]))),
                             .y.ticks = c(0, ifelse(roi.loop == 'SMC', 30, 10)),
                             # .y.ticks = c(0, ceiling(max(.y.max$sum) / 5) * 5), # ifelse(metric.loop == 'sum', ceiling(.y.max[[metric.loop]] / 5) * 5, max(.y.max[[metric.loop]]))),
                             .x.ticks = c(-1000, 0, 500),
                             .colors = noun.colors$active[c('noun1','noun2'),'hex'],
                             .y.label = '',
                             .x.label = '',
                             .title = ifelse(plot.grid.loop, role.loop, ''),
                             .background = rgb(1,1,1,0),
                             .zoom = 1.8,
                             .theme = theme.loop)
            
            # Save
            if(! plot.grid.loop){dev.off()}
            
          }; rm(role.loop)
          
          # Save
          if(plot.grid.loop){dev.off()}
          
        }; rm(plot.grid.loop)
        
      }#; rm(metric.loop)
      
    }; rm(roi.loop)
      # }#; rm(theme.loop)
    
  }; rm(task.loop)
  
  
  
  
  
  
  
}#; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


