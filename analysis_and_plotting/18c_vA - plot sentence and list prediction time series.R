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

### Load ROIs and stages
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
# for(model.loop in model.types){ ## UNCOMMENT
model.loop = model.types[1] ## UNCOMMENT

### Loop thru rois
for(roi.loop in rois){ ## UNCOMMENT
  # roi.loop = rois[3] ## UNCOMMENT
  
  ### Loop thru stages
  # For clean up
  stage.keep <- c(ls(), 'stage.keep', 'stage.loop')
  
  # Loop!
  for(stage.loop in names(stage.sample.labels)){ ## UNCOMMENT
    # stage.loop = names(stage.sample.labels)[17] ## UNCOMMENT
    
    # Clean up
    rm(list = ls()[! ls() %in% stage.keep])
    gc()
    
    # Load SVO, list1, list2 stats ("combined.sl.sig.windows","combined.sl.stats","half.n.smoothing.samples.sl","patient.sl.sig.windows","patient.sl.stats")
    load(paste0(output.path, 
                'data/10b_vE - get sentence stats from permutations - 50ms/',
                model.loop,'/',
                roi.loop,'/',
                stage.loop,'/patient and combined sentence stats - 1000 shuffles.RData'))
    
    # Load non-sentence-initial function word stats
    load(paste0(output.path, 
                'data/13b_vA - get stats - non-sentence-initial function words/',
                model.loop,'/',
                roi.loop,'/',
                stage.loop,'/patient and combined sentence stats - 1000 shuffles.RData'))
    
    # Load non-sentence-initial function word stats
    load(paste0(output.path, 
                'data/14b_vA - get stats - sentence-initial function words/',
                model.loop,'/',
                roi.loop,'/',
                stage.loop,'/patient and combined sentence stats - 1000 shuffles.RData'))
    
    
    ### Merge patient stats
    patient.stats <- list()
    patient.sig.windows <- list()
    for(patient in names(patient.sl.stats)){
      # patient = names(patient.sl.stats)[1]
      
      patient.stats[[patient]] <- list()
      patient.sig.windows[[patient]] <- list()
      
      # SL stats
      for(role.loop in names(patient.sl.stats[[patient]])){
        # role.loop = names(patient.sl.stats[[patient]])[1]
        patient.stats[[patient]][[role.loop]] <- patient.sl.stats[[patient]][[role.loop]]
        patient.sig.windows[[patient]][[role.loop]] <- patient.sl.sig.windows[[patient]][[role.loop]]
      }; rm(role.loop)
      
      # FW.NS1 stats (First word, non-sentence-initial)
      for(role.loop in names(patient.fw.ns1.stats[[patient]])){
        # role.loop = names(patient.fw.ns1.stats[[patient]])[1]
        patient.stats[[patient]][[role.loop]] <- patient.fw.ns1.stats[[patient]][[role.loop]]
        patient.sig.windows[[patient]][[role.loop]] <- patient.fw.ns1.sig.windows[[patient]][[role.loop]]
      }; rm(role.loop)
      
      # FW.S1 stats (First word, sentence-initial)
      for(role.loop in names(patient.fw.s1.stats[[patient]])){
        # role.loop = names(patient.fw.s1.stats[[patient]])[1]
        patient.stats[[patient]][[role.loop]] <- patient.fw.s1.stats[[patient]][[role.loop]]
        patient.sig.windows[[patient]][[role.loop]] <- patient.fw.s1.sig.windows[[patient]][[role.loop]]
      }; rm(role.loop)
      
    }; rm(patient)
    
    
    ### Merge combined stats
    combined.stats <- list()
    combined.sig.windows <- list()
    role.sample.labels <- list()
    
    # SL stats
    for(role.loop in names(combined.sl.stats)){
      # role.loop = names(patient.sl.stats[[patient]])[1]
      combined.stats[[role.loop]] <- combined.sl.stats[[role.loop]]
      combined.sig.windows[[role.loop]] <- combined.sl.sig.windows[[role.loop]]
      role.sample.labels[[role.loop]] <- combined.sl.stats[[role.loop]][[1]]$sample.label
    }; rm(role.loop)
    
    # FW.NS1 stats (First word, non-sentence-initial)
    for(role.loop in names(combined.fw.ns1.stats)){
      # role.loop = names(combined.fw.ns1.stats)[1]
      combined.stats[[role.loop]] <- combined.fw.ns1.stats[[role.loop]]
      combined.sig.windows[[role.loop]] <- combined.fw.ns1.sig.windows[[role.loop]]
      role.sample.labels[[role.loop]] <- combined.fw.ns1.stats[[role.loop]][[1]]$sample.label
    }; rm(role.loop)
    
    # FW.S1 stats (First word, sentence-initial)
    for(role.loop in names(combined.fw.s1.stats)){
      # role.loop = names(combined.fw.s1.stats)[1]
      combined.stats[[role.loop]] <- combined.fw.s1.stats[[role.loop]]
      combined.sig.windows[[role.loop]] <- combined.fw.s1.sig.windows[[role.loop]]
      role.sample.labels[[role.loop]] <- combined.fw.s1.stats[[role.loop]][[1]]$sample.label
    }; rm(role.loop)
    
    
    # Ranges
    role.time.ranges <- lapply(role.sample.labels, function(x){
      range(time.convert(x, 'sample.labels', 'times', 256))
    })
    
    
    ##
    ## Metadata
    ##
    
    save.plot.dir.main <- paste0(output.path, 'figures/18c_vA - plot sentence and list prediction time series - 50ms/')
    zoom <- 1.8
    
    # What data to plot: sentences? lists?
    lx.comparisons <- list(
      # 'sentences' = c('subject.det',
      #                 'subject.noun',
      #                 'aux',
      #                 'verb',
      #                 'object.det',
      #                 'object.noun'), # from left to right in grid
      'actives' = c(#'subject.det.active',
                    'subject.noun.active',
                    'aux.active',
                    'verb.active',
                    # 'object.det.active',
                    'object.noun.active'),
      'passives' = c(#'subject.det.passive',
                     'subject.noun.passive',
                     'aux.passive',
                     'passive.be',
                     'verb.passive',
                     'passive.by',
                     # 'object.det.passive',
                     'object.noun.passive'),
      'lists' = c('list1',
                  'list2')
    )
    
    
    
    # Plots will appear in rows from first word decoded to last, so make the first one in each row "special" (have y-axis labels, show earlier time points)
    # leftmost.grid.plots <- unlist(sapply(lx.comparisons, function(x){x[grepl('subject',x)]}))
    # leftmost.grid.plots <- c(leftmost.grid.plots, 'list1', 'list2')
    
    plot.parameters <- list()
    plot.parameters[['subject.noun.active']] <- list('time.limits' = c(-1100, 500),
                                                     'x.ticks' = c(-1000, 0, 500),
                                                     'x.tick.labels' = c(-1000, 0, 500),
                                                     'y.axis' = TRUE)
    plot.parameters[['aux.active']] <- list('time.limits' = c(-500, 500),
                                            'x.ticks' = c(-500, 0, 500),
                                            'x.tick.labels' = c("-500", "", "500"),
                                            'y.axis' = FALSE)
    plot.parameters[['verb.active']] <-  list('time.limits' = c(-500, 500),
                                              'x.ticks' = c(-500, 0, 500),
                                              'x.tick.labels' = c(-500, 0, 500),
                                              'y.axis' = FALSE)
    plot.parameters[['object.noun.active']] <- plot.parameters$verb.active
    plot.parameters[['subject.noun.passive']] <- plot.parameters$subject.noun.active
    plot.parameters[['aux.passive']] <- plot.parameters$aux.active
    plot.parameters[['passive.be']] <- plot.parameters$aux.active
    plot.parameters[['verb.passive']] <- plot.parameters$verb.active
    plot.parameters[['passive.by']] <- plot.parameters$aux.active
    plot.parameters[['object.noun.passive']] <- plot.parameters$object.noun.active
    plot.parameters[['list1']] <- plot.parameters$subject.noun.active
    plot.parameters[['list2']] <- list('time.limits' = c(-1000, 500),
                                       'x.ticks' = c(-1000, 0, 500),
                                       'x.tick.labels' = c(-1000, 0, 500),
                                       'y.axis' = FALSE)
    
    
    
    
    ##
    ## Plot: Means of patients
    ##
    
    for(theme.loop in c('white','black')[1]){
      # theme.loop = c('white','black')[1]
      
      for(comparison.loop in names(lx.comparisons)){ ## UNCOMMENT
        # comparison.loop = names(lx.comparisons)[2] ## UNCOMMENT
        
        # Set up storage
        ys <- list()
        xs <- list()
        sig.windows <- list()
        shuffle.mean <- list()
        shuffle.error.lower <- list()
        shuffle.error.upper <- list()
        
        # Get plot data
        for(lx.loop in lx.comparisons[[comparison.loop]]){
          # lx.loop = lx.comparisons[[comparison.loop]][1]
          
          ys[[lx.loop]] <- lapply(combined.stats[[lx.loop]], function(x){x$accuracy})
          xs[[lx.loop]] <- lapply(combined.stats[[lx.loop]], function(x){x$sample.label})
          sig.windows[[lx.loop]] <- combined.sig.windows[[lx.loop]]
          shuffle.mean[[lx.loop]] <- combined.stats[[lx.loop]]$noun1$shuffle.mean
          shuffle.error.lower[[lx.loop]] <- 
            combined.stats[[lx.loop]]$noun1$shuffle.mean - 
            combined.stats[[lx.loop]]$noun1$shuffle.95ci.lower
          shuffle.error.upper[[lx.loop]] <- 
            combined.stats[[lx.loop]]$noun1$shuffle.95ci.upper - 
            combined.stats[[lx.loop]]$noun1$shuffle.mean
        }; rm(lx.loop)
        
        # Y limits
        all.y.vals <- c(unlist(ys),
                        unlist(shuffle.mean) + unlist(shuffle.error.upper),
                        unlist(shuffle.mean) - unlist(shuffle.error.lower))
        all.y.vals <- all.y.vals[! is.na(all.y.vals)]
        y.limits <- range(all.y.vals)
        y.limits[1] <- floor(100 * y.limits[1] / 2) / 100 * 2
        y.limits[2] <- ceiling(100 * y.limits[2] / 2) / 100 * 2
        y.ticks <- y.limits
        
        # For manuscript plot manually set this
        if(roi.loop == 'SMC' & stage.loop == 'window_50_to_100ms'){
          y.limits <- c(.09, .3)
          y.ticks <- c(.1, .3)
        }
        
        
        ### Grids and individual plots
        for(grid.loop in c(TRUE, FALSE)){
          # grid.loop = c(TRUE, FALSE)[1]
          
          # Save dir
          save.means.dir <- paste0(save.plot.dir.main, 
                                   'mean accuracy/',
                                   theme.loop,'/',
                                   ifelse(grid.loop, 'grids/', 'individual plots for publication/'),
                                   model.loop,'/',
                                   comparison.loop,'/')
          dir.create(save.means.dir, showWarnings = FALSE, recursive = TRUE)
          
          # Separate grids with significant findings from those with none
          current.sig <- any(as.vector(sapply(sig.windows[lx.comparisons[[comparison.loop]]],
                                              function(x){sapply(x, function(y){nrow(y) > 0})})))
          
          # Just don't plot if not sig
          if(current.sig){
            
            # Save grid plot
            if(grid.loop){
              
              current.save.means.sig.dir <- 
                paste0(save.means.dir,
                       ifelse(current.sig,'significant/', 'not significant/'))
              dir.create(current.save.means.sig.dir, showWarnings = FALSE, recursive = TRUE)
              
              pdf(paste0(current.save.means.sig.dir,
                         comparison.loop,' - ',
                         roi.loop,' - ',
                         stage.loop,'.pdf'),
                  width = 7.5 + (5 * length(lx.comparisons[[comparison.loop]])), 
                  height = 4.5)
              par(mfrow = c(1, length(lx.comparisons[[comparison.loop]])))
            } # if(! grid.loop){
            
            ### Plot
            for(lx.loop in lx.comparisons[[comparison.loop]]){
              # lx.loop = lx.comparisons[[comparison.loop]][1]
              
              # if(grid.loop){
              # x.limits <- c(-1000, 500)
              # }else{
              # x.limits <- c(ifelse(lx.loop %in% leftmost.grid.plots, -1000, -500), 500)
              # }
              x.limits <- c(-1000, 500)
              
              # plot.x.lims <- c(ifelse(lx.loop %in% leftmost.grid.plots, -1100, -600), 600)
              plot.x.lims <- plot.parameters[[lx.loop]]$time.limits
              plot.x.lims.sample.labels <- time.convert(
                time.convert(plot.x.lims[1] - 50, 'times', 'samples', 256):
                  time.convert(plot.x.lims[2] + 50, 'times', 'samples', 256),
                'samples', 'sample.labels', 256)
              plot.indices <- which(xs[[lx.loop]][['noun1']] %in% plot.x.lims.sample.labels)
              current.ys <- lapply(ys[[lx.loop]], function(.noun){.noun[plot.indices]})
              current.xs <- lapply(xs[[lx.loop]], function(.noun){.noun[plot.indices]})
              
              # Save individual plot
              if(! grid.loop){
                pdf(paste0(save.means.dir, 
                           comparison.loop,' - ',
                           roi.loop,' - ',
                           stage.loop,' - ',
                           lx.loop,'.pdf'),
                    width = 7, height = 4.5)
              } # if(! grid.loop){
              
              plot.time.series(.y.values = current.ys,
                               .x.values = current.xs,
                               .sampling.rate = 256,
                               # .colors = noun.colors$active[c('noun1','noun2'),'hex'],
                               # .colors = rep(rgb(0,0,0), times = 2),
                               .colors = colorspace::darken(noun.colors$active[c('noun1','noun2'),'hex'], amount = .45, method = "relative"),
                               .y.lwd = c(rep(3, times = length(patient.stats)), 6),
                               .y.limits = y.limits,
                               .y.ticks = y.ticks,
                               .y.tick.labels = c(ifelse(plot.parameters[[lx.loop]]$y.axis, y.ticks[1], ''),
                                                  ifelse(plot.parameters[[lx.loop]]$y.axis, y.ticks[2], '')),
                               .y.label = '',
                               .x.label = '',
                               .x.limits = x.limits,
                               # .x.ticks = seq(x.limits[1], x.limits[2], by = 500),
                               .x.ticks = plot.parameters[[lx.loop]]$x.ticks,
                               .x.tick.labels = plot.parameters[[lx.loop]]$x.tick.labels,
                               .shuffle.dist.mean = shuffle.mean[[lx.loop]][plot.indices],
                               .shuffle.dist.error.bars.lower = shuffle.error.lower[[lx.loop]][plot.indices],
                               .shuffle.dist.error.bars.upper = shuffle.error.upper[[lx.loop]][plot.indices],
                               .shuffle.dist.color = 'grey',
                               .polygons.x = unlist(stage.ranges[stage.loop,c('start.time','end.time')]),
                               .polygons.color = ifelse(theme.loop == 'white', 'black', 'white'),
                               # show.polygons = ! grepl('verb',lx.loop),
                               show.polygons = ! (grepl('verb',lx.loop) | grepl('aux',lx.loop) | grepl('.by',lx.loop,fixed = TRUE) | grepl('.be',lx.loop,fixed = TRUE) | grepl('.det',lx.loop,fixed = TRUE)),
                               .sig.windows = sig.windows[[lx.loop]],
                               # .sig.color = colors[[theme.loop]]$rainbow_bright['purple','hex'],
                               # .sig.color = colors[[theme.loop]]$rainbow['pink','hex'],#ifelse(theme.loop == 'white', 'black', 'white'),##noun.colors$passive[c('noun1','noun2'),'hex'],
                               .sig.color = noun.colors$active[c('noun1','noun2'),'hex'],#colors[[theme.loop]]$rainbow['pink','hex'],#ifelse(theme.loop == 'white', 'black', 'white'),##noun.colors$passive[c('noun1','noun2'),'hex'],
                               show.y.axis = plot.parameters[[lx.loop]]$y.axis,
                               show.sig.bars = FALSE,
                               show.sig.highlights = TRUE,
                               # .sig.highlights.in.front = TRUE,
                               .sig.highlights.lwd = 3 * 3 * zoom,
                               .center.sig.bars.vertically = FALSE,
                               .title = '',
                               .theme = theme.loop,
                               .background = rgb(1,1,1,0),
                               .margin = c(3,3,1.5,1.5),
                               .zoom = zoom)
              
              # Save
              if(! grid.loop){dev.off()}
              
            }; rm(lx.loop)
            
            # Save
            if(grid.loop){dev.off()}
            
          } # if(current.sig)
        }; rm(grid.loop)
      }; rm(comparison.loop)
      
      
      
      ##
      ## Individual patients
      ##
      
      for(patient in names(patient.stats)){ ## UNCOMMENT
        # patient = names(patient.stats)[1] ## UNCOMMENT
        
        for(comparison.loop in names(lx.comparisons)){ ## UNCOMMENT
          # comparison.loop = names(lx.comparisons)[4] ## UNCOMMENT
          
          # Set up storage
          ys <- list()
          xs <- list()
          sig.windows <- list()
          shuffle.mean <- list()
          shuffle.error.lower <- list()
          shuffle.error.upper <- list()
          
          # Get plot data
          for(lx.loop in lx.comparisons[[comparison.loop]]){
            # lx.loop = lx.comparisons[[comparison.loop]][1]
            
            ys[[lx.loop]] <- lapply(patient.stats[[patient]][[lx.loop]], function(x){x$accuracy})
            xs[[lx.loop]] <- lapply(patient.stats[[patient]][[lx.loop]], function(x){x$sample.label})
            sig.windows[[lx.loop]] <- patient.sig.windows[[patient]][[lx.loop]]
            shuffle.mean[[lx.loop]] <- patient.stats[[patient]][[lx.loop]]$noun1$shuffle.mean
            shuffle.error.lower[[lx.loop]] <- 
              shuffle.mean[[lx.loop]] - 
              patient.stats[[patient]][[lx.loop]]$noun1$shuffle.95ci.lower
            shuffle.error.upper[[lx.loop]] <- 
              patient.stats[[patient]][[lx.loop]]$noun1$shuffle.95ci.upper - 
              shuffle.mean[[lx.loop]]
          }; rm(lx.loop)
          
          # Y limits
          all.y.vals <- c(unlist(ys),
                          unlist(shuffle.mean) + unlist(shuffle.error.upper),
                          unlist(shuffle.mean) - unlist(shuffle.error.lower))
          all.y.vals <- all.y.vals[! is.na(all.y.vals)]
          y.limits <- range(all.y.vals)
          y.limits[1] <- 0#floor(100 * y.limits[1] / 2) / 100 * 2
          y.limits[2] <- ceiling(100 * y.limits[2] / 5) / 100 * 5
            # ceiling(100 * y.limits[2] / 2) / 100 * 2
          y.ticks <- y.limits
          
          
          ### Grids and individual plots
          for(grid.loop in c(TRUE, FALSE)){
            # grid.loop = c(TRUE, FALSE)[2]
            
            # Save dir
            save.patients.dir <-
              paste0(save.plot.dir.main,
                     'individual patients/',
                     theme.loop,'/',
                     ifelse(grid.loop, 'grids/', 'individual plots for publication/'),
                     model.loop,'/',
                     comparison.loop,'/')
            dir.create(save.patients.dir, showWarnings = FALSE, recursive = TRUE)
            
            # Sort grids into folders by whether any significant results this patient/roi/stage
            current.sig <- any(as.vector(sapply(sig.windows[lx.comparisons[[comparison.loop]]],
                                                function(x){sapply(x, function(y){nrow(y) > 0})})))
            
            # Just don't plot if not sig
            if(current.sig){
              
              # Save grid plot
              if(grid.loop){
                current.save.patients.sig.dir <- 
                  paste0(save.patients.dir,
                         ifelse(current.sig,'significant/', 'not significant/'))
                dir.create(current.save.patients.sig.dir, showWarnings = FALSE, recursive = TRUE)
                
                # Begin figure save
                pdf(paste0(current.save.patients.sig.dir,
                           comparison.loop,' - ',
                           roi.loop,' - ',
                           patient,' - ',
                           stage.loop,'.pdf'),
                    width = 6 + (4 * (length(lx.comparisons[[comparison.loop]]) - 1)), 
                    height = 4.5)
                par(mfrow = c(1, length(lx.comparisons[[comparison.loop]])))
              } # if(! grid.loop){
              
              ### Plot
              for(lx.loop in lx.comparisons[[comparison.loop]]){
                # lx.loop = lx.comparisons[[comparison.loop]][1]
                
                # Skip patients with no data
                if(! any(sapply(ys[[lx.loop]], function(x){all(is.na(x))}))){
                  
                  # X limits
                  # if(grid.loop){
                    # x.limits <- c(-500, 500)
                  # }else{
                    # x.limits <- c(ifelse(lx.loop %in% leftmost.grid.plots, -1000, -500), 500)
                  # }
                  x.limits <- c(-1000, 500)
                  
                  # plot.x.lims <- c(ifelse(lx.loop %in% leftmost.grid.plots, -1100, -600), 600)
                  plot.x.lims <- plot.parameters[[lx.loop]]$time.limits
                  plot.x.lims.sample.labels <- time.convert(
                    time.convert(plot.x.lims[1] - 50, 'times', 'samples', 256):
                    time.convert(plot.x.lims[2] + 50, 'times', 'samples', 256),
                    'samples', 'sample.labels', 256)
                  plot.indices <- which(xs[[lx.loop]][['noun1']] %in% plot.x.lims.sample.labels)
                  current.ys <- lapply(ys[[lx.loop]], function(.noun){.noun[plot.indices]})
                  current.xs <- lapply(xs[[lx.loop]], function(.noun){.noun[plot.indices]})
                  
                  # Save individual plot
                  if(! grid.loop){
                    pdf(paste0(save.patients.dir, 
                               comparison.loop,' - ',
                               roi.loop,' - ',
                               patient,' - ',
                               stage.loop,' - ',
                               lx.loop,'.pdf'),
                        width = 7, 
                        height = 4.5)
                  } # if(! grid.loop){
                  
                  plot.time.series(.y.values = current.ys,
                                   .x.values = current.xs,
                                   .sampling.rate = 256,
                                   # .colors = noun.colors$active[c('noun1','noun2'),'hex'],
                                   # .colors = rep(rgb(0,0,0), times = 2),
                                   .colors = colorspace::darken(noun.colors$active[c('noun1','noun2'),'hex'], amount = .45, method = "relative"),
                                   .y.lwd = c(rep(3, times = length(patient.stats)), 6),
                                   .y.limits = y.limits,
                                   .y.ticks = c(y.limits[1], y.limits[2]),
                                   .y.tick.labels = c(ifelse(plot.parameters[[lx.loop]]$y.axis, y.ticks[1], ''),
                                                      ifelse(plot.parameters[[lx.loop]]$y.axis, y.ticks[2], '')),
                                   .y.label = '',
                                   .x.label = '',
                                   .x.limits = x.limits,
                                   # .x.ticks = unique(c(ifelse(lx.loop %in% leftmost.grid.plots, -1000, -500), -500, 0, 500)),
                                   .x.ticks = plot.parameters[[lx.loop]]$x.ticks,
                                   .x.tick.labels = plot.parameters[[lx.loop]]$x.tick.labels,
                                   .shuffle.dist.mean = shuffle.mean[[lx.loop]][plot.indices],
                                   .shuffle.dist.error.bars.lower = shuffle.error.lower[[lx.loop]][plot.indices],
                                   .shuffle.dist.error.bars.upper = shuffle.error.upper[[lx.loop]][plot.indices],
                                   .shuffle.dist.color = 'grey',
                                   .polygons.x = unlist(stage.ranges[stage.loop,c('start.time','end.time')]),
                                   .polygons.color = ifelse(theme.loop == 'white', 'black', 'white'),
                                   show.polygons = (! (grepl('verb',lx.loop) | grepl('aux',lx.loop) | grepl('.by',lx.loop,fixed = TRUE) | grepl('.be',lx.loop,fixed = TRUE) | grepl('.det',lx.loop,fixed = TRUE))),
                                   # show.polygons = FALSE,
                                   # .highlight.x.axis.segments = unlist(stage.ranges[stage.loop,c('start.time','end.time')]),
                                   # .highlight.x.axis.segments.color = ifelse((! (grepl('verb',lx.loop) | grepl('aux',lx.loop) | grepl('.by',lx.loop,fixed = TRUE) | grepl('.be',lx.loop,fixed = TRUE) | grepl('.det',lx.loop,fixed = TRUE))), 'black', rgb(1,1,1,0)),
                                   # show.y.axis = (lx.loop %in% leftmost.grid.plots) & (lx.loop != 'list2'),
                                   show.y.axis = plot.parameters[[lx.loop]]$y.axis,
                                   .sig.windows = sig.windows[[lx.loop]],
                                   .sig.color = noun.colors$active[c('noun1','noun2'),'hex'],#colors[[theme.loop]]$rainbow['pink','hex'],#ifelse(theme.loop == 'white', 'black', 'white'),##noun.colors$passive[c('noun1','noun2'),'hex'],
                                   show.sig.bars = FALSE,
                                   show.sig.highlights = TRUE,
                                   # .sig.highlights.in.front = TRUE,
                                   .sig.highlights.lwd = 3 * 3 * zoom,
                                   .center.sig.bars.vertically = FALSE,
                                   .title = '',
                                   .theme = theme.loop,
                                   .background = rgb(1,1,1,0),
                                   .margin = c(3,3,1.5,1.5),
                                   .zoom = zoom)
                  
                  # Save
                  if(! grid.loop){dev.off()}
                  
                } # Skip patients with no data
              }; rm(lx.loop)
              
              # Save
              if(grid.loop){dev.off()}
              
            } # if(current.sig)
            
          }; rm(grid.loop)
          
        }#; rm(comparison.loop)
        
      }#; rm(patient)
      
    }#; rm(theme.loop)
    
    
    # Progress update
    message('Completed: ',
            model.loop,' - ',
            band.loop,' - ',
            roi.loop,' - ',
            stage.loop,'.',
            Sys.time())
    
  }; rm(stage.loop)
}; rm(roi.loop)


# }; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


