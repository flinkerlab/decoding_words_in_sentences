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
output.path <- paste0(path, 'analysis/R/track words during sentences/output/',band.loop,'/')

### Load elec info
load(paste0(output.path, 'data/0a - definitions/elec info/elec info.RData'))

### Load data
# loads "median.rt.samples","patient.trial.info","sampling.rate","warped.data","warped.sample.labels"
# load(paste0(path,'analysis/R/downsample data/warped data/high_gamma/warped high_gamma data at 256 Hz.RData')) 
# summary(patient.trial.info$NY869$sp$case)



### Load clusters
nmf.window.loop <- "50ms_post_stim-to-200ms_post_prod"
distance.metric <- "bayesian.t.test"
winning.rank <- list('high_gamma' = 9, 'beta' = 6)[[band.loop]]
cluster.assignments <- 
  read.csv(paste0(output.path,
                  'data/3a_vD - cluster RDMs - stack sig RDMs/data for brain plots/RDMs - ',
                  distance.metric,
                  '/stage=',nmf.window.loop,
                  '/cluster_assignments_rank=',winning.rank,'.csv'))
rownames(cluster.assignments) <- cluster.assignments$elec
cluster.elecs <- split(cluster.assignments$elec, cluster.assignments$cluster)
clusters <- names(cluster.elecs)
clusters <- clusters[order(as.numeric(gsub("NMF_","",clusters)))]


### Load change points
load(paste0(output.path, 'data/7a - multivariate change point detection - 50ms/breakpoints.RData'))


### Load RTs
load(paste0(path, 'analysis/R/downsample data/warped data/high_gamma/median.rt.samples 256 Hz.RData'))
median.rt.times <- time.convert(median.rt.samples, "samples", "times", 256)
rm(median.rt.samples) # risky to keep since dealing with 256 and 512 Hz datasets


### Load colors
load(paste0(output.path,'data/0a - definitions/noun colors/noun colors.RData'))
load(paste0(path,'analysis/R/color palettes/output/all palettes.RData'))


### Loop thru models
model.types <- list.files(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/'))
for(model.loop in model.types){ ## UNCOMMENT
  # model.loop = model.types[1] ## UNCOMMENT
  
  ### Loop thru clusters
  for(cluster.loop in rev(clusters)){ ## UNCOMMENT
    # cluster.loop = rev(clusters[8]) ## UNCOMMENT
    
    ### Loop thru stages
    # For clean up
    stage.keep <- c(ls(), 'stage.keep', 'stage.loop')
    
    # Loop!
    for(stage.loop in rev(names(stage.sample.labels[[cluster.loop]]))){ ## UNCOMMENT
      # stage.loop = names(stage.sample.labels[[cluster.loop]])[11] ## UNCOMMENT
      
      # Clean up
      rm(list = ls()[! ls() %in% stage.keep])
      gc()
      
      # Load stats ("combined.sl.sig.windows","combined.sl.stats","half.n.smoothing.samples.sl","patient.sl.sig.windows","patient.sl.stats")
      load(paste0(output.path, 
                  'data/10b_vE - get sentence stats from permutations - 50ms/',
                  model.loop,'/',
                  cluster.loop,'/',
                  stage.loop,'/patient and combined sentence stats - 1000 shuffles.RData'))
      
      
      ##
      ## Metadata
      ##
      
      save.plot.dir.main <- paste0(output.path, 'figures/10c - plot sentence and list prediction time series - 50ms/')
      zoom <- 1.8
      
      # What data to plot: sentences? lists?
      lx.comparisons <- list(
        'sentences' = c('subject.noun','verb','object.noun'), # from left to right in grid
        'actives' = c('subject.noun.active','verb.active','object.noun.active'),
        'passives' = c('subject.noun.passive','verb.passive','object.noun.passive'),
        'lists' = c('list1','list2')
      )
      
      # Plots will appear in rows from first word decoded to last, so make the first one in each row "special" (have y-axis labels, show earlier time points)
      leftmost.grid.plots <- sapply(lx.comparisons, function(x){x[1]})
      
      
      ##
      ## Plot: Means of patients
      ##
      
      for(theme.loop in c('white','black')[1]){
        # theme.loop = c('white','black')[1]
        
        for(comparison.loop in names(lx.comparisons)){ ## UNCOMMENT
          # comparison.loop = names(lx.comparisons)[1] ## UNCOMMENT
          
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
            
            ys[[lx.loop]] <- lapply(combined.sl.stats[[lx.loop]], function(x){x$accuracy})
            xs[[lx.loop]] <- lapply(combined.sl.stats[[lx.loop]], function(x){x$sample.label})
            sig.windows[[lx.loop]] <- combined.sl.sig.windows[[lx.loop]]
            shuffle.mean[[lx.loop]] <- combined.sl.stats[[lx.loop]]$noun1$shuffle.mean
            shuffle.error.lower[[lx.loop]] <- 
              combined.sl.stats[[lx.loop]]$noun1$shuffle.mean - 
              combined.sl.stats[[lx.loop]]$noun1$shuffle.95ci.lower
            shuffle.error.upper[[lx.loop]] <- 
              combined.sl.stats[[lx.loop]]$noun1$shuffle.95ci.upper - 
              combined.sl.stats[[lx.loop]]$noun1$shuffle.mean
          }; rm(lx.loop)
          
          # Y limits
          all.y.vals <- c(unlist(ys),
                          unlist(shuffle.mean) + unlist(shuffle.error.upper),
                          unlist(shuffle.mean) - unlist(shuffle.error.lower))
          y.limits <- range(all.y.vals)
          y.limits[1] <- floor(100 * y.limits[1] / 2) / 100 * 2
          y.limits[2] <- ceiling(100 * y.limits[2] / 2) / 100 * 2
          
          
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
            
            # Save grid plot
            if(grid.loop){
              
              # Separate grids with significant findings from those with none
              current.sig <- any(as.vector(sapply(sig.windows[lx.comparisons[[comparison.loop]]],
                                                  function(x){sapply(x, function(y){nrow(y) > 0})})))
              current.save.means.sig.dir <- 
                paste0(save.means.dir,
                       ifelse(current.sig,'significant/', 'not significant/'))
              dir.create(current.save.means.sig.dir, showWarnings = FALSE, recursive = TRUE)
              
              pdf(paste0(current.save.means.sig.dir,
                         comparison.loop,' - ',
                         cluster.loop,' - ',
                         stage.loop,'_vE.pdf'),
                  width = 7.5 + (5 * length(lx.comparisons[[comparison.loop]])), 
                  height = 4.5)
              par(mfrow = c(1, length(lx.comparisons[[comparison.loop]])))
            } # if(! grid.loop){
            
            ### Plot
            for(lx.loop in lx.comparisons[[comparison.loop]]){
              # lx.loop = lx.comparisons[[comparison.loop]][1]
              
              # if(grid.loop){
              x.limits <- c(-1000, 500)
              # }else{
              # x.limits <- c(ifelse(lx.loop %in% leftmost.grid.plots, -1000, -500), 500)
              # }
              
              # Save individual plot
              if(! grid.loop){
                pdf(paste0(save.means.dir, 
                           comparison.loop,' - ',
                           cluster.loop,' - ',
                           stage.loop,' - ',
                           lx.loop,'_vE.pdf'),
                    width = 7, height = 4.5)
              } # if(! grid.loop){
              
              plot.time.series(.y.values = ys[[lx.loop]],
                               .x.values = xs[[lx.loop]],
                               .sampling.rate = 256,
                               .colors = noun.colors$active[c('noun1','noun2'),'hex'],
                               .y.lwd = c(rep(3, times = length(patient.sl.stats)), 6),
                               .y.limits = y.limits,
                               .y.ticks = c(y.limits[1], y.limits[2]),
                               .y.tick.labels = c(ifelse(lx.loop %in% leftmost.grid.plots, y.limits[1], ''),
                                                  ifelse(lx.loop %in% leftmost.grid.plots, y.limits[2], '')),
                               .y.label = '',
                               .x.label = '',
                               .x.limits = x.limits,
                               .x.ticks = seq(x.limits[1], x.limits[2], by = 500),
                               .shuffle.dist.mean = shuffle.mean[[lx.loop]],
                               .shuffle.dist.error.bars.lower = shuffle.error.lower[[lx.loop]],
                               .shuffle.dist.error.bars.upper = shuffle.error.upper[[lx.loop]],
                               .shuffle.dist.color = 'grey',
                               .polygons.x = stage.times.range[[cluster.loop]][[stage.loop]],
                               .polygons.color = ifelse(theme.loop == 'white', 'black', 'white'),
                               show.polygons = ! grepl('verb',lx.loop),
                               .sig.windows = sig.windows[[lx.loop]],
                               .sig.color = colors[[theme.loop]]$rainbow_bright['pink','hex'],
                               show.sig.bars = FALSE,
                               show.sig.highlights = TRUE,
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
            
          }; rm(grid.loop)
        }; rm(comparison.loop)
        
        
        
        ##
        ## Individual patients
        ##
        
        for(patient in names(patient.sl.stats)){ ## UNCOMMENT
          # patient = names(patient.sl.stats)[1] ## UNCOMMENT
          
          for(comparison.loop in names(lx.comparisons)){ ## UNCOMMENT
            # comparison.loop = names(lx.comparisons)[1] ## UNCOMMENT
            
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
              
              ys[[lx.loop]] <- lapply(patient.sl.stats[[patient]][[lx.loop]], function(x){x$accuracy})
              xs[[lx.loop]] <- lapply(patient.sl.stats[[patient]][[lx.loop]], function(x){x$sample.label})
              sig.windows[[lx.loop]] <- patient.sl.sig.windows[[patient]][[lx.loop]]
              shuffle.mean[[lx.loop]] <- patient.sl.stats[[patient]][[lx.loop]]$noun1$shuffle.mean
              shuffle.error.lower[[lx.loop]] <- 
                shuffle.mean[[lx.loop]] - 
                patient.sl.stats[[patient]][[lx.loop]]$noun1$shuffle.95ci.lower
              shuffle.error.upper[[lx.loop]] <- 
                patient.sl.stats[[patient]][[lx.loop]]$noun1$shuffle.95ci.upper - 
                shuffle.mean[[lx.loop]]
            }; rm(lx.loop)
            
            # Y limits
            all.y.vals <- c(unlist(ys),
                            unlist(shuffle.mean) + unlist(shuffle.error.upper),
                            unlist(shuffle.mean) - unlist(shuffle.error.lower))
            y.limits <- range(all.y.vals)
            y.limits[1] <- floor(100 * y.limits[1] / 2) / 100 * 2
            y.limits[2] <- ceiling(100 * y.limits[2] / 2) / 100 * 2
            
            
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
              
              # Save grid plot
              if(grid.loop){
                # Sort grids into folders by whether any significant results this patient/cluster/stage
                current.sig <- any(as.vector(sapply(sig.windows[lx.comparisons[[comparison.loop]]],
                                                    function(x){sapply(x, function(y){nrow(y) > 0})})))
                current.save.patients.sig.dir <- 
                  paste0(save.patients.dir,
                         ifelse(current.sig,'significant/', 'not significant/'))
                dir.create(current.save.patients.sig.dir, showWarnings = FALSE, recursive = TRUE)
                
                # Begin figure save
                pdf(paste0(current.save.patients.sig.dir,
                           comparison.loop,' - ',
                           cluster.loop,' - ',
                           patient,' - ',
                           stage.loop,'_vE.pdf'),
                    width = 7.5 + (5 * length(lx.comparisons[[comparison.loop]])), 
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
                  x.limits <- c(-1000, 500)
                  # }else{
                  #   x.limits <- c(ifelse(lx.loop %in% leftmost.grid.plots, -1000, -500), 500)
                  # }
                  
                  # Save individual plot
                  if(! grid.loop){
                    pdf(paste0(save.patients.dir, 
                               comparison.loop,' - ',
                               cluster.loop,' - ',
                               patient,' - ',
                               stage.loop,' - ',
                               lx.loop,'_vE.pdf'),
                        width = 7, height = 4.5)
                  } # if(! grid.loop){
              
              plot.time.series(.y.values = ys[[lx.loop]],
                               .x.values = xs[[lx.loop]],
                               .sampling.rate = 256,
                               .colors = noun.colors$active[c('noun1','noun2'),'hex'],
                               .y.lwd = c(rep(3, times = length(patient.sl.stats)), 6),
                               .y.limits = y.limits,
                               .y.ticks = c(y.limits[1], y.limits[2]),
                               .y.tick.labels = c(ifelse(lx.loop %in% leftmost.grid.plots, y.limits[1], ''),
                                                  ifelse(lx.loop %in% leftmost.grid.plots, y.limits[2], '')),
                               .y.label = '',
                               .x.label = '',
                               .x.limits = x.limits,
                               .x.ticks = seq(x.limits[1], x.limits[2], by = 500),
                               .shuffle.dist.mean = shuffle.mean[[lx.loop]],
                               .shuffle.dist.error.bars.lower = shuffle.error.lower[[lx.loop]],
                               .shuffle.dist.error.bars.upper = shuffle.error.upper[[lx.loop]],
                               .shuffle.dist.color = 'grey',
                               .polygons.x = stage.times.range[[cluster.loop]][[stage.loop]],
                               .polygons.color = ifelse(theme.loop == 'white', 'black', 'white'),
                               show.polygons = ! grepl('verb',lx.loop),
                               .sig.windows = sig.windows[[lx.loop]],
                               .sig.color = colors[[theme.loop]]$rainbow_bright['pink','hex'],
                               show.sig.bars = FALSE,
                               show.sig.highlights = TRUE,
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
              
            }; rm(grid.loop)
            
          }#; rm(comparison.loop)
          
        }#; rm(patient.loop)
        
      }#; rm(theme.loop)
      
      
      # Progress update
      message('Completed: ',
              model.loop,' - ',
              band.loop,' - ',
              cluster.loop,' - ',
              stage.loop,'.',
              Sys.time())
      
    }; rm(stage.loop)
  }; rm(cluster.loop)
  
  
}; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


