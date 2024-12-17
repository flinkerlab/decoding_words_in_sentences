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
library('binom') # for binom.confint()


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
load(paste0(output.path,'data/0a - definitions/roi metadata/ROIs and temporal lobe splits.RData'))


### Loop thru models
model.types <- list.files(paste0(output.path,'data/9a - train classifiers - get variable importances - 50ms/weighted by aov.bfs/'))
# for(model.loop in model.types){ ## UNCOMMENT
model.loop = model.types[1] ## UNCOMMENT


### Get all significant windows and z-scored accuracies
combined.sig <- list()
patient.sig <- list()
message(band.loop,' - ',model.loop,' - Attaching significant windows... ',Sys.time())
for(roi.loop in rois){ ## UNCOMMENT
  # roi.loop = rois[1] ## UNCOMMENT
  
  ### Loop thru stages
  # Storage
  combined.sig[[roi.loop]] <- list()
  patient.sig[[roi.loop]] <- list()
  
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
    read.filename <- list.files(read.data.dir)[grepl("patient and combined sentence stats", list.files(read.data.dir))]
    if(read.filename != "patient and combined sentence stats - 1000 shuffles.RData"){message("Warning: ",model.loop,' - ',roi.loop,' - ',stage.loop,": Not 1000 shuffles!!!!")}
    attach(paste0(read.data.dir, read.filename)) 
    combined.sig[[roi.loop]][[stage.loop]] <- combined.sl.sig.windows
    patient.sig[[roi.loop]][[stage.loop]] <- patient.sl.sig.windows
    detach()
    
    # Clean up
    rm(read.filename)
    
  }; rm(stage.loop)
}; rm(roi.loop)
message('...detached! ', Sys.time()) # < 1min


### Plot grids
task.grids <- list(
  'active' = c('subject.noun.active','verb.active','object.noun.active'),
  'passive' = c('subject.noun.passive','verb.passive','object.noun.passive'),
  'sentence' = c('subject.noun','verb','object.noun'),
  'list' = c('list1','list2'))
task.grid.length <- lapply(task.grids, length)


###
### Define times and get significances
###

##
## For all individual classifiers (patient-specific)
##

# Storage
patient.sig.counts <- list()
patient.sig.proportions <- list()
patient.sig.proportions.95ci.upper <- list()
patient.sig.proportions.95ci.lower <- list()
patient.n.classifiers <- list()

# Fill in
for(task.loop in c('active','passive','sentence')){
  # task.loop = c('active','passive','sentence')[1]
  
  # Get current subject, verb, and object
  current.roles <- task.grids[[task.loop]]
  current.subject <- current.roles[grepl("subject", current.roles)]
  current.object <- current.roles[grepl("object", current.roles)]
  current.verb <- current.roles[grepl("verb", current.roles)]
  
  # Storage
  patient.sig.counts[[task.loop]] <- list()
  patient.sig.proportions[[task.loop]] <- list()
  patient.sig.proportions.95ci.upper[[task.loop]] <- list()
  patient.sig.proportions.95ci.lower[[task.loop]] <- list()
  patient.n.classifiers[[task.loop]] <- list()
  
  for(roi.loop in rois){
    # roi.loop = rois[3]
    
    # Storage
    patient.sig.counts[[task.loop]][[roi.loop]] <- list()
    patient.sig.proportions[[task.loop]][[roi.loop]] <- list()
    patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]] <- list()
    patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]] <- list()
    patient.n.classifiers[[task.loop]][[roi.loop]] <- list()
    
    for(role.loop in current.roles[! grepl("verb", current.roles)]){
      # role.loop = current.roles[! grepl("verb", current.roles)][2]
      
      # Storage
      patient.sig.counts[[task.loop]][[roi.loop]][[role.loop]] <- c()
      patient.sig.proportions[[task.loop]][[roi.loop]][[role.loop]] <- c()
      patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]][[role.loop]] <- c()
      patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]][[role.loop]] <- c()
      patient.n.classifiers[[task.loop]][[roi.loop]][[role.loop]] <- c()
      
      for(noun.loop in names(patient.sig[[1]][[1]][[1]][[1]])){
        # noun.loop = names(patient.sig[[1]][[1]][[1]][[1]])[2]
        
        ### For subjects:
        if(grepl("subject", role.loop)){ # If it's a subject
          
          # Get whether any significant detection per classifier
          current.sig <- list()
          for(train.loop in names(patient.sig[[roi.loop]])){
            # train.loop = names(patient.sig[[roi.loop]])[1]
            
            current.data <- patient.sig[[roi.loop]][[train.loop]]
            
            current.sig[[train.loop]] <- 
              lapply(current.data, function(.patient){
                
                # Significant if anything significant before 100ms post-onset or onset train time...
                .start.by <- max(100, stage.ranges[train.loop,'start.time'])
                .current.n.rows <- nrow(subset(.patient[[current.subject]][[noun.loop]], start.time < .start.by))
                # ...or anything significant before verb onset
                .current.n.rows <- .current.n.rows + 
                  nrow(subset(.patient[[current.verb]][[noun.loop]], end.time < 0))
                
                return(as.numeric(.current.n.rows > 0))
              })
            
            rm(current.data)
            
          }#; rm(train.loop)
        } # if(grepl("subject", role.loop)){
        
        
        ### For objects:
        if(grepl("object", role.loop)){ # If it's an object
          
          # Get whether any significant detection per classifier
          current.sig <- list()
          for(train.loop in names(patient.sig[[roi.loop]])){
            # train.loop = names(patient.sig[[roi.loop]])[1]
            
            current.data <- patient.sig[[roi.loop]][[train.loop]]
            
            # Significant if anything significant after 100ms pre-onset...
            current.sig[[train.loop]] <- 
              lapply(current.data, function(.patient){
                
                .current.n.rows <- nrow(subset(.patient[[current.object]][[noun.loop]], end.time > -100))
                # ...or anything significant after verb onset
                .current.n.rows <- .current.n.rows + 
                  nrow(subset(.patient[[current.verb]][[noun.loop]], start.time > 0))
                
                # Store
                return(as.numeric(.current.n.rows > 0))
              })
            
            rm(current.data)
            
          }#; rm(train.loop)
        } # if(grepl("object", role.loop)){
        
        # Stats
        n.sig.classifiers <- sum(unlist(current.sig))
        n.classifiers <- length(unlist(current.sig))
        
        # Store
        patient.sig.counts[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          n.sig.classifiers
        patient.sig.proportions[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          n.sig.classifiers / n.classifiers
        patient.n.classifiers[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          n.classifiers
        
        # Calculate 95% confidence interval
        current.ci <- binom.confint(x = n.sig.classifiers,
                                    n = n.classifiers,
                                    conf.level = .95, 
                                    methods = "wilson")
        patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          current.ci$upper - (n.sig.classifiers / n.classifiers)
        patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          (n.sig.classifiers / n.classifiers) - current.ci$lower
        
        # Clean up
        rm(n.sig.classifiers, n.classifiers, current.ci, current.sig)
        
      }#; rm(noun.loop)
    }#; rm(role.loop)
    
    # Reorganize
    patient.sig.counts[[task.loop]][[roi.loop]] <- 
      do.call(cbind, patient.sig.counts[[task.loop]][[roi.loop]])
    patient.sig.proportions[[task.loop]][[roi.loop]] <- 
      do.call(cbind, patient.sig.proportions[[task.loop]][[roi.loop]])
    patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]] <- 
      do.call(cbind, patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]])
    patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]] <- 
      do.call(cbind, patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]])
    patient.n.classifiers[[task.loop]][[roi.loop]] <- 
      do.call(cbind, patient.n.classifiers[[task.loop]][[roi.loop]])
    
  }#; rm(roi.loop)
}#; rm(task.loop)


### For lists:
task.loop <- 'list'
current.roles <- task.grids[[task.loop]]

# Storage
patient.sig.counts[[task.loop]] <- list()
patient.sig.proportions[[task.loop]] <- list()
patient.sig.proportions.95ci.upper[[task.loop]] <- list()
patient.sig.proportions.95ci.lower[[task.loop]] <- list()
patient.n.classifiers[[task.loop]] <- list()


for(roi.loop in rois){
  # roi.loop = rois[1]
  
  # Storage
  patient.sig.counts[[task.loop]][[roi.loop]] <- list()
  patient.sig.proportions[[task.loop]][[roi.loop]] <- list()
  patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]] <- list()
  patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]] <- list()
  patient.n.classifiers[[task.loop]][[roi.loop]] <- list()
  
  for(role.loop in current.roles){
    # role.loop = current.roles[1]
    
    # Storage
    patient.sig.counts[[task.loop]][[roi.loop]][[role.loop]] <- c()
    patient.sig.proportions[[task.loop]][[roi.loop]][[role.loop]] <- c()
    patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]][[role.loop]] <- c()
    patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]][[role.loop]] <- c()
    patient.n.classifiers[[task.loop]][[roi.loop]][[role.loop]] <- c()
    
    for(noun.loop in names(patient.sig[[1]][[1]][[1]][[1]])){
      # noun.loop = names(patient.sig[[1]][[1]][[1]][[1]])[1]
      
      # Get whether any significant detection per classifier
      current.sig <- list()
      for(train.loop in names(patient.sig[[roi.loop]])){
        # train.loop = names(patient.sig[[roi.loop]])[1]
        
        current.sig[[train.loop]] <- 
          lapply(patient.sig[[roi.loop]][[train.loop]], function(.patient){
            
            ## For list1
            if(role.loop == 'list1'){
              # Significant if anything significant before 0 or the onset of the training window if that's later
              .start.by <- max(0, stage.ranges[train.loop,'start.time'])
              .current.sig <- nrow(subset(.patient[['list1']][[noun.loop]], start.time < .start.by))
            }
            
            ## For list2
            if(role.loop == 'list2'){
              # Significant if sig after onset of list2... 
              .current.sig <- nrow(subset(.patient[['list2']][[noun.loop]], end.time > 0))
              # ...or after 350ms of list1
              .current.sig <- .current.sig +
                nrow(subset(.patient[['list1']][[noun.loop]], end.time > 350))
            }
            
            # Sig?
            .current.sig <- as.numeric(.current.sig > 0)
            
            # Output
            return(.current.sig)
          })
        
      }#; rm(train.loop)
      
      # Stats
      n.sig.classifiers <- sum(unlist(current.sig))
      n.classifiers <- length(unlist(current.sig))
      
      # Store
      patient.sig.counts[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
        n.sig.classifiers
      patient.sig.proportions[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
        n.sig.classifiers / n.classifiers
      patient.n.classifiers[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <-
        n.classifiers
      
      # Calculate 95% confidence interval
      current.ci <- binom.confint(x = n.sig.classifiers,
                                  n = n.classifiers,
                                  conf.level = .95, 
                                  methods = "wilson")
      patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
        current.ci$upper - (n.sig.classifiers / n.classifiers)
      patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
        (n.sig.classifiers / n.classifiers) - current.ci$lower
      
      # Clean up
      rm(n.sig.classifiers, n.classifiers, current.ci, current.sig)
      
    }#; rm(noun.loop)
  }#; rm(role.loop)
  
  # Reorganize
  patient.sig.counts[[task.loop]][[roi.loop]] <- 
    do.call(cbind, patient.sig.counts[[task.loop]][[roi.loop]])
  patient.sig.proportions[[task.loop]][[roi.loop]] <- 
    do.call(cbind, patient.sig.proportions[[task.loop]][[roi.loop]])
  patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]] <- 
    do.call(cbind, patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]])
  patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]] <- 
    do.call(cbind, patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]])
  patient.n.classifiers[[task.loop]][[roi.loop]] <- 
    do.call(cbind, patient.n.classifiers[[task.loop]][[roi.loop]])
  
}#; rm(roi.loop)

# Clean up
rm(task.loop, current.roles)



##
## Collapsing across patients
##

# Storage
combined.sig.counts <- list()
combined.sig.proportions <- list()
combined.sig.proportions.95ci.upper <- list()
combined.sig.proportions.95ci.lower <- list()
combined.n.classifiers <- list()

# Fill in
for(task.loop in c('active','passive','sentence')){
  # task.loop = c('active','passive','sentence')[2]
  
  # Get current subject, verb, and object
  current.roles <- task.grids[[task.loop]]
  current.subject <- current.roles[grepl("subject", current.roles)]
  current.object <- current.roles[grepl("object", current.roles)]
  current.verb <- current.roles[grepl("verb", current.roles)]
  
  # Storage
  combined.sig.counts[[task.loop]] <- list()
  combined.sig.proportions[[task.loop]] <- list()
  combined.sig.proportions.95ci.upper[[task.loop]] <- list()
  combined.sig.proportions.95ci.lower[[task.loop]] <- list()
  combined.n.classifiers[[task.loop]] <- list()
  
  for(roi.loop in rois){
    # roi.loop = rois[1]
    
    # Storage
    combined.sig.counts[[task.loop]][[roi.loop]] <- list()
    combined.sig.proportions[[task.loop]][[roi.loop]] <- list()
    combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]] <- list()
    combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]] <- list()
    combined.n.classifiers[[task.loop]][[roi.loop]] <- list()
    
    for(role.loop in current.roles[! grepl("verb", current.roles)]){
      # role.loop = current.roles[! grepl("verb", current.roles)][1]
      
      # Storage
      combined.sig.counts[[task.loop]][[roi.loop]][[role.loop]] <- c()
      combined.sig.proportions[[task.loop]][[roi.loop]][[role.loop]] <- c()
      combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]][[role.loop]] <- c()
      combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]][[role.loop]] <- c()
      combined.n.classifiers[[task.loop]][[roi.loop]][[role.loop]] <- c()
      
      for(noun.loop in names(combined.sig[[1]][[1]][[1]])){
        # noun.loop = names(combined.sig[[1]][[1]][[1]])[1]
        
        ### For subjects:
        if(grepl("subject", role.loop)){ # If it's a subject
          
          # Get whether any significant detection per classifier
          current.sig <- list()
          for(train.loop in names(combined.sig[[roi.loop]])){
            # train.loop = names(combined.sig[[roi.loop]])[8]
            
            current.data <- combined.sig[[roi.loop]][[train.loop]]
            
            # Significant if anything significant before 100ms post-onset or onset train time...
            .start.by <- max(100, stage.ranges[train.loop,'start.time'])
            current.sig[[train.loop]] <- 
              nrow(subset(current.data[[current.subject]][[noun.loop]], start.time < .start.by))
            # ...or anything significant before verb onset
            current.sig[[train.loop]] <- current.sig[[train.loop]] + 
              nrow(subset(current.data[[current.verb]][[noun.loop]], end.time < 0))
            
            # Store
            current.sig[[train.loop]] <- as.numeric(current.sig[[train.loop]] > 0)
            
            rm(current.data, .start.by)
            
          }#; rm(train.loop)
        } # if(grepl("subject", role.loop)){
        
        
        ### For objects:
        if(grepl("object", role.loop)){ # If it's an object
          
          # Get whether any significant detection per classifier
          current.sig <- list()
          for(train.loop in names(combined.sig[[roi.loop]])){
            # train.loop = names(combined.sig[[roi.loop]])[1]
            
            current.data <- combined.sig[[roi.loop]][[train.loop]]
            
            # Significant if anything significant after 100ms pre-onset...
            current.sig[[train.loop]] <- 
              nrow(subset(current.data[[current.object]][[noun.loop]], end.time > -100))
            # ...or anything significant after verb onset
            current.sig[[train.loop]] <- current.sig[[train.loop]] + 
              nrow(subset(current.data[[current.verb]][[noun.loop]], start.time > 0))
            
            # Store
            current.sig[[train.loop]] <- as.numeric(current.sig[[train.loop]] > 0)
            
            rm(current.data)
            
          }#; rm(train.loop)
        } # if(grepl("object", role.loop)){
        
        # Stats
        n.sig.classifiers <- sum(unlist(current.sig))
        n.classifiers <- length(unlist(current.sig))
        
        # Store
        combined.sig.counts[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          n.sig.classifiers
        combined.sig.proportions[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          n.sig.classifiers / n.classifiers
        combined.n.classifiers[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <-
          n.classifiers
        
        # Calculate 95% confidence interval
        current.ci <- binom.confint(x = n.sig.classifiers,
                                    n = n.classifiers,
                                    conf.level = .95, 
                                    methods = "wilson")
        combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          current.ci$upper - (n.sig.classifiers / n.classifiers)
        combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
          (n.sig.classifiers / n.classifiers) - current.ci$lower
        
        # Clean up
        rm(n.sig.classifiers, n.classifiers, current.ci, current.sig)
        
      }#; rm(noun.loop)
    }#; rm(role.loop)
    
    # Reorganize
    combined.sig.counts[[task.loop]][[roi.loop]] <- 
      do.call(cbind, combined.sig.counts[[task.loop]][[roi.loop]])
    combined.sig.proportions[[task.loop]][[roi.loop]] <- 
      do.call(cbind, combined.sig.proportions[[task.loop]][[roi.loop]])
    combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]] <- 
      do.call(cbind, combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]])
    combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]] <- 
      do.call(cbind, combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]])
    combined.n.classifiers[[task.loop]][[roi.loop]] <- 
      do.call(cbind, combined.n.classifiers[[task.loop]][[roi.loop]])
    
  }#; rm(roi.loop)
}#; rm(task.loop)


### For lists:
task.loop <- 'list'
current.roles <- task.grids[[task.loop]]

# Storage
combined.sig.counts[[task.loop]] <- list()
combined.sig.proportions[[task.loop]] <- list()
combined.sig.proportions.95ci.upper[[task.loop]] <- list()
combined.sig.proportions.95ci.lower[[task.loop]] <- list()
combined.n.classifiers[[task.loop]] <- list()

for(roi.loop in rois){
  # roi.loop = rois[3]
  
  # Storage
  combined.sig.counts[[task.loop]][[roi.loop]] <- list()
  combined.sig.proportions[[task.loop]][[roi.loop]] <- list()
  combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]] <- list()
  combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]] <- list()
  combined.n.classifiers[[task.loop]][[roi.loop]] <- list()
  
  for(role.loop in current.roles){
    # role.loop = current.roles[1]
    
    # Storage
    combined.sig.counts[[task.loop]][[roi.loop]][[role.loop]] <- c()
    combined.sig.proportions[[task.loop]][[roi.loop]][[role.loop]] <- c()
    combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]][[role.loop]] <- c()
    combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]][[role.loop]] <- c()
    combined.n.classifiers[[task.loop]][[roi.loop]][[role.loop]] <- list()
    
    for(noun.loop in names(combined.sig[[1]][[1]][[1]])){
      # noun.loop = names(combined.sig[[1]][[1]][[1]])[1]
      
      # Get whether any significant detection per classifier
      current.sig <- list()
      for(train.loop in names(combined.sig[[roi.loop]])){
        # train.loop = names(combined.sig[[roi.loop]])[1]
        
        current.data <- combined.sig[[roi.loop]][[train.loop]]
        
        ## For list1
        if(role.loop == 'list1'){
          # Significant if anything significant before 0 or the onset of the training window if that's later
          .start.by <- max(0, stage.ranges[train.loop,'start.time'])
          current.sig[[train.loop]] <- nrow(subset(current.data[['list1']][[noun.loop]], start.time < .start.by))
        }
        
        ## For list2
        if(role.loop == 'list2'){
          # Significant if sig after onset of list2... 
          current.sig[[train.loop]] <- nrow(subset(current.data[['list2']][[noun.loop]], end.time > 0))
          # ...or after 350ms of list1
          current.sig[[train.loop]] <- current.sig[[train.loop]] +
            nrow(subset(current.data[['list1']][[noun.loop]], end.time > 350))
        }
        
        # Sig?
        current.sig[[train.loop]] <- as.numeric(current.sig[[train.loop]] > 0)
        
        # Clean up
        rm(current.data)
        
      }#; rm(train.loop)
      
      # Stats
      n.sig.classifiers <- sum(unlist(current.sig))
      n.classifiers <- length(unlist(current.sig))
      
      # Store
      combined.sig.counts[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
        n.sig.classifiers
      combined.sig.proportions[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
        n.sig.classifiers / n.classifiers
      combined.n.classifiers[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <-
        n.classifiers
      
      # Calculate 95% confidence interval
      current.ci <- binom.confint(x = n.sig.classifiers,
                                  n = n.classifiers,
                                  conf.level = .95, 
                                  methods = "wilson")
      combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
        current.ci$upper - (n.sig.classifiers / n.classifiers)
      combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]][[role.loop]][noun.loop] <- 
        (n.sig.classifiers / n.classifiers) - current.ci$lower
      
      # Clean up
      rm(n.sig.classifiers, n.classifiers, current.ci, current.sig)
      
    }#; rm(noun.loop)
  }#; rm(role.loop)
  
  # Reorganize
  combined.sig.counts[[task.loop]][[roi.loop]] <- 
    do.call(cbind, combined.sig.counts[[task.loop]][[roi.loop]])
  combined.sig.proportions[[task.loop]][[roi.loop]] <- 
    do.call(cbind, combined.sig.proportions[[task.loop]][[roi.loop]])
  combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]] <- 
    do.call(cbind, combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]])
  combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]] <- 
    do.call(cbind, combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]])
  combined.n.classifiers[[task.loop]][[roi.loop]] <- 
    do.call(cbind, combined.n.classifiers[[task.loop]][[roi.loop]])
  
}#; rm(roi.loop)

# Clean up
rm(task.loop, current.roles)


###
### Plots
###

text.size.small <- 1.4
text.size.med <- 1.6
zoom <- 1

##
## Patient stats
##

for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  # Save dir
  bar.plot.dir.rois <- paste0(output.path, 
                              'figures/15a - barplots - proportion classifiers sig/',
                              theme.loop,
                              '/all patients - by roi/')
  dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)
  
  for(roi.loop in rois){
    # roi.loop = rois[3]
    
    pdf(paste0(bar.plot.dir.rois, 'proportion of classifiers sig - ',roi.loop,' - ',model.loop,'.pdf'),
        width = 6 * 4, height = 6)
    par(oma = c(0,0,2,0),
        mar = c(3,3,1,1) * zoom,
        mfrow = c(1, 4),
        lwd = 2 * zoom,
        bg = rgb(1,1,1,0))
    
    for(task.loop in c('active','passive','list','sentence')){
      # task.loop = names(task.grids)[1]
      
      ys <- patient.sig.proportions[[task.loop]][[roi.loop]]
      current.errors.upper <- patient.sig.proportions.95ci.upper[[task.loop]][[roi.loop]]
      current.errors.lower <- patient.sig.proportions.95ci.lower[[task.loop]][[roi.loop]]
      
      xs <- barplot(ys,
                    ylim = c(0, .3),
                    beside = TRUE,
                    col = noun.colors[[task.loop]][rownames(ys),'hex'],
                    border = theme.loop,
                    cex.names = text.size.small * zoom,
                    yaxt = 'n',
                    las = 1)
      # Error bars
      for(row.loop in 1:nrow(xs)){
        for(col.loop in 1:ncol(xs)){
          arrows(x0 = xs[row.loop, col.loop],
                 y0 = ys[row.loop, col.loop] - current.errors.lower[row.loop, col.loop],
                 y1 = ys[row.loop, col.loop] + current.errors.upper[row.loop, col.loop],
                 angle = 90,
                 code = 3,
                 length = 0, 
                 col = rgb(.7,.7,.7))
        }; rm(col.loop)
      }; rm(row.loop)
      
      axis(side = 2,
           at = c(0, .3),
           # labels = .y.tick.labels,
           las = 0,
           tck = -.025 * zoom, # length of tick
           padj = -.45 * zoom, # distance between tick and label
           lwd = 1.5 * zoom,
           lwd.ticks = 1.5 * zoom,
           cex.axis = text.size.small * zoom,
           col = ifelse(theme.loop == 'white', 'black', 'white'),
           col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
      add.text.line.multiple.colors(text.segments = rownames(ys),
                                    text.colors = noun.colors[[task.loop]][rownames(ys),'hex'],
                                    .side = 3,
                                    .line = 1,
                                    .outer = TRUE)
      
    }#; rm(task.loop)
    
    dev.off()  
    
  }#; rm(roi.loop)
  
}#; rm(theme.loop)



##
## Patient stats - just active and passive
##

for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  # Save dir
  bar.plot.dir.rois <- paste0(output.path, 
                              'figures/15a - barplots - proportion classifiers sig/',
                              theme.loop,
                              '/all patients - by roi - just actives and passives/')
  dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)
  
  for(roi.loop in rois){
    # roi.loop = rois[3]
    
    pdf(paste0(bar.plot.dir.rois, 'proportion of classifiers sig - ',roi.loop,' - ',model.loop,'.pdf'),
        width = 6 * 2, height = 8)
    par(oma = c(0,0,2,0),
        mar = c(3,3,1,1) * zoom,
        mfrow = c(1, 2),
        lwd = 2 * zoom,
        bg = rgb(1,1,1,0))
    
    current.tasks <- c('active','passive')
    ys <- lapply(patient.sig.proportions[current.tasks], function(x){x[[roi.loop]]})
    current.errors.upper <- lapply(patient.sig.proportions.95ci.upper[current.tasks], function(x){x[[roi.loop]]})
    current.errors.lower <- lapply(patient.sig.proportions.95ci.lower[current.tasks], function(x){x[[roi.loop]]})
    y.max <- max(unlist(ys) + unlist(current.errors.upper))
    
    for(task.loop in current.tasks){
      # task.loop = current.tasks[1]
      
      xs <- barplot(ys[[task.loop]],
                    ylim = c(0, ifelse(roi.loop == 'SMC', .25, .15)),
                    beside = TRUE,
                    col = noun.colors[[task.loop]][rownames(ys[[task.loop]]),'hex'],
                    border = theme.loop,
                    cex.names = text.size.small * zoom,
                    yaxt = 'n',
                    las = 1)
      # Error bars
      for(row.loop in 1:nrow(xs)){
        for(col.loop in 1:ncol(xs)){
          arrows(x0 = xs[row.loop, col.loop],
                 y0 = ys[[task.loop]][row.loop, col.loop] - current.errors.lower[[task.loop]][row.loop, col.loop],
                 y1 = ys[[task.loop]][row.loop, col.loop] + current.errors.upper[[task.loop]][row.loop, col.loop],
                 angle = 90,
                 code = 3,
                 length = 0, 
                 col = rgb(.7,.7,.7))
        }; rm(col.loop)
      }; rm(row.loop)
      
      axis(side = 2,
           # at = c(0, max(.1, floor(100 * y.max / 5) * 5 / 100)),
           at = c(0, ifelse(roi.loop == 'SMC', .25, .15)),
           # labels = .y.tick.labels,
           las = 0,
           tck = -.025 * zoom, # length of tick
           padj = -.45 * zoom, # distance between tick and label
           lwd = 1.5 * zoom,
           lwd.ticks = 1.5 * zoom,
           cex.axis = text.size.small * zoom,
           col = ifelse(theme.loop == 'white', 'black', 'white'),
           col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
      add.text.line.multiple.colors(text.segments = rownames(ys[[task.loop]]),
                                    text.colors = noun.colors[[task.loop]][rownames(ys[[task.loop]]),'hex'],
                                    .side = 3,
                                    .line = 1,
                                    .outer = TRUE)
      
    }#; rm(task.loop)
    
    dev.off()  
    
  }#; rm(roi.loop)
  
}#; rm(theme.loop)


##
## Patient stats, collapsing across ROI
##

for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  # Save dir
  bar.plot.dir.rois <- paste0(output.path, 
                              'figures/15a - barplots - proportion classifiers sig/',
                              theme.loop,
                              '/all patients - collapsed rois/')
  dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)
  
  pdf(paste0(bar.plot.dir.rois, 'proportion of classifiers sig - ',model.loop,'.pdf'),
      width = 6 * 4, height = 6)
  par(oma = c(0,0,2,0),
      mar = c(3,3,1,1) * zoom,
      mfrow = c(1, 4),
      lwd = 2 * zoom,
      bg = rgb(1,1,1,0))
  
  for(task.loop in c('active','passive','list','sentence')){
    # task.loop = names(task.grids)[2]
    
    ys.sum <- elementwise.matrix.apply(patient.sig.counts[[task.loop]], .function = "sum")
    ys.denom <- elementwise.matrix.apply(patient.n.classifiers[[task.loop]], .function = "sum")
    ys <- ys.sum / ys.denom
    current.errors <- binom.confint(x = as.vector(ys.sum),
                                    n = as.vector(ys.denom),
                                    conf.level = .95, 
                                    methods = "wilson")
    current.errors.upper <- matrix(data = current.errors$upper, nrow = 2) - ys
    current.errors.lower <- ys - matrix(data = current.errors$lower, nrow = 2)
    
    xs <- barplot(ys,
                  ylim = c(0, .3),
                  beside = TRUE,
                  col = noun.colors[[task.loop]][rownames(ys),'hex'],
                  border = theme.loop,
                  cex.names = text.size.small * zoom,
                  yaxt = 'n',
                  las = 1)
    # Error bars
    for(row.loop in 1:nrow(xs)){
      for(col.loop in 1:ncol(xs)){
        arrows(x0 = xs[row.loop, col.loop],
               y0 = ys[row.loop, col.loop] - current.errors.lower[row.loop, col.loop],
               y1 = ys[row.loop, col.loop] + current.errors.upper[row.loop, col.loop],
               angle = 90,
               code = 3,
               length = 0, 
               col = rgb(.7,.7,.7))
      }; rm(col.loop)
    }; rm(row.loop)
    
    axis(side = 2,
         at = c(0, .3),
         # labels = .y.tick.labels,
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.small * zoom,
         col = ifelse(theme.loop == 'white', 'black', 'white'),
         col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
    add.text.line.multiple.colors(text.segments = rownames(ys),
                                  text.colors = noun.colors[[task.loop]][rownames(ys),'hex'],
                                  .side = 3,
                                  .line = 1,
                                  .outer = TRUE)
    
  }#; rm(task.loop)
  
  dev.off()  
  
}#; rm(theme.loop)


##
## Patient stats, collapsing across ROI - COUNT
##

for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  ## Data and y-limits
  ys <- ys.denom <- current.errors.upper <- current.errors.lower <- list()
  for(task.loop in c('active','passive','list','sentence')){
    # task.loop = names(task.grids)[1]
    
    ys[[task.loop]] <- elementwise.matrix.apply(patient.sig.counts[[task.loop]], .function = "sum")
    ys.denom[[task.loop]] <- elementwise.matrix.apply(patient.n.classifiers[[task.loop]], .function = "sum")
    current.errors <- binom.confint(x = as.vector(ys[[task.loop]]),
                                    n = as.vector(ys.denom[[task.loop]]),
                                    conf.level = .95, 
                                    methods = "wilson")
    current.errors.upper[[task.loop]] <- 
      (matrix(data = current.errors$upper * ys.denom[[task.loop]], nrow = 2) - ys[[task.loop]])
    current.errors.lower[[task.loop]] <- 
      (ys[[task.loop]] - matrix(data = current.errors$lower * ys.denom[[task.loop]], nrow = 2))
  }#; rm(task.loop)
  y.limits <- c(0, ceiling(max(unlist(ys) + unlist(current.errors.upper)) / 5) * 5)
  
  # Save dir
  bar.plot.dir.rois <- paste0(output.path, 
                              'figures/15a - barplots - proportion classifiers sig/',
                              theme.loop,
                              '/all patients - collapsed rois - counts/')
  dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)
  
  pdf(paste0(bar.plot.dir.rois, 'proportion of classifiers sig - ',model.loop,' - n classifiers total=',unique(unlist(ys.denom)),'.pdf'),
      width = 6 * 4, height = 6)
  par(oma = c(0,0,2,0),
      mar = c(3,3,1,1) * zoom,
      mfrow = c(1, 4),
      lwd = 2 * zoom,
      bg = rgb(1,1,1,0))
  
  for(task.loop in c('active','passive','list','sentence')){
    # task.loop = names(task.grids)[2]
    
    xs <- barplot(ys[[task.loop]],
                  ylim = y.limits,
                  beside = TRUE,
                  col = noun.colors[[task.loop]][rownames(ys[[task.loop]]),'hex'],
                  border = theme.loop,
                  cex.names = text.size.small * zoom,
                  yaxt = 'n',
                  las = 1)
    # Error bars
    for(row.loop in 1:nrow(xs)){
      for(col.loop in 1:ncol(xs)){
        arrows(x0 = xs[row.loop, col.loop],
               y0 = ys[[task.loop]][row.loop, col.loop] - current.errors.lower[[task.loop]][row.loop, col.loop],
               y1 = ys[[task.loop]][row.loop, col.loop] + current.errors.upper[[task.loop]][row.loop, col.loop],
               angle = 90,
               code = 3,
               length = 0, 
               col = rgb(.7,.7,.7))
      }; rm(col.loop)
    }; rm(row.loop)
    
    axis(side = 2,
         at = y.limits,
         # labels = .y.tick.labels,
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.small * zoom,
         col = ifelse(theme.loop == 'white', 'black', 'white'),
         col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
    add.text.line.multiple.colors(text.segments = rownames(ys),
                                  text.colors = noun.colors[[task.loop]][rownames(ys),'hex'],
                                  .side = 3,
                                  .line = 1,
                                  .outer = TRUE)
    
  }#; rm(task.loop)
  
  dev.off()  
  
}#; rm(theme.loop)



##
## Patient stats, collapsing across ROI - grouped by congruent/incongruent - COUNT
##

zoom <- 1.8

## Data and y-limits
ys <- ys.denom <- current.errors.upper <- current.errors.lower <- list()
for(task.loop in c(#'list',
  'active','passive')){
  # task.loop = names(task.grids)[1]
  
  ys[[task.loop]] <- elementwise.matrix.apply(patient.sig.counts[[task.loop]], .function = "sum")
  ys.denom[[task.loop]] <- elementwise.matrix.apply(patient.n.classifiers[[task.loop]], .function = "sum")
  current.errors <- binom.confint(x = as.vector(ys[[task.loop]]),
                                  n = as.vector(ys.denom[[task.loop]]),
                                  conf.level = .95, 
                                  methods = "wilson")
  current.errors.upper[[task.loop]] <- 
    (matrix(data = current.errors$upper * ys.denom[[task.loop]], nrow = 2) - ys[[task.loop]])
  current.errors.lower[[task.loop]] <- 
    (ys[[task.loop]] - matrix(data = current.errors$lower * ys.denom[[task.loop]], nrow = 2))
}#; rm(task.loop)
y.limits <- c(0, ceiling(max(unlist(ys) + unlist(current.errors.upper))))
y.ticks <- c(0, floor(max(unlist(ys) + unlist(current.errors.upper)) / 10) * 10)


counts <- list()
uppers <- list()
lowers <- list()

counts[['congruent']] <- 
  cbind(ld2a(lapply(ys, as.data.frame))[1,1,],
        ld2a(lapply(ys, as.data.frame))[2,2,])
uppers[['congruent']] <- 
  cbind(ld2a(lapply(current.errors.upper, as.data.frame))[1,1,],
        ld2a(lapply(current.errors.upper, as.data.frame))[2,2,])
lowers[['congruent']] <- 
  cbind(ld2a(lapply(current.errors.lower, as.data.frame))[1,1,],
        ld2a(lapply(current.errors.lower, as.data.frame))[2,2,])

counts[['incongruent']] <- 
  cbind(ld2a(lapply(ys, as.data.frame))[1,2,],
        ld2a(lapply(ys, as.data.frame))[2,1,])
uppers[['incongruent']] <- 
  cbind(ld2a(lapply(current.errors.upper, as.data.frame))[1,2,],
        ld2a(lapply(current.errors.upper, as.data.frame))[2,1,])
lowers[['incongruent']] <- 
  cbind(ld2a(lapply(current.errors.lower, as.data.frame))[1,2,],
        ld2a(lapply(current.errors.lower, as.data.frame))[2,1,])

### Stats
## Comparisons:
##  - incongruent passive during subject to incongruent active during subject
prop.test(x = c(counts$incongruent['passive',1], 
                counts$incongruent['active',1]),
          n = c(ys.denom$passive['noun2','subject.noun.passive'],
                ys.denom$active['noun2','subject.noun.active']),
          alternative = 'greater')
##  - incongruent passive during object to incongruent active during object
prop.test(x = c(counts$incongruent['passive',2], 
                counts$incongruent['active',2]),
          n = c(ys.denom$passive['noun1','object.noun.passive'],
                ys.denom$active['noun1','object.noun.active']),
          alternative = 'greater')
##  - incongruent passive during subject to congruent passive during subject
prop.test(x = c(counts$incongruent['passive',2], 
                counts$congruent['passive',1]),
          n = c(ys.denom$passive['noun2','subject.noun.passive'],
                ys.denom$passive['noun1','subject.noun.passive']),
          alternative = 'greater')





for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  # Save dir
  bar.plot.dir.rois <- paste0(output.path, 
                              'figures/15a - barplots - proportion classifiers sig/',
                              theme.loop,
                              '/all patients - collapsed rois - by congruent noun - counts/')
  dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)
  
  for(cong.loop in names(counts)){
    # cong.loop = names(counts)[1]
    
    pdf(paste0(bar.plot.dir.rois, 'proportion of classifiers sig - ',model.loop,' - n classifiers total=',unique(unlist(ys.denom)),' - ',cong.loop,'.pdf'),
        width = 7, height = 4.5) # was 6 and 5.5
    par(oma = c(0,0,0,0),
        mar = c(1,3,1,3) * zoom,
        mfrow = c(1, 1),
        lwd = 2 * zoom,
        bg = rgb(1,1,1,0))
    
    xs <- barplot(counts[[cong.loop]],
                  ylim = y.limits,
                  col = c(rgb(.25,.25,.25), rgb(.75,.75,.75)),
                  beside = TRUE,
                  # col = noun.colors[[task.loop]][rownames(ys[[task.loop]]),'hex'],
                  border = theme.loop,
                  lwd = 2 * zoom,
                  cex.names = text.size.small * zoom,
                  yaxt = 'n',
                  las = 1)
    # Error bars
    for(row.loop in 1:nrow(xs)){
      for(col.loop in 1:ncol(xs)){
        arrows(x0 = xs[row.loop, col.loop],
               y0 = counts[[cong.loop]][row.loop, col.loop] - lowers[[cong.loop]][row.loop, col.loop],
               y1 = counts[[cong.loop]][row.loop, col.loop] + uppers[[cong.loop]][row.loop, col.loop],
               lwd = 2 * zoom,
               angle = 90,
               code = 3,
               length = 0, 
               col = rgb(.5,.5,.5))
      }; rm(col.loop)
    }; rm(row.loop)
    
    axis(side = 2,
         at = y.ticks,
         # labels = .y.tick.labels,
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.small * zoom,
         col = ifelse(theme.loop == 'white', 'black', 'white'),
         col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
    add.text.line.multiple.colors(text.segments = rownames(ys),
                                  text.colors = noun.colors[[task.loop]][rownames(ys),'hex'],
                                  .side = 3,
                                  .line = 1,
                                  .outer = TRUE)
    
    dev.off()  
    
  }#; rm(cong.loop)
  
}#; rm(theme.loop)




##
## Patient stats, stacking ROIs - grouped by congruent/incongruent - COUNT
##

zoom <- 1.8

## Data and y-limits
y.roi.counts <- ys <- ys.denom <- current.errors.upper <- current.errors.lower <- list()
for(task.loop in c(#'list',
  'active','passive')){
  # task.loop = names(task.grids)[1]
  
  y.roi.counts[[task.loop]] <- ld2a(lapply(patient.sig.counts[[task.loop]], data.frame))
  ys[[task.loop]] <- elementwise.matrix.apply(patient.sig.counts[[task.loop]], .function = "sum")
  ys.denom[[task.loop]] <- elementwise.matrix.apply(patient.n.classifiers[[task.loop]], .function = "sum")
  current.errors <- binom.confint(x = as.vector(ys[[task.loop]]),
                                  n = as.vector(ys.denom[[task.loop]]),
                                  conf.level = .95, 
                                  methods = "wilson")
  current.errors.upper[[task.loop]] <- 
    (matrix(data = current.errors$upper * ys.denom[[task.loop]], nrow = 2) - ys[[task.loop]])
  current.errors.lower[[task.loop]] <- 
    (ys[[task.loop]] - matrix(data = current.errors$lower * ys.denom[[task.loop]], nrow = 2))
}#; rm(task.loop)
y.limits <- c(0, ceiling(max(unlist(ys) + unlist(current.errors.upper))))
y.ticks <- c(0, floor(max(unlist(ys) + unlist(current.errors.upper)) / 10) * 10)

roi.counts <- list()
counts <- list()
uppers <- list()
lowers <- list()

roi.counts[['congruent']] <- 
  roi.counts[['incongruent']] <- list()
rois.in.order <- rev(dimnames(y.roi.counts$active)[[3]])
rois.in.order <- c(rois.in.order[-which(rois.in.order == 'MFG')], 'MFG')
for(roi.loop in rois.in.order){
  roi.counts[['congruent']][[roi.loop]] <- 
    cbind(ld2a(lapply(y.roi.counts, function(.voice){as.data.frame(.voice[,,roi.loop])}))[1,1,],
          ld2a(lapply(y.roi.counts, function(.voice){as.data.frame(.voice[,,roi.loop])}))[2,2,])
  roi.counts[['incongruent']][[roi.loop]] <- 
    cbind(ld2a(lapply(y.roi.counts, function(.voice){as.data.frame(.voice[,,roi.loop])}))[1,2,],
          ld2a(lapply(y.roi.counts, function(.voice){as.data.frame(.voice[,,roi.loop])}))[2,1,])
}; rm(roi.loop)

counts[['congruent']] <- 
  cbind(ld2a(lapply(ys, as.data.frame))[1,1,],
        ld2a(lapply(ys, as.data.frame))[2,2,])
uppers[['congruent']] <- 
  cbind(ld2a(lapply(current.errors.upper, as.data.frame))[1,1,],
        ld2a(lapply(current.errors.upper, as.data.frame))[2,2,])
lowers[['congruent']] <- 
  cbind(ld2a(lapply(current.errors.lower, as.data.frame))[1,1,],
        ld2a(lapply(current.errors.lower, as.data.frame))[2,2,])

counts[['incongruent']] <- 
  cbind(ld2a(lapply(ys, as.data.frame))[1,2,],
        ld2a(lapply(ys, as.data.frame))[2,1,])
uppers[['incongruent']] <- 
  cbind(ld2a(lapply(current.errors.upper, as.data.frame))[1,2,],
        ld2a(lapply(current.errors.upper, as.data.frame))[2,1,])
lowers[['incongruent']] <- 
  cbind(ld2a(lapply(current.errors.lower, as.data.frame))[1,2,],
        ld2a(lapply(current.errors.lower, as.data.frame))[2,1,])

prop.stats <- data.frame(matrix(nrow = length(rois), ncol = 2,
                                dimnames = list(rois,c('subject','object'))))
prop.stats <- list('congruent' = list('active' = prop.stats, 'passive' = prop.stats),
                   'incongruent' = list('active' = prop.stats, 'passive' = prop.stats))
for(cong.loop in names(roi.counts)){
  # cong.loop = names(counts)[2]
  for(syntax.loop in c('active','passive')){
    # syntax.loop = c('active','passive')[1]
    for(role.loop in c(1,2)){
      # role.loop = 1
      
      current.n.sig <- sapply(roi.counts[[cong.loop]], function(.roi){unname(.roi[syntax.loop, role.loop])})
      current.n.total <- sapply(patient.n.classifiers[[syntax.loop]], function(.roi){.roi[1,1]})
      
      for(roi.loop in rois){
        # roi.loop = rois[1]
        
        current.roi.n.sig <- current.n.sig[roi.loop]
        current.roi.n.total <- current.n.total[roi.loop]
        current.other.n.sig <- sum(current.n.sig) - current.roi.n.sig
        current.other.n.total <- sum(current.n.total) - current.roi.n.total
        
        prop.stats[[cong.loop]][[syntax.loop]][roi.loop,ifelse(role.loop == 1, 'subject','object')] <-
          prop.test(x = c(current.roi.n.sig, current.other.n.sig),
                    n = c(current.roi.n.total, current.other.n.total),
                    alternative = "greater")$p.value
        
      }; rm(roi.loop)
      
      # FDR Correct p-values across ROIs
      prop.stats[[cong.loop]][[syntax.loop]][,role.loop] <- 
        p.adjust(prop.stats[[cong.loop]][[syntax.loop]][,role.loop], method = 'fdr')
      
    }; rm(role.loop)
  }; rm(syntax.loop)
}; rm(cong.loop)


for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  # Save dir
  bar.plot.dir.rois.stacked <- paste0(output.path, 
                                      'figures/15a - barplots - proportion classifiers sig/',
                                      theme.loop,
                                      '/all patients - rois stacked - by congruent noun - counts/')
  dir.create(bar.plot.dir.rois.stacked, showWarnings = FALSE, recursive = TRUE)
  
  for(cong.loop in names(roi.counts)){
    # cong.loop = names(counts)[2]
    
    ##
    ## ALL ROIS
    ##
    
    pdf(paste0(bar.plot.dir.rois.stacked, 'count of classifiers sig - ',model.loop,' - n classifiers total=',unique(unlist(ys.denom)),' - ',cong.loop,'.pdf'),
        width = 6, height = 5.5)
    par(oma = c(0,0,0,0),
        mar = c(1,3,1,3) * zoom,
        mfrow = c(1, 1),
        lwd = 2 * zoom,
        bg = rgb(1,1,1,0))
    
    xs <- barplot(counts[[cong.loop]],
                  ylim = y.limits,
                  col = rgb(1,1,1,0),
                  beside = TRUE,
                  border = theme.loop,
                  lwd = 2 * zoom,
                  cex.names = text.size.small * zoom,
                  yaxt = 'n',
                  las = 1)
    
    top.vals <- counts[[cong.loop]]
    for(roi.loop in names(roi.counts[[cong.loop]])){
      # roi.loop = names(roi.counts[[cong.loop]])[1]
      
      barplot(top.vals,
              # col = c(rgb(.25,.25,.25), rgb(.75,.75,.75)),
              col = c(roi.colors$rostral.caudal[roi.loop],
                      colorspace::darken(roi.colors$rostral.caudal[roi.loop], 
                                         amount = .1, method = "relative")),
              beside = TRUE,
              border = theme.loop,
              lwd = 2 * zoom,
              # cex.names = text.size.small * zoom,
              yaxt = 'n',
              las = 1,
              add = TRUE)
      
      top.vals <- top.vals - roi.counts[[cong.loop]][[roi.loop]]
      
    }; rm(roi.loop)
    
    # Error bars
    for(row.loop in 1:nrow(xs)){
      for(col.loop in 1:ncol(xs)){
        arrows(x0 = xs[row.loop, col.loop],
               y0 = counts[[cong.loop]][row.loop, col.loop] - lowers[[cong.loop]][row.loop, col.loop],
               y1 = counts[[cong.loop]][row.loop, col.loop] + uppers[[cong.loop]][row.loop, col.loop],
               lwd = 2 * zoom,
               angle = 90,
               code = 3,
               length = 0, 
               col = ifelse(theme.loop == 'white', 'black', 'white'))
      }; rm(col.loop)
    }; rm(row.loop)
    
    axis(side = 2,
         at = y.ticks,
         # labels = .y.tick.labels,
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.small * zoom,
         col = ifelse(theme.loop == 'white', 'black', 'white'),
         col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
    add.text.line.multiple.colors(text.segments = rownames(ys),
                                  text.colors = noun.colors[[task.loop]][rownames(ys),'hex'],
                                  .side = 3,
                                  .line = 1,
                                  .outer = TRUE)
    
    dev.off()  
    
    
    ##
    ## Just IFG and MFG
    ##
    
    pdf(paste0(bar.plot.dir.rois.stacked, 'count of classifiers sig - ',model.loop,' - n classifiers total=',unique(unlist(ys.denom)),' - just prefrontal - ',cong.loop,'.pdf'),
        width = 6, height = 5.5)
    par(oma = c(0,0,0,0),
        mar = c(1,3,1,3) * zoom,
        mfrow = c(1, 1),
        lwd = 2 * zoom,
        bg = rgb(1,1,1,0))
    
    xs <- barplot(counts[[cong.loop]],
                  ylim = y.limits,
                  col = c(rgb(.25,.25,.25), rgb(.75,.75,.75)),
                  beside = TRUE,
                  border = theme.loop,
                  lwd = 2 * zoom,
                  cex.names = text.size.small * zoom,
                  yaxt = 'n',
                  las = 1)
    
    top.vals <- counts[[cong.loop]]
    for(roi.loop in names(roi.counts[[cong.loop]])){
      # roi.loop = names(roi.counts[[cong.loop]])[1]
      
      current.cols <- rep(rgb(1,1,1,0), times = 2)
      if(roi.loop == 'MFG'){current.cols <- rep(roi.colors$rostral.caudal[roi.loop], times = 2)}
      if(roi.loop == 'IFG'){current.cols <- rep(roi.colors$rostral.caudal[roi.loop], times = 2)}
      
      barplot(top.vals,
              # col = c(rgb(.25,.25,.25), rgb(.75,.75,.75)),
              # col = c(roi.colors$rostral.caudal[roi.loop],
              #         colorspace::darken(roi.colors$rostral.caudal[roi.loop], 
              #                            amount = .1, method = "relative")),
              col = current.cols,
              beside = TRUE,
              border = FALSE,
              lwd = 2 * zoom,
              # cex.names = text.size.small * zoom,
              yaxt = 'n',
              las = 1,
              add = TRUE)
      
      top.vals <- top.vals - roi.counts[[cong.loop]][[roi.loop]]
      
    }; rm(roi.loop)
    
    # Replot original to overlay borders
    barplot(counts[[cong.loop]],
            ylim = y.limits,
            col = rgb(1,1,1,0),
            beside = TRUE,
            border = theme.loop,
            lwd = 2 * zoom,
            cex.names = text.size.small * zoom,
            yaxt = 'n',
            las = 1,
            add = TRUE)
    
    # Error bars
    for(row.loop in 1:nrow(xs)){
      for(col.loop in 1:ncol(xs)){
        arrows(x0 = xs[row.loop, col.loop],
               y0 = counts[[cong.loop]][row.loop, col.loop] - lowers[[cong.loop]][row.loop, col.loop],
               y1 = counts[[cong.loop]][row.loop, col.loop] + uppers[[cong.loop]][row.loop, col.loop],
               lwd = 2 * zoom,
               angle = 90,
               code = 3,
               length = 0, 
               col = ifelse(theme.loop == 'white', 'black', 'white'))
      }; rm(col.loop)
    }; rm(row.loop)
    
    axis(side = 2,
         at = y.ticks,
         # labels = .y.tick.labels,
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.small * zoom,
         col = ifelse(theme.loop == 'white', 'black', 'white'),
         col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
    add.text.line.multiple.colors(text.segments = rownames(ys),
                                  text.colors = noun.colors[[task.loop]][rownames(ys),'hex'],
                                  .side = 3,
                                  .line = 1,
                                  .outer = TRUE)
    
    dev.off()  
    
    
    ##
    ## Pie charts - just MFG and IFG
    ##
    
    # Save dir
    pie.plot.dir.rois.stacked <- paste0(output.path, 
                                        'figures/15a - barplots - proportion classifiers sig/',
                                        theme.loop,
                                        '/all patients - roi pie charts - by congruent noun - counts/')
    dir.create(pie.plot.dir.rois.stacked, showWarnings = FALSE, recursive = TRUE)
    
    
    pie.cols <- c()
    for(roi.loop in names(roi.counts[[cong.loop]])){
      if(! roi.loop %in% c('IFG','MFG')){
        pie.cols <- c(pie.cols,
                      c('active' = rgb(.25,.25,.25), 'passive' = rgb(.75,.75,.75))['passive'])}
      if(roi.loop == 'MFG'){pie.cols <- c(pie.cols, roi.colors$rostral.caudal[roi.loop])}
      if(roi.loop == 'IFG'){pie.cols <- c(pie.cols, roi.colors$rostral.caudal[roi.loop])}  
    }#; rm(roi.loop)
    
    for(syntax.loop in c('active','passive')){
      # syntax.loop = 'active'
      
      for(role.loop in c(1,2)){
        # role.loop = 2
        
        # P-values
        temp.ps <- t(prop.stats[[cong.loop]][[syntax.loop]])[role.loop,]
        temp.ps <- paste(paste0(names(temp.ps), '=', round(temp.ps, 3)), collapse = '_')
        
        ### All rois
        pdf(paste0(pie.plot.dir.rois.stacked, 'count of classifiers sig - ',model.loop,' - n classifiers total=',unique(unlist(ys.denom)),' - just prefrontal - ',cong.loop,' - ',syntax.loop,' - ',c('during subject','during object')[role.loop],'_',temp.ps,'.pdf'),
            width = 5, height = 5)
        
        par(bg = rgb(1,1,1,0))
        pie(sapply(roi.counts[[cong.loop]], function(.roi){.roi['passive',role.loop]}),
            col = roi.colors$rostral.caudal[names(roi.counts[[cong.loop]])],
            # col = pie.cols,
            border = NA,
            labels = '')  
        
        dev.off()
        
        ### Outline of IFG and MFG
        for(roi.loop in c('IFG','MFG')){
          # roi.loop = 'IFG'
          
          current.roi.counts <- sapply(roi.counts[[cong.loop]], function(.roi){unname(.roi[syntax.loop, role.loop])})
          current.roi.counts <- c(sum(current.roi.counts) - current.roi.counts[roi.loop], 
                                  current.roi.counts[roi.loop])
          
          if(sum(current.roi.counts)){
        pdf(paste0(pie.plot.dir.rois.stacked, 'count of classifiers sig - ',model.loop,' - n classifiers total=',unique(unlist(ys.denom)),' - just prefrontal - ',cong.loop,' - ',syntax.loop,' - ',c('during subject','during object')[role.loop],' - ',roi.loop,' border.pdf'),
            width = 5, height = 5)
        
        # Add boarder around prefrontal regions
        par(bg = rgb(1,1,1,0),
            lwd = zoom * 8)
        
        pie(current.roi.counts,
            col = rgb(1,1,1,0),
            border = c(rgb(1,1,1,0), ifelse(theme.loop == 'white','black','white')),
            labels = '')
        
        dev.off()
          }# if(sum(current.roi.counts)){
        
        }#; rm(roi.loop)
      }#; rm(role.loop)
      
    }; rm(syntax.loop)
    
  }#; rm(cong.loop)
  
}#; rm(theme.loop)



##
## Patient stats, collapsing across ROI - grouped by syntactic role - COUNT
##

zoom <- 1.8

counts <- list()
uppers <- list()
lowers <- list()

counts[['subject']] <- ld2a(lapply(ys, as.data.frame))[,1,]
counts[['object']] <- ld2a(lapply(ys, as.data.frame))[,2,]

uppers[['subject']] <- ld2a(lapply(current.errors.upper, as.data.frame))[,1,]
uppers[['object']] <- ld2a(lapply(current.errors.upper, as.data.frame))[,2,]

lowers[['subject']] <- ld2a(lapply(current.errors.lower, as.data.frame))[,1,]
lowers[['object']] <- ld2a(lapply(current.errors.lower, as.data.frame))[,2,]


for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  # Save dir
  bar.plot.dir.rois <- paste0(output.path, 
                              'figures/15a - barplots - proportion classifiers sig/',
                              theme.loop,
                              '/all patients - collapsed rois - by role - counts/')
  dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)
  
  for(role.loop in names(counts)){
    # role.loop = names(counts)[1]
    
    pdf(paste0(bar.plot.dir.rois, 'proportion of classifiers sig - ',model.loop,' - n classifiers total=',unique(unlist(ys.denom)),' - ',role.loop,'.pdf'),
        width = 6, height = 5.5)
    par(oma = c(0,0,0,0),
        mar = c(1,3,1,3) * zoom,
        mfrow = c(1, 1),
        lwd = 2 * zoom,
        bg = rgb(1,1,1,0))
    
    xs <- barplot(counts[[role.loop]],
                  ylim = y.limits,
                  # col = c(rgb(.25,.25,.25), rgb(.75,.75,.75)),
                  beside = TRUE,
                  col = noun.colors[['active']][rownames(counts[[role.loop]]),'hex'],
                  density = list('subject' = c(-1, 10, -1, 10),
                                 'object' = c(10, -1, 10, -1))[[role.loop]],
                  # border = theme.loop,
                  border = noun.colors[['active']][rownames(counts[[role.loop]]),'hex'],
                  lwd = 2 * zoom,
                  cex.names = text.size.small * zoom,
                  yaxt = 'n',
                  names.arg = rep('', 2),
                  space = c(0.15, 1),
                  las = 1)
    # Error bars
    for(row.loop in 1:nrow(xs)){
      for(col.loop in 1:ncol(xs)){
        arrows(x0 = xs[row.loop, col.loop],
               y0 = counts[[role.loop]][row.loop, col.loop] - lowers[[role.loop]][row.loop, col.loop],
               y1 = counts[[role.loop]][row.loop, col.loop] + uppers[[role.loop]][row.loop, col.loop],
               lwd = 2 * zoom,
               angle = 90,
               code = 3,
               length = 0, 
               col = rgb(.5,.5,.5))
      }; rm(col.loop)
    }; rm(row.loop)
    
    axis(side = 2,
         at = y.ticks,
         # labels = .y.tick.labels,
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.small * zoom,
         col = ifelse(theme.loop == 'white', 'black', 'white'),
         col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
    
    dev.off()  
    
  }#; rm(role.loop)
  
}#; rm(theme.loop)

# zoom <- 1



##
## Combined stats
##

for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  # Save dir
  bar.plot.dir.rois <- paste0(output.path, 
                              'figures/15a - barplots - proportion classifiers sig/',
                              theme.loop,
                              '/average across patients - by roi/')
  dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)
  
  for(roi.loop in rois){
    # roi.loop = rois[3]
    
    pdf(paste0(bar.plot.dir.rois, 'proportion of classifiers sig - ',roi.loop,' - ',model.loop,'.pdf'),
        width = 6 * 4, height = 6)
    par(oma = c(0,0,2,0),
        mar = c(3,3,1,1) * zoom,
        mfrow = c(1, 4),
        lwd = 2 * zoom,
        bg = rgb(1,1,1,0))
    
    for(task.loop in c('active','passive','list','sentence')){
      # task.loop = names(task.grids)[1]
      
      ys <- combined.sig.proportions[[task.loop]][[roi.loop]]
      current.errors.upper <- combined.sig.proportions.95ci.upper[[task.loop]][[roi.loop]]
      current.errors.lower <- combined.sig.proportions.95ci.lower[[task.loop]][[roi.loop]]
      
      xs <- barplot(ys,
                    ylim = c(0, .5),
                    beside = TRUE,
                    col = noun.colors[[task.loop]][rownames(ys),'hex'],
                    border = theme.loop,
                    cex.names = text.size.small * zoom,
                    yaxt = 'n',
                    las = 1)
      # Error bars
      for(row.loop in 1:nrow(xs)){
        for(col.loop in 1:ncol(xs)){
          arrows(x0 = xs[row.loop, col.loop],
                 y0 = ys[row.loop, col.loop] - current.errors.lower[row.loop, col.loop],
                 y1 = ys[row.loop, col.loop] + current.errors.upper[row.loop, col.loop],
                 angle = 90,
                 code = 3,
                 length = 0, 
                 col = rgb(.7,.7,.7))
        }; rm(col.loop)
      }; rm(row.loop)
      
      axis(side = 2,
           at = c(0, .5),
           # labels = .y.tick.labels,
           las = 0,
           tck = -.025 * zoom, # length of tick
           padj = -.45 * zoom, # distance between tick and label
           lwd = 1.5 * zoom,
           lwd.ticks = 1.5 * zoom,
           cex.axis = text.size.small * zoom,
           col = ifelse(theme.loop == 'white', 'black', 'white'),
           col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
      add.text.line.multiple.colors(text.segments = rownames(ys),
                                    text.colors = noun.colors[[task.loop]][rownames(ys),'hex'],
                                    .side = 3,
                                    .line = 1,
                                    .outer = TRUE)
      
    }#; rm(task.loop)
    
    dev.off()  
    
  }#; rm(roi.loop)
  
}#; rm(theme.loop)


##
## Combined stats, collapsing across ROI
##

for(theme.loop in c('black','white')){
  # theme.loop = c('black','white')[2]
  
  # Save dir
  bar.plot.dir.rois <- paste0(output.path, 
                              'figures/15a - barplots - proportion classifiers sig/',
                              theme.loop,
                              '/average across patients - collapsed rois/')
  dir.create(bar.plot.dir.rois, showWarnings = FALSE, recursive = TRUE)
  
  pdf(paste0(bar.plot.dir.rois, 'proportion of classifiers sig - ',model.loop,'.pdf'),
      width = 6 * 4, height = 6)
  par(oma = c(0,0,2,0),
      mar = c(3,3,1,1) * zoom,
      mfrow = c(1, 4),
      lwd = 2 * zoom,
      bg = rgb(1,1,1,0))
  
  for(task.loop in c('active','passive','list','sentence')){
    # task.loop = names(task.grids)[2]
    
    ys.sum <- elementwise.matrix.apply(combined.sig.counts[[task.loop]], .function = "sum")
    ys.denom <- elementwise.matrix.apply(combined.n.classifiers[[task.loop]], .function = "sum")
    ys <- ys.sum / ys.denom
    current.errors <- binom.confint(x = as.vector(ys.sum),
                                    n = as.vector(ys.denom),
                                    conf.level = .95, 
                                    methods = "wilson")
    current.errors.upper <- matrix(data = current.errors$upper, nrow = 2) - ys
    current.errors.lower <- ys - matrix(data = current.errors$lower, nrow = 2)
    
    xs <- barplot(ys,
                  ylim = c(0, .3),
                  beside = TRUE,
                  col = noun.colors[[task.loop]][rownames(ys),'hex'],
                  border = theme.loop,
                  cex.names = text.size.small * zoom,
                  yaxt = 'n',
                  las = 1)
    # Error bars
    for(row.loop in 1:nrow(xs)){
      for(col.loop in 1:ncol(xs)){
        arrows(x0 = xs[row.loop, col.loop],
               y0 = ys[row.loop, col.loop] - current.errors.lower[row.loop, col.loop],
               y1 = ys[row.loop, col.loop] + current.errors.upper[row.loop, col.loop],
               angle = 90,
               code = 3,
               length = 0, 
               col = rgb(.7,.7,.7))
      }; rm(col.loop)
    }; rm(row.loop)
    
    axis(side = 2,
         at = c(0, .3),
         # labels = .y.tick.labels,
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.small * zoom,
         col = ifelse(theme.loop == 'white', 'black', 'white'),
         col.axis = ifelse(theme.loop == 'white', 'black', 'white'))
    add.text.line.multiple.colors(text.segments = rownames(ys),
                                  text.colors = noun.colors[[task.loop]][rownames(ys),'hex'],
                                  .side = 3,
                                  .line = 1,
                                  .outer = TRUE)
    
  }#; rm(task.loop)
  
  dev.off()  
  
}#; rm(theme.loop)






# }#; rm(model.loop)
# }; rm(band.loop)








# Finish!
message('Script completed successfully. ',Sys.time())


