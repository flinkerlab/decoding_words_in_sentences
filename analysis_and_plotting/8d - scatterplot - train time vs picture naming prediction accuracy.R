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
model.types <- c('multinom','nnet')
for(model.loop in model.types){ ## UNCOMMENT
# model.loop = model.types[1] ## UNCOMMENT


###
### Get data: Peak accuracy vs. central train time
###

template <- data.frame('roi' = NA,
                       'stage' = NA,
                       'patient' = NA,
                       'mean.train.time' = NA,
                       'peak.accuracy' = NA)
plot.data <- list()

### Loop thru rois
for(roi.loop in rev(rois)){ ## UNCOMMENT
  # roi.loop = rois[3] ## UNCOMMENT
  
  ### Loop thru stages
  # For clean up
  stage.keep <- c(ls(), 'stage.keep', 'stage.loop')
  
  # Loop!
  for(stage.loop in rev(names(stage.sample.labels))){ ## UNCOMMENT
    # stage.loop = names(stage.sample.labels)[12] ## UNCOMMENT
    
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
    
    
    for(patient in names(patient.pn.stats)){
      # patient = names(patient.pn.stats)[1]
      
      ## If any sig windows
      if(nrow(patient.pn.sig.windows[[patient]]) > 0){
        
        ## Get windows to draw peak from
        .sig.samples <- c()
        for(i in 1:nrow(patient.pn.sig.windows[[patient]])){
          .start.sample <- time.convert(patient.pn.sig.windows[[patient]]$start.time[i], 'times', 'samples', 256)
          .end.sample <- time.convert(patient.pn.sig.windows[[patient]]$end.time[i], 'times', 'samples', 256)
          .sig.samples <- c(.sig.samples,
                            time.convert(.start.sample:.end.sample, 'samples', 'sample.labels', 256))
          rm(.start.sample, .end.sample)
        }; rm(i)
        
        ## Get data
        .index <- length(plot.data) + 1
        plot.data[[.index]] <- template
        plot.data[[.index]]$roi <- roi.loop
        plot.data[[.index]]$stage <- stage.loop
        plot.data[[.index]]$mean.train.time <- mean(unlist(stage.ranges[stage.loop, c('start.time','end.time')]))
        plot.data[[.index]]$start.train.time <- stage.ranges[stage.loop, 'start.time']
        plot.data[[.index]]$patient <- patient
        plot.data[[.index]]$peak.accuracy <- max(patient.pn.stats[[patient]][.sig.samples,'accuracy'])
        plot.data[[.index]]$peak.accuracy.time <- time.convert(.sig.samples[which.max(patient.pn.stats[[patient]][.sig.samples,'accuracy'])[1]], "sample.labels", "times", 256)
        plot.data[[.index]]$earliest.sig.time <- min(time.convert(.sig.samples, 'sample.labels', 'times', 256))
        plot.data[[.index]]$roi.color <- roi.colors$rostral.caudal[roi.loop]
        
      } # if(nrow(patient.pn.sig.windows[[patient]]) > 0){
      
    }; rm(patient)
  }; rm(stage.loop)
}; rm(roi.loop)
message('Done getting peaks! ',Sys.time())## < 1 min


### Clean up
plot.data <- bind_rows(plot.data)


###
save.plot.dir.main <- paste0(output.path, 'figures/8d - scatterplot - train time vs picture naming prediction accuracy/')
zoom <- 1.8
text.size.big <- 1.8
text.size.med <- 1.6
text.size.small <- 1.4


##
## Plot: Accuracy
##

for(theme.loop in c('white','black')){ ## UNCOMMENT
  # theme.loop = c('white','black')[1] ## UNCOMMENT
  
  for(dv.loop in c('peak.accuracy', 'peak.accuracy.time', 'earliest.sig.time')){
    # dv.loop = c('peak.accuracy', 'peak.accuracy.time', 'earliest.sig.time')[2]
    
    for(color.loop in c("roi colored","purple")){
      # color.loop = c("roi colored","purple")[1]
      

    
    # Save dir
    plot.dir.same.axes.square <- paste0(save.plot.dir.main,'same axes - square/',theme.loop,'/')
    dir.create(plot.dir.same.axes.square, showWarnings = FALSE, recursive = TRUE)
    
    # Model
    model <- lm(plot.data[,dv.loop] ~ plot.data$mean.train.time)
    current.p <- summary(model)$coefficients[2,4]
    
    # Y limits
    if(dv.loop == 'peak.accuracy'){
      xs <- plot.data$mean.train.time + 
        sample(-15:15, size = nrow(plot.data), replace = TRUE) # jitter train times to reduce overlap/increase visibility
      y.range <- range(plot.data[,dv.loop])
      x.range <- range(plot.data$mean.train.time)
    }else{
      xs <- plot.data$start.train.time + 
        sample(-15:15, size = nrow(plot.data), replace = TRUE) # jitter train times to reduce overlap/increase visibility
      x.range <- range(xs)
      y.range <- range(plot.data[,dv.loop])
      # Make them the same
      x.range <- y.range <- range(c(x.range, y.range))
    }
    
    # Colors
    if(color.loop == "roi colored"){
      current.colors <- adjust.transparency(plot.data$roi.color, .8)
    }else{
      current.colors <- adjust.transparency(colors[[theme.loop]]$rainbow['purple','hex'], .5)
    }
    
    
    # Plot
    pdf(paste0(plot.dir.same.axes.square, dv.loop, ' - ', color.loop,' - p=',round(current.p, 4),'.pdf'),
        width = 5, height = 5)
    par(bg = rgb(1,1,1,0),
        pty = 's',
        # mar = c(6,6,3,1) * zoom,
        mar = c(3,3,0,1) * zoom,
        oma = c(0,0,1,0))
    plot(x = xs,
         y = plot.data[,dv.loop],
         xlab = '',
         xlim = x.range,
         ylim = y.range,
         ylab = '',
         pch = 20,
         cex = 1.4 * zoom,
         col = current.colors,
         bty = 'n',
         xaxt = 'n',
         yaxt = 'n')
    abline(model,
           lwd = 3 * zoom)
    if(dv.loop != 'peak.accuracy'){
      # abline(h = 0,
             # lty = 2, 
             # lwd = 1.5 * zoom)
      abline(a = 0, 
             b = 1, 
             lty = 2,
             lwd = 1.5 * zoom)
    }else{
      abline(v = 0,
             lty = 2,
             lwd = 1.5 * zoom)
    }
    ## Axes
    
    if(dv.loop == "peak.accuracy"){
      x.range <- c(ceiling(min(plot.data$mean.train.time) / 50) * 50,
                   floor(max(plot.data$mean.train.time) / 50) * 50)
      y.range <- c(ceiling(min(plot.data[,dv.loop]) * 20) / 20,
                   floor(max(plot.data[,dv.loop]) * 20) / 20)
      if(y.range[1] == y.range[2]){
        y.range <- c(y.range[1] - .5,
                     y.range[2] + .5)
      } # if(y.range[1] == y.range[2]){
    }else{ # if(dv.loop == "peak.accuracy"){
      y.range <- x.range <- c(ceiling(min(x.range) / 50) * 50,
                   floor(max(x.range) / 50) * 50)
    } # if(dv.loop == "peak.accuracy"){}else{
    axis(side = 1,
         at = sort(c(x.range, ifelse(dv.loop == "peak.accuracy", numeric(0), 0))),
         labels = sort(c(x.range, ifelse(dv.loop == "peak.accuracy", numeric(0), 0))),
         las = 1,
         tck = -.025 * zoom, # length of tick
         padj = .6 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.med * zoom,
         col = ifelse(theme.loop == 'black', 'white', 'black'),
         col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
    axis(side = 2,
         at = sort(c(y.range, ifelse(dv.loop == "peak.accuracy", numeric(0), 0))),
         labels = sort(c(y.range, ifelse(dv.loop == "peak.accuracy", numeric(0), 0))),
         las = 0,
         tck = -.025 * zoom, # length of tick
         padj = -.45 * zoom, # distance between tick and label
         lwd = 1.5 * zoom,
         lwd.ticks = 1.5 * zoom,
         cex.axis = text.size.med * zoom,
         col = ifelse(theme.loop == 'black', 'white', 'black'),
         col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
    
    if(color.loop == "roi colored"){
      add.text.line.multiple.colors(rois, text.colors = roi.colors$rostral.caudal[rois], .side = 1, .line = 5)
    }
    
    dev.off()
  
    }#; rm(color.loop)  
    }; rm(dv.loop)
  }#; rm(theme.loop)
  
  

##
## By ROI
##

for(theme.loop in c('white','black')){ ## UNCOMMENT
  # theme.loop = c('white','black')[1] ## UNCOMMENT
  
  for(dv.loop in c('peak.accuracy', 'peak.accuracy.time', 'earliest.sig.time')){
    # dv.loop = c('peak.accuracy', 'peak.accuracy.time', 'earliest.sig.time')[2]
    
    for(color.loop in c("roi colored","purple")){
      # color.loop = c("roi colored","purple")[1]
      
      for(roi.loop in rois){
      
      # Current data
      current.plot.data <- plot.data[plot.data$roi == roi.loop,]
      
      # Save dir
      plot.dir.same.axes.square <- paste0(save.plot.dir.main,'same axes - square - by ROI/',theme.loop,'/',dv.loop,'/',color.loop,'/')
      dir.create(plot.dir.same.axes.square, showWarnings = FALSE, recursive = TRUE)
      
      # Model
      model <- lm(current.plot.data[,dv.loop] ~ current.plot.data$mean.train.time)
      current.p <- summary(model)$coefficients[2,4]
      
      # Y limits
      if(dv.loop == 'peak.accuracy'){
        xs <- current.plot.data$mean.train.time + 
          sample(-15:15, size = nrow(current.plot.data), replace = TRUE) # jitter train times to reduce overlap/increase visibility
        y.range <- range(current.plot.data[,dv.loop])
        x.range <- range(current.plot.data$mean.train.time)
      }else{
        xs <- current.plot.data$start.train.time + 
          sample(-15:15, size = nrow(current.plot.data), replace = TRUE) # jitter train times to reduce overlap/increase visibility
        x.range <- range(xs)
        y.range <- range(current.plot.data[,dv.loop])
        # Make them the same
        x.range <- y.range <- range(c(x.range, y.range))
      }
      
      # Colors
      if(color.loop == "roi colored"){
        current.colors <- adjust.transparency(current.plot.data$roi.color, .8)
      }else{
        current.colors <- adjust.transparency(colors[[theme.loop]]$rainbow['purple','hex'], .5)
      }
      
      
      # Plot
      pdf(paste0(plot.dir.same.axes.square, roi.loop,' - ',dv.loop,' - p=',round(current.p, 4),'.pdf'),
          width = 5, height = 5)
      par(bg = rgb(1,1,1,0),
          pty = 's',
          # mar = c(6,6,3,1) * zoom,
          mar = c(3,3,0,1) * zoom,
          oma = c(0,0,1,0))
      plot(x = xs,
           y = current.plot.data[,dv.loop],
           xlab = '',
           xlim = x.range,
           ylim = y.range,
           ylab = '',
           pch = 20,
           cex = 1.4 * zoom,
           col = current.colors,
           bty = 'n',
           xaxt = 'n',
           yaxt = 'n')
      abline(model,
             lwd = 3 * zoom)
      if(dv.loop != 'peak.accuracy'){
        # abline(h = 0,
        # lty = 2, 
        # lwd = 1.5 * zoom)
        abline(a = 0, 
               b = 1, 
               lty = 2,
               lwd = 1.5 * zoom)
      }else{
        abline(v = 0,
               lty = 2,
               lwd = 1.5 * zoom)
      }
      ## Axes
      
      if(dv.loop == "peak.accuracy"){
        x.range <- c(ceiling(min(current.plot.data$mean.train.time) / 50) * 50,
                     floor(max(current.plot.data$mean.train.time) / 50) * 50)
        y.range <- c(ceiling(min(current.plot.data[,dv.loop]) * 20) / 20,
                     floor(max(current.plot.data[,dv.loop]) * 20) / 20)
        if(y.range[1] == y.range[2]){
          y.range <- c(y.range[1] - .5,
                       y.range[2] + .5)
        } # if(y.range[1] == y.range[2]){
      }else{ # if(dv.loop == "peak.accuracy"){
        y.range <- x.range <- c(ceiling(min(x.range) / 50) * 50,
                                floor(max(x.range) / 50) * 50)
      } # if(dv.loop == "peak.accuracy"){}else{
      axis(side = 1,
           at = sort(c(x.range, ifelse(dv.loop == "peak.accuracy", numeric(0), 0))),
           labels = sort(c(x.range, ifelse(dv.loop == "peak.accuracy", numeric(0), 0))),
           las = 1,
           tck = -.025 * zoom, # length of tick
           padj = .6 * zoom, # distance between tick and label
           lwd = 1.5 * zoom,
           lwd.ticks = 1.5 * zoom,
           cex.axis = text.size.med * zoom,
           col = ifelse(theme.loop == 'black', 'white', 'black'),
           col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
      axis(side = 2,
           at = sort(c(y.range, ifelse(dv.loop == "peak.accuracy", numeric(0), 0))),
           labels = sort(c(y.range, ifelse(dv.loop == "peak.accuracy", numeric(0), 0))),
           las = 0,
           tck = -.025 * zoom, # length of tick
           padj = -.45 * zoom, # distance between tick and label
           lwd = 1.5 * zoom,
           lwd.ticks = 1.5 * zoom,
           cex.axis = text.size.med * zoom,
           col = ifelse(theme.loop == 'black', 'white', 'black'),
           col.axis = ifelse(theme.loop == 'black', 'white', 'black'))
      
      if(color.loop == "roi colored"){
        add.text.line.multiple.colors(rois, text.colors = roi.colors$rostral.caudal[rois], .side = 1, .line = 5)
      }
      
      dev.off()
  
      }#; rm(roi.loop)
    }#; rm(color.loop)  
  }; rm(dv.loop)
}#; rm(theme.loop)


  
  }; rm(model.loop)
  # }; rm(band.loop)
  
  
  
  
  
  
  
  
  # Finish!
  message('Script completed successfully. ',Sys.time())
  
  
  