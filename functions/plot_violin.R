### Plot time series - general function
### Fall 2023 (!!!!!!)
### adam.milton.morgan@gmail.com

### TO DO:
# - specify distribution type ("gamma"? for each column as well as column-specific lower limits for gamma distributions)
# - fix error bars: allow for choice of SD vs. SD vs. manual;
# - fix error bars: allow for different lengths for upper and lower bars
# - allow for manual transparency adjustment
# - add subcolumn labels on bottom axis; change super column labels to top axis?
# - add option to multiply jitter vals by density height so the dots have more of a bell curve distribution
# - add option for manual points (these are noise distributions, want real data points too)

plot.violin <- function(
    .values,
    .y.limits = NULL,
    .x.limits = NULL,
    .starting.x.val = NULL,
    .y.label = '',
    .x.label = '',
    .title = '',
    .colors = NULL,
    .dot.colors = NULL,
    .violin.colors = NULL,
    .violin.alpha = NULL,
    .dot.alpha = NULL,
    .jitter = TRUE,
    .jitter.mean = 0,
    .jitter.width = 1,
    .scale.jitter.by.density = TRUE,
    .zoom = 1.2,
    .point.size = 1.6,
    .mean.point.size = 1.6,
    .mean.pch = 20,
    .mean.colors = NULL,
    .show.violin = TRUE,
    .smoothing.factor = 3,
    .show.dots = TRUE,
    .show.mean = TRUE,
    .show.bars = FALSE,
    .bar.spacing = .2,
    .bar.colors = NULL,
    .bar.border.colors = NA,
    # .show.se = FALSE,
    # .show.sd = TRUE,
    .error.bar.type = c('sd','se')[1],
    .error.bar.caps = FALSE,
    .error.bar.colors = NULL,
    .multiply.error.bars.by = NULL, # e.g., qnorm(.975)
    .error.bar.upper.length = NULL,
    .error.bar.lower.length = NULL,
    .n.sds.omit.outliers = 3,
    .cutoff.upper.val = NULL,
    .cutoff.lower.val = NULL,
    .cutoff.upper.percentile = NULL,
    .cutoff.lower.percentile = NULL,
    .trim.base.to.remove.values.greater.than = NA,
    .trim.base.to.remove.values.less.than = NA,
    .violin.side = c('left','right','both')[1],
    .dots.side = c('right','left','both')[1],
    .behind = c("violins", "dots")[1],
    .extra.points = NULL,
    .extra.point.pch = 18,
    .extra.point.size = 2,
    .extra.point.color = rgb(.8,.8,.1),
    .show.x.axis = TRUE,
    .show.y.axis = TRUE,
    .y.ticks = NULL,
    .y.tick.labels = NULL,
    .y.axis.color = NULL,
    .x.ticks = NULL,
    .x.tick.labels = NULL,
    .x.tick.labels.vertical = FALSE,
    .horizontal.line.at = NA,
    .theme = 'white',
    .background = NULL,
    .margin = NULL,
    .midline.offset = .1,
    .max.density.y = .35,
    .distribution = c('normal','gamma')[1], # if the data are gamma distributed, get rid of taper at lower end
    .gamma.dist.lower.limit = NULL
){
  
  ### Lemmas and packages
  source(paste0(path,'/analysis/R/functions/min_max.R'))
  source(paste0(path,'/analysis/R/functions/adjust_transparency.R'))
  source(paste0(path,'/analysis/R/functions/get_density.R'))
  library(pals)
  
  
  # Text sizes
  text.size.big <- 1.8
  text.size.med <- 1.6
  text.size.small <- 1.4
  
  # Adjust size stuff
  .midline.offset <- .midline.offset * .zoom
  
  if(class(.values) != "list"){message("Input '.values' must be a list of vectors or a list of lists of vectors.")}
  
  # Make it a list of lists if it's not already
  if(class(.values[[1]]) != "list"){.values <- list(.values)}
  if(! is.null(.extra.points)){
    if(class(.extra.points[[1]]) != "list"){.extra.points <- list(.extra.points)}  
  }
  
  # Number of bins
  n.super.bins <- length(.values)
  n.sub.bins <- length(.values[[1]])
  
  ### Variable names
  # Outer loop ("super")
  if(! is.null(names(.values))){
    super.labels <- names(.values)
  }else{
    if(length(.values) > 1){print('List entries not labeled, using numbers.')}
    super.labels <- paste0("super_", 1:length(.values))
    names(.values) <- super.labels
  }
  # Inner loop ("sub")
  if(! is.null(names(.values[[1]]))){
    sub.labels <- names(.values[[1]])
  }else{
    print('List sub-entries not labeled, using numbers.')
    sub.labels <- paste0("sub_", 1:n.sub.bins)
    for(i in 1:n.super.bins){names(.values[[i]]) <- sub.labels}; rm(i)
  }
  
  ### Assign these names to .extra.points too
  if(! is.null(.extra.points)){
    # Outer loop ("super")
    names(.extra.points) <- super.labels
    # Inner loop ("sub")
    for(i in 1:n.sub.bins){
      # i = 1
      .extra.points <- 
        lapply(.extra.points, function(x){
          names(x)[i] <- sub.labels[i]
          return(x)
        })
    }; rm(i)
  } # if(! is.null(.extra.points)){
  
  ### Get colors
  if(is.null(.dot.colors)){
    if(!is.null(.colors)){
      .dot.colors <- .colors
    }else{
      colors.df <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',.theme,'/rainbow_bright.csv'), row.names = 1)
      .dot.colors <- cubicl(n.sub.bins)  
    }
  }
  if(length(.dot.colors) == 1){.dot.colors <- rep(.dot.colors, times = n.sub.bins)}
  if(is.null(names(.dot.colors))){
    names(.dot.colors) <- sub.labels  
  }
  if(is.null(.violin.colors)){
    .violin.colors <- .dot.colors
  }
  if(length(.violin.colors) == 1){.violin.colors <- rep(.violin.colors, times = n.sub.bins)}
  if(is.null(names(.violin.colors))){
    names(.violin.colors) <- sub.labels  
  }
  
  # Error bar colors
  if(is.null(.error.bar.colors)){
    if(!is.null(.colors)){
      .error.bar.colors <- .colors
    }else{
      colors.df <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',.theme,'/rainbow_bright.csv'), row.names = 1)
      .error.bar.colors <- cubicl(n.sub.bins)  
    }
  }
  if(length(.error.bar.colors) == 1){.error.bar.colors <- rep(.error.bar.colors, times = n.sub.bins)}
  if(is.null(names(.error.bar.colors))){
    names(.error.bar.colors) <- sub.labels  
  }
  
  # Mean colors
  # if(.show.mean){
    if(is.null(.mean.colors)){
      if(!is.null(.colors)){
        .mean.colors <- .colors
      }else{
        colors.df <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',.theme,'/rainbow_bright.csv'), row.names = 1)
        .mean.colors <- cubicl(n.sub.bins)  
      }
    }
    if(length(.mean.colors) == 1){.mean.colors <- rep(.mean.colors, times = n.sub.bins)}
    if(is.null(names(.mean.colors))){
      names(.mean.colors) <- sub.labels  
    }
  # }
  
  # Error bar colors
  if(is.null(.error.bar.colors)){
    if(!is.null(.mean.colors)){
      .error.bar.colors <- .mean.colors
    }else{
      colors.df <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',.theme,'/rainbow_bright.csv'), row.names = 1)
      .error.bar.colors <- cubicl(n.sub.bins)  
    }
  }
  if(length(.error.bar.colors) == 1){.error.bar.colors <- rep(.error.bar.colors, times = n.sub.bins)}
  if(is.null(names(.error.bar.colors))){
    names(.error.bar.colors) <- sub.labels  
  }
  
  ## Barplot colors
  # Bar fill colors
  if(is.null(.bar.colors)){
    if(!is.null(.colors)){
      .bar.colors <- .colors
    }else{
      colors.df <- read.csv(paste0(path,'analysis/R/color palettes/output/theme_',.theme,'/rainbow_bright.csv'), row.names = 1)
      .bar.colors <- cubicl(n.sub.bins)  
    }
  }
  if(length(.bar.colors) == 1){.bar.colors <- rep(.bar.colors, times = n.sub.bins)}
  if(is.null(names(.bar.colors))){
    names(.bar.colors) <- sub.labels  
  }
  # Bar border colors
  if(length(.bar.border.colors) == 1){.bar.border.colors <- rep(.bar.border.colors, times = n.sub.bins)}
  if(is.null(names(.bar.border.colors))){
    names(.bar.border.colors) <- sub.labels  
  }
  
  ## Axis colors
  if(is.null(.y.axis.color)){
    .y.axis.color <- ifelse(.theme == 'black', 'white', 'black')
  }
  if(is.null(.y.axis.color)){
    .y2.axis.color <- ifelse(.theme == 'black', 'white', 'black')
  }
  
  ## Set up error bars
  if(is.na(.error.bar.type)){
    .error.bar.type <- ''
  }
  
  ## Dataframe of plot metadata
  plot.df <- data.frame(x = NA,
                        mean = NA,
                        sd = NA,
                        se = NA,
                        bootstrap.ci = NA,
                        sub.label = NA,
                        super.label = NA,
                        extra.points = NA)[0,]
  if(is.null(.starting.x.val)){
    x <- ifelse(length(.values) > 1, -.5, -1.5)  
  }else{
    x <- .starting.x.val - 1
  }
  
  row <- 0
  for(super.loop in 1:n.super.bins){
    # super.loop = 1
    x <- x + 1
    for(sub.loop in 1:n.sub.bins){
      # sub.loop = 1
      x <- x + 1 
      row <- row + 1
      plot.df[row, 'x'] <- x
      if(length(.values[[super.loop]][[sub.loop]]) > 0){
        plot.df[row, 'mean'] <- mean(.values[[super.loop]][[sub.loop]], na.rm = TRUE)
        plot.df[row, 'color'] <- .dot.colors[sub.labels[sub.loop]]
        plot.df[row, 'error.bar.color'] <- .error.bar.colors[sub.labels[sub.loop]]
        plot.df[row, 'bar.color'] <- .bar.colors[sub.labels[sub.loop]]
        plot.df[row, 'mean.color'] <- .mean.colors[sub.labels[sub.loop]]
        plot.df[row, 'bar.border.color'] <- .bar.border.colors[sub.labels[sub.loop]]
        plot.df[row, 'sd'] <- sd(.values[[super.loop]][[sub.loop]], na.rm = TRUE)
        plot.df[row, 'se'] <- plot.df$sd[row] / sqrt(length(.values[[super.loop]][[sub.loop]][! is.na(.values[[super.loop]][[sub.loop]])]))
        plot.df[row, 'super.label'] <- super.labels[super.loop]
        plot.df[row, 'sub.label'] <- sub.labels[sub.loop]
        plot.df[row, 'bootstrap.ci'] <- sd(replicate(1000,
                                                     mean(sample(.values[[super.loop]][[sub.loop]],
                                                                 size = length(.values[[super.loop]][[sub.loop]]),
                                                                 replace = TRUE)))) * qnorm(.975)  
        if(! is.null(.extra.points)){
          plot.df[row, 'extra.points'] <- .extra.points[[super.loop]][[sub.loop]]
        } # if(! is.null(.extra.points)){
      } # if any values
    }; rm(sub.loop)
  }; rm(super.loop, row, x)
  
  
  ## Get densities
  densities <- list()
  for(super.loop in super.labels){
    # super.loop = super.labels[1]
    densities[[super.loop]] <- list()
    for(sub.loop in sub.labels){
      # sub.loop = sub.labels[1]
      
      ## Current vals to get density of
      current.density.values <- .values[[super.loop]][[sub.loop]]
      # Get rid of outliers by percentiles first
      if(! is.null(.cutoff.upper.percentile)){
        current.density.values <- current.density.values[current.density.values < quantile(current.density.values, .cutoff.upper.percentile)]
      }
      if(! is.null(.cutoff.lower.percentile)){
        current.density.values <- current.density.values[current.density.values > quantile(current.density.values, .cutoff.lower.percentile)]
      }
      # Get rid of outliers by value
      if(! is.null(.cutoff.upper.val)){
        current.density.values <- current.density.values[current.density.values < .cutoff.upper.val]
      }
      if(! is.null(.cutoff.lower.val)){
        current.density.values <- current.density.values[current.density.values > .cutoff.lower.val]
      }
      # Get rid of outliers by SD
      if(! is.null(.n.sds.omit.outliers)){
        current.mean <- mean(current.density.values)
        current.sd <- sd(current.density.values)
        upper.sd.cutoff <- current.mean + (.n.sds.omit.outliers * current.sd)
        lower.sd.cutoff <- current.mean - (.n.sds.omit.outliers * current.sd)
        current.density.values <- current.density.values[current.density.values < upper.sd.cutoff]
        current.density.values <- current.density.values[current.density.values > lower.sd.cutoff]
        rm(current.mean, current.sd, upper.sd.cutoff, lower.sd.cutoff)
      }
      
      ## Calculate density
      densities[[super.loop]][[sub.loop]] <- 
        get.density(.data = current.density.values,
                    .distribution = .distribution,
                    .smoothing.factor = .smoothing.factor)
      
      ## Trim if specified
      if(! is.na(.trim.base.to.remove.values.greater.than)){
        # Get indices to remove
        remove.indices <- which(densities[[super.loop]][[sub.loop]]$x > 
                                  .trim.base.to.remove.values.greater.than)
        if(length(remove.indices) > 0){
          # Get the index of the x-val at which to start removing distribution
          max.x.index <- min(remove.indices)
          # Get the y-value there
          min.y.val <- densities[[super.loop]][[sub.loop]]$y[max.x.index]
          # Exclude
          densities[[super.loop]][[sub.loop]] <-
            densities[[super.loop]][[sub.loop]][1:max.x.index,]
          # Subtract that y-value and set any resuling negative values to 0
          densities[[super.loop]][[sub.loop]]$y <- 
            densities[[super.loop]][[sub.loop]]$y - min.y.val
          densities[[super.loop]][[sub.loop]]$y[
            densities[[super.loop]][[sub.loop]]$y < 0] <- 0
        }
      } # if(! is.na(.trim.base.to.remove.values.greater.than)){
      if(! is.na(.trim.base.to.remove.values.less.than)){
        # Get indices to remove
        remove.indices <- which(densities[[super.loop]][[sub.loop]]$x < 
                                  .trim.base.to.remove.values.less.than)
        if(length(remove.indices) > 0){
          # Get the index of the x-val at which to start removing distribution
          min.x.index <- max(remove.indices)
          # Get the y-value there
          min.y.val <- densities[[super.loop]][[sub.loop]]$y[min.x.index]
          # Exclude
          densities[[super.loop]][[sub.loop]] <-
            densities[[super.loop]][[sub.loop]][min.x.index:nrow(densities[[super.loop]][[sub.loop]]),]
          # Subtract that y-value and set any resuling negative values to 0
          densities[[super.loop]][[sub.loop]]$y <- 
            densities[[super.loop]][[sub.loop]]$y - min.y.val
          densities[[super.loop]][[sub.loop]]$y[
            densities[[super.loop]][[sub.loop]]$y < 0] <- 0 
        } # if(length(remove.indices) > 0){
      } # if(! is.na(.trim.base.to.remove.values.less.than)){
      
      ## Scale
      # Scale to 1
      densities[[super.loop]][[sub.loop]]$y <- 
        densities[[super.loop]][[sub.loop]]$y / 
        max(densities[[super.loop]][[sub.loop]]$y)
      # Scale to .max.density.y
      densities[[super.loop]][[sub.loop]]$y <- densities[[super.loop]][[sub.loop]]$y * .max.density.y
    }; rm(super.loop)
  }; rm(sub.loop)
  
  ## Jitter .values
  if(.jitter){
    .jitter.values <- list()
    for(super.loop in super.labels){
      # super.loop = super.labels[1]
      .jitter.values[[super.loop]] <- list()
      for(sub.loop in sub.labels){
        # sub.loop = sub.labels[1]
        .jitter.values[[super.loop]][[sub.loop]] <- 
          # rnorm(length(.values[[super.loop]][[sub.loop]]),
          #       mean = .jitter.mean,
          #       sd = .jitter.sd)
          runif(length(.values[[super.loop]][[sub.loop]]),
                min = -.85 * .max.density.y * .jitter.width, # 85% of the maximum of density plot looks good
                max = .85 * .max.density.y * .jitter.width)
        
        if(.scale.jitter.by.density){
          # Find closest x-value
          closest.density.indices <- 
            sapply(.values[[super.loop]][[sub.loop]], function(z){
              which.min(abs(z - densities[[super.loop]][[sub.loop]]$x))[1]})
          # Scale by the density y-value at that x-value
          .jitter.values[[super.loop]][[sub.loop]] <-
            .jitter.values[[super.loop]][[sub.loop]] * 
            (densities[[super.loop]][[sub.loop]]$y[closest.density.indices] / 
               .max.density.y) # unscale from .max.density.y (i.e., scale to 1)
          # Make top and bottom points low
          .jitter.values[[super.loop]][[sub.loop]][c(which.max(.values[[super.loop]][[sub.loop]]),
                                                     which.min(.values[[super.loop]][[sub.loop]]))] <- 
            min(abs(.jitter.values[[super.loop]][[sub.loop]]))
        }
        
        if(.dots.side == 'right'){
          # Make all values positive (for right) and offset a bit to make room for density plot
          .jitter.values[[super.loop]][[sub.loop]] <- 
            abs(.jitter.values[[super.loop]][[sub.loop]]) + .midline.offset
        } # if(.dots.side == 'right'){
        if(.dots.side == 'left'){
          # Make all values positive (for right) and offset a bit to make room for density plot
          .jitter.values[[super.loop]][[sub.loop]] <- 
            -(abs(.jitter.values[[super.loop]][[sub.loop]]) + .midline.offset)
        } # if(.dots.side == 'right'){
      }; rm(super.loop)
    }; rm(sub.loop)
  } # if(.jitter)
  
  ## Add manually specified lower error bar if unspecified
  if(is.null(.error.bar.lower.length) & (! is.null(.error.bar.upper.length))){
    .error.bar.lower.length <- .error.bar.upper.length
  }
  
  ## Don't adjust error bars if not specified
  if(is.null(.multiply.error.bars.by)){.multiply.error.bars.by <- 1}
  
  ## Get x limits
  if(is.null(.x.limits)){
    # .x.limits <- c(0, max(plot.df$x) + ifelse(length(.values) > 1, 1, 0.5))  
    .x.limits <- c(0, max(plot.df$x) + .midline.offset + .max.density.y)  
  } # if(is.null(.x.limits)){
  
  ## Get y limits
  if(is.null(.y.limits)){
    all.y.values <- unlist(.values)
    if(.error.bar.type == 'se'){all.y.values <- 
      c(all.y.values, ((.multiply.error.bars.by * plot.df$se) + plot.df$mean))}
    if(.error.bar.type == 'sd'){all.y.values <- 
      c(all.y.values, ((.multiply.error.bars.by * plot.df$sd) + plot.df$mean))}
    if(! is.null(.error.bar.lower.length)){
      all.y.values <- c(all.y.values, (plot.df$mean - .error.bar.lower.length))}
    if(! is.null(.error.bar.upper.length)){
      all.y.values <- c(all.y.values, (plot.df$mean + .error.bar.upper.length))}
    if(.show.violin){all.y.values <- c(all.y.values, unlist(lapply(densities, function(x){lapply(x, function(y){y$x})})))}
    if(! is.null(.extra.points)){
      all.y.values <- c(all.y.values, unlist(.extra.points))
    }
    .y.limits <- min.max(all.y.values, .na.rm = TRUE)
  } # if(is.null(.y.limits)){
  
  ## Background color
  if(is.null(.background)){
    .background <- .theme
  }
  
  ## Plot area
  if(is.null(.margin)){
    par(mar = c(5, 5.5, 5, 1) * .zoom,
        bg = .background)
  }else{
    par(mar = .margin * .zoom,
        bg = .background)  
  }
  
  ## Set up plot
  plot(x = NULL,
       y = NULL,
       xlim = .x.limits,
       ylim = .y.limits,
       xaxt = 'n',
       yaxt = 'n',
       xlab = '',
       ylab = '',
       main = '', 
       bty = 'n')
  title(.title,
        cex.main = text.size.big * .zoom,
        col.main = ifelse(.theme == 'black', 'white', 'black'))
  title(ylab = .y.label,
        cex.lab = text.size.big * .zoom,
        line = 3.75 * .zoom,
        #col.lab = ifelse(.theme == 'black', 'white', 'black')
        col.lab = .y.axis.color)
  title(xlab = .x.label,
        cex.lab = text.size.big * .zoom,
        line = 3.5 * .zoom,
        col.lab = ifelse(.theme == 'black', 'white', 'black'))
  
  ### Plot any background lines
  if(! is.na(.horizontal.line.at)){
    for(i in 1:length(.horizontal.line.at)){
      abline(h = .horizontal.line.at[i],
             col = ifelse(.theme == 'black', 'white', 'black'),
             lwd = 1.5 * .zoom,
             lty = 2) 
    }
  } # if(! is.na(.horizontal.line.at)){
  
  
  ###
  ### Plot data!
  ###
  
  if(.behind == "dots"){plot.order <- c("dots", "violins")}
  if(.behind %in% c("violin","violins")){plot.order <- c("violins","dots")}
  
  for(plot.type.loop in plot.order){
    # plot.type.loop = plot.order[1]
    
    ### Barplot
    if(.show.bars){
      for(row.loop in 1:nrow(plot.df)){
        # row.loop = 1
        polygon(x = c(plot.df$x[row.loop] - .5 + .bar.spacing,
                      plot.df$x[row.loop] + .5 - .bar.spacing,
                      plot.df$x[row.loop] + .5 - .bar.spacing,
                      plot.df$x[row.loop] - .5 + .bar.spacing),
                y = rep(c(0, plot.df$mean[row.loop]), each = 2),
                border = plot.df$bar.border.color[row.loop],
                lwd = 2 * .zoom,
                col = plot.df$bar.color[row.loop])
      }; rm(row.loop)
    }
    
    ### Plot points
    for(row.loop in 1:nrow(plot.df)){
      # row.loop = 1
      super.label <- plot.df$super.label[row.loop]
      sub.label <- plot.df$sub.label[row.loop]
      
      ### Plot points!
      if(.show.dots & (plot.type.loop == "dots")){
        if(.jitter){
          current.x <- plot.df$x[row.loop] + .jitter.values[[super.label]][[sub.label]]
        }else{
          current.x <- rep(plot.df$x[row.loop], times = length(.values[[super.label]][[sub.label]]))
        }
        points(x = current.x,
               y = .values[[super.label]][[sub.label]],
               pch = 20,
               cex = .point.size * .zoom,
               col = adjust.transparency(.dot.colors[sub.label], 
                                         alpha = ifelse(is.null(.dot.alpha), .2, .dot.alpha)))
      } # if(.show.dots & (plot.type.loop == "dots")){
      
      ### Plot density!
      if(.show.violin & (plot.type.loop == "violins")){
        if(.violin.side == 'both'){
          polygon(x = c(plot.df$x[row.loop] - densities[[super.label]][[sub.label]]$y, 
                        rev(plot.df$x[row.loop] + densities[[super.label]][[sub.label]]$y)),
                  y = c(densities[[super.label]][[sub.label]]$x, 
                        rev(densities[[super.label]][[sub.label]]$x)),
                  border = "NA",
                  col = adjust.transparency(.violin.colors[sub.label], 
                                            alpha = ifelse(is.null(.violin.alpha), .85, .violin.alpha)))
        }else{ # if(.violin.side == 'both'){
          density.x.direction <- ifelse(.violin.side == 'right', 1, -1)
          polygon(x = c(plot.df$x[row.loop] + density.x.direction * ((.85 * .midline.offset) + densities[[super.label]][[sub.label]]$y), 
                        plot.df$x[row.loop] + density.x.direction * ((.85 * .midline.offset) + densities[[super.label]][[sub.label]]$y[1])),
                  y = c(densities[[super.label]][[sub.label]]$x, 
                        densities[[super.label]][[sub.label]]$x[1]),
                  border = "NA",
                  col = adjust.transparency(.violin.colors[sub.label], 
                                            alpha = ifelse(is.null(.violin.alpha), .85, .violin.alpha)))
        } # if(.violin.side == 'both'){}else{
      } # if(.show.violin & (plot.type.loop == "violins")){
    }; rm(row.loop)
  }; rm(plot.type.loop)
  
  for(row.loop in 1:nrow(plot.df)){
    # row.loop = 1
    
    ### Mean and error bars
    if(sum(as.numeric(.error.bar.type %in% c('sd','se')), (! is.null(.error.bar.upper.length))) > 1){
      message("Error! More than one type of error bar specified ('.error.bar.type', and '.error.bar.upper.length'). Choosing 'error.bar.upper.length'")
      .error.bar.type <- ''
    }
    if(.error.bar.type == 'se'){
      current.upper.error <- plot.df$mean[row.loop] + (.multiply.error.bars.by * plot.df$se[row.loop])
      current.lower.error <- plot.df$mean[row.loop] - (.multiply.error.bars.by * plot.df$se[row.loop])
    }
    if(.error.bar.type == 'sd'){
      current.upper.error <- plot.df$mean[row.loop] + (.multiply.error.bars.by * plot.df$sd[row.loop])
      current.lower.error <- plot.df$mean[row.loop] - (.multiply.error.bars.by * plot.df$sd[row.loop])
    }
    if(! is.null(.error.bar.upper.length)){
      if(.multiply.error.bars.by != 1){
        message("Warning!! Specified both '.multiply.error.bars.by' and '.error.bar.upper.length'.  Probably by mistake?  Ignoring '.multiply.error.bars.by'.")
      }
      
      # FINISH THIS!
      message('FINISH THIS CODE -- turn things into lists etc.')
      current.upper.error <- 0
      current.lower.error <- 0
    }
    .current.error.color <- plot.df$error.bar.color[row.loop]
    if(.error.bar.type %in% c('se','sd')){
      segments(x0 = plot.df$x[row.loop],
               y0 = current.lower.error,
               y1 = current.upper.error,
               lwd = 2 * .zoom,
               col = .current.error.color)
      if(.error.bar.caps){
        segments(x0 = plot.df$x[row.loop] - .025,
                 x1 = plot.df$x[row.loop] + .025,
                 y0 = current.upper.error,
                 lwd = 2 * .zoom,
                 col = .current.error.color)
        segments(x0 = plot.df$x[row.loop] - .025,
                 x1 = plot.df$x[row.loop] + .025,
                 y0 = current.lower.error,
                 lwd = 2 * .zoom,
                 col = .current.error.color)
      }
    } # if(.error.bar.type %in% c('se','sd')){
    if(.show.mean){
      points(x = plot.df$x[row.loop],
             y = plot.df$mean[row.loop],
             pch = ifelse(is.null(.mean.pch), 20, .mean.pch),  
             cex = .mean.point.size * 1.2 * .zoom,
             col = plot.df$mean.color[row.loop])
    } # if(.show.mean)
    
    if(! is.null(.extra.points)){
      points(x = plot.df$x[row.loop],
             y = plot.df$extra.points[row.loop],
             pch = .extra.point.pch,  
             cex = .extra.point.size * .zoom,
             col = .extra.point.color)
    } # if(! is.null(.extra.points)){
    
  }; rm(row.loop)
  
  if(.show.x.axis){
    if(is.null(.x.ticks)){
      if(! is.null(.x.tick.labels)){message("Disregarding '.x.tick.labels' because '.x.ticks' not specified.")}
      axis(side = 1,
           at = plot.df$x,
           labels = plot.df$sub.label,
           las = ifelse(.x.tick.labels.vertical, 2, 1),
           tck = -.025 * .zoom, # length of tick
           padj = ifelse(.x.tick.labels.vertical, .3, .6) * .zoom, # distance between tick and label
           lwd = 1.5 * .zoom,
           lwd.ticks = 1.5 * .zoom,
           cex.axis = text.size.med * .zoom,
           col = rgb(0,0,0,0),
           col.axis = ifelse(.theme == 'black', 'white', 'black'))  
      if(length(super.labels) > 1){
        axis(side = 3,
             at = sort(sapply(split(plot.df$x, plot.df$super.label), mean)),
             labels = names(sort(sapply(split(plot.df$x, plot.df$super.label), mean))),
             las = 1,
             tck = -.025 * .zoom, # length of tick
             padj = -.6 * .zoom, # distance between tick and label
             lwd = 1.5 * .zoom,
             lwd.ticks = 1.5 * .zoom,
             cex.axis = text.size.med * .zoom,
             col = rgb(0,0,0,1),
             col.axis = ifelse(.theme == 'black', 'white', 'black'))    
      }
    }else{ # if(is.null(.x.ticks)){
      if(is.null(.x.tick.labels)){
        .x.tick.labels <- .x.ticks
      }
      axis(side = 1,
           at = .x.ticks,
           labels = .x.tick.labels,
           las = ifelse(.x.tick.labels.vertical, 2, 1),
           tck = -.025 * .zoom, # length of tick
           padj = ifelse(.x.tick.labels.vertical, .3, .6) * .zoom, # distance between tick and label
           lwd = 1.5 * .zoom,
           lwd.ticks = 1.5 * .zoom,
           cex.axis = text.size.med * .zoom,
           col = rgb(0,0,0,0),
           col.axis = ifelse(.theme == 'black', 'white', 'black'))
    } # if(is.null(.x.ticks)){}else{
  } # if(.show.x.axis){
  if(.show.y.axis){
    if(is.null(.y.ticks)){
      if(! is.null(.y.tick.labels)){message("Disregarding '.y.tick.labels' because '.y.ticks' not specified.")}
      axis(side = 2,
           las = 0,
           tck = -.025 * .zoom, # length of tick
           padj = -.45 * .zoom, # distance between tick and label
           lwd = 1.5 * .zoom,
           lwd.ticks = 1.5 * .zoom,
           cex.axis = text.size.med * .zoom,
           col = .y.axis.color,
           col.axis = .y.axis.color)  
    }else{ # if(is.null(.y.ticks)){
      if(is.null(.y.tick.labels)){
        .y.tick.labels <- .y.ticks
      }
      axis(side = 2,
           at = .y.ticks,
           labels = .y.tick.labels,
           las = 0,
           tck = -.025 * .zoom, # length of tick
           padj = -.45 * .zoom, # distance between tick and label
           lwd = 1.5 * .zoom,
           lwd.ticks = 1.5 * .zoom,
           cex.axis = text.size.med * .zoom,
           col = .y.axis.color,
           col.axis = .y.axis.color)
    } # if(is.null(.y.ticks)){}else{
  } # if(.show.y.axis){
}