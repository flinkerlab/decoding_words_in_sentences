### Adjust a pheatmap object so when you save it as a PDF there aren't white borders everywhere
### Summer 2023
### adam.milton.morgan@gmail.com
### Based on a script by Mark Gardener 2015

touchup.heatmap <- function(.heatmap, 
                            .theme = 'white',
                            .zoom = 1){
  
  ### Change the cell borders from white to the cell color
  # Get the names of the grobs corresponding to the individual cells (rectangles)
  grob.classes <- purrr::map(.heatmap$gtable$grobs, class)
  cell.indices <- which(purrr::map_lgl(grob.classes, function(cl) 'gTree' %in% cl))
  
  for(index in cell.indices){
    # Find the rectangle grobs
    grob_names <- names(.heatmap$gtable$grobs[[index]]$children)
    idx_rect <- grob_names[grep('rect', grob_names)][1]
    
    # Replace the color parameter for the border with the color parameter for the fill
    .heatmap$gtable$grobs[[index]]$children[[idx_rect]]$gp$col <-
      .heatmap$gtable$grobs[[index]]$children[[idx_rect]]$gp$fill
  }; rm(index)
  
  # And the color labels grob:
  rectangle.indices <- which(purrr::map_lgl(grob.classes, function(cl) 'rect' %in% cl))
  for(index in rectangle.indices){
    .heatmap$gtable$grobs[[index]]$gp$col <- 
      .heatmap$gtable$grobs[[index]]$gp$fill
  }; rm(index)
  
  
  ### Change background and font colors
  # Set colors
  text.color <- ifelse(.theme == 'white', 'black', 'white')
  
  # Get text indices
  text.indices <- which(purrr::map_lgl(grob.classes, function(cl) 'text' %in% cl))
  for(index in 1:length(grob.classes)){
    # index = 7
    
    # If this is a text grob, change color
    if(index %in% text.indices){
      .heatmap$gtable$grobs[[index]]$gp$col <- text.color
    }else{ # if(index %in% text.indices){
      # Otherwise look for text grobs
      if(! is.null(.heatmap$gtable$grobs[[index]]$children)){
        text.subindices <- which(purrr::map_lgl(purrr::map(.heatmap$gtable$grobs[[index]]$children, class), function(cl) 'text' %in% cl))
        for(index2 in text.subindices){
          .heatmap$gtable$grobs[[index]]$children[[index2]]$gp$col <- text.color
        }; rm(index2)
      } # if(! is.null(.heatmap$gtable$grobs[[index]]$children)){
    } # if(index %in% text.indices){}else{
  }; rm(index)
  
  # And the colorbar grob
  if(length(.heatmap$gtable$grobs) >= 5){
    colorbar.text.indices <- grep("text", names(.heatmap$gtable$grobs[[5]]$children))
    for(index in colorbar.text.indices){
      .heatmap$gtable$grobs[[5]]$children[[index]]$gp$col <- text.color
    }; rm(index)
  } # if(length(.heatmap$gtable$grobs) >= 5){
  
  
  ### Change dendrogram line colors and widths
  .heatmap$gtable$grobs[[1]]$gp <- gpar(lwd = 1 * .zoom, 
                                        col = text.color,
                                        lineend = 'butt',
                                        linejoin = 'mitre')
  .heatmap$gtable$grobs[[2]]$gp <- gpar(lwd = 1 * .zoom, 
                                        col = text.color,
                                        lineend = 'butt',
                                        linejoin = 'mitre')
  
  
  ### Return
  return(.heatmap)
} # touchup.heatmap()
