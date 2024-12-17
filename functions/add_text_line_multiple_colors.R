### Plot time series - general function
### Winter 2023 (!!!!!!)
### adam.milton.morgan@gmail.com

add.text.line.multiple.colors <- 
  function(text.segments,
           text.colors,
           .side = 3,
           .line = 0,
           .outer = FALSE,
           .cex = 1){
    
    # Check inputs
    if(length(text.segments) != length(text.colors)){
      message("Number of text segments and colors don't match!")
    }
    
    # Loop thru text segments
    for(segment.loop in 1:length(text.segments)){
      # segment.loop = 7
      
      # Get the text before and after the current segment
      if(segment.loop == 1){
        previous.segments <- c()
      }else{
        previous.segments <- 1:(segment.loop - 1)
      }
      if(segment.loop == length(text.segments)){
        subsequent.segments <- c()
      }else{
        subsequent.segments <- (segment.loop + 1):length(text.segments)
      }
      
      # Define text before and after current text segment
      previous.text <- paste(paste0(text.segments[previous.segments],' '), collapse = '')
      subsequent.text <- paste(paste0(text.segments[subsequent.segments],' '), collapse = '')
      current.text <- text.segments[segment.loop]
      
      # Get rid of leading and trailing spaces
      if(previous.text == ' '){previous.text <- ''}
      if(subsequent.text == ' '){subsequent.text <- ''}
      
      # Write just the current word in its color
      mtext(bquote(phantom(.(previous.text)) * 
                     .(current.text) * 
                     phantom(.(subsequent.text))), 
            side = .side,
            col = text.colors[segment.loop],
            cex = .cex,
            outer = .outer,
            line = .line)
    }; rm(segment.loop)
  }