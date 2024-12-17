### Adjust transparency of color
### Summer 2021
### adam.milton.morgan@gmail.com
### Based on a script by Mark Gardener 2015

adjust.transparency <- function(.colors, 
                                alpha = 120,
                                max.alpha = NA){
  
  # Get RGB values for named color
  rgb.vals <- col2rgb(.colors)
  
  # Convert alpha to the 1 to 255 scale
  if(is.na(max.alpha)){
    if(alpha <= 1){
      max.alpha <- 1
    }else{
      max.alpha <- 255
    }
  }
  if(max.alpha != 255){
    alpha <- 255 * alpha / max.alpha
  }
  
  # Repeat alpha if only one val
  if(length(alpha) == 1){
    alpha <- rep(alpha, times = length(.colors))
  }
  
  # Make new color using input color as base and alpha set by transparency
  output.colors <- c()
  for(i in 1:length(.colors)){ # i = 2
    output.colors <- c(output.colors,
                       rgb(rgb.vals[1,i], rgb.vals[2,i], rgb.vals[3,i],
                           max = 255,
                           alpha = alpha[i]))
  }
  
  # Return the color
  return(output.colors)
}

