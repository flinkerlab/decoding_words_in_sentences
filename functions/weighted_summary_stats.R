### Weighted mean and standard error
### June 2024
### adam.milton.morgan@gmail.com

## Function for weighted average
weighted.average <- function(x, w){
  ## Sum of the weights 
  sum.w <- sum(w, na.rm = TRUE)
  ## Sum of the weighted $x_i$ 
  xw <- sum(w*x, na.rm = TRUE)
  
  ## Return the weighted average 
  return(xw/sum.w)
}

## Function for getting weighted squared error (thx to Alex Stephenson)
weighted.se.mean <- function(x, w, na.rm = TRUE){
  if(na.rm){
    keep.these <- which(! is.na(x))
    x <- x[keep.these]
    w <- w[keep.these]
  }
  
  ## Calculate effective N and correction factor
  n_eff <- (sum(w))^2/(sum(w^2))
  correction = n_eff/(n_eff-1)
  
  ## Get weighted variance 
  numerator = sum(w*(x-weighted.average(x,w))^2)
  denominator = sum(w)
  
  ## get weighted standard error of the mean 
  se_x = sqrt((correction * (numerator/denominator))/n_eff)
  return(se_x)
}