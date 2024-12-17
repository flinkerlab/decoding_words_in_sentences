### Convert from logits to proportions and back
### Summer 2024
### adam.milton.morgan@gmail.com

### Logit transform (for doing stats on proportions)
proportion.to.logit <- function(proportions, 
                                define.min.max.from.distribution = TRUE,
                                .min = .01, 
                                .max = .99){
  # Can't logit proportions of 0 or 1, so set to reasonable low/high proportions
  if(all(proportions %in% c(0, 1))){
    # message('WARNING (proportion.to.logit()): Weird... all proportions are 0 or 1. That right?')
    # For cases like the patient who only had one passive object trial where, so all results were either 0 or 1
    proportions[proportions == 0] <- .min
    proportions[proportions == 1] <- .max
  }else{
    if(define.min.max.from.distribution & (length(proportions) > 10)){
      .min <- min(proportions[proportions != 0]) / 2
      .temp.max <- (1 - max(proportions[proportions != 1])) / 2
      # Make symmetrical
      .min <- .temp.max <- min(c(.min, .temp.max))
      .max <- 1 - .temp.max
    }
  }
  proportions[proportions == 0] <- .min
  proportions[proportions == 1] <- .max
  logits <- log(proportions / (1 - proportions))
  return(logits)
} # proportion.to.logit()


### Inverse logit 
logit.to.proportion <- function(logits){
  # Just to be safe, curb very high and very low values (I don't think there shouldn't be any values like this)
  logits[logits > 10] <- 10
  logits[logits < -10] <- -10
  # Calculate!
  proportions <- exp(logits) / (1 + exp(logits))
  return(proportions)
} # logit.to.proportion()