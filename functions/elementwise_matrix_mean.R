### Average across dataframes/matrices by element 
### Spring 2022
### adam.milton.morgan@gmail.com

elementwise.matrix.mean <- function(list.of.matrices){
  
  input.class <- class(list.of.matrices[[1]])[1]
  
  output <- Reduce("+", list.of.matrices) / length(list.of.matrices) 
  
  if(input.class == "data.frame"){output <- data.frame(output)}
  
  # End
  return(output)
}