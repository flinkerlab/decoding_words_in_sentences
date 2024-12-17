### Functions for subsetting data to just good trials
### Spring 2023
### adam.milton.morgan@gmail.com

###
### Picture naming
###

just.good.pic.naming <- function(.data,
                                 .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                 .remove.outliers = TRUE,
                                 .remove.over.3000ms = TRUE,
                                 .rt.quantile.range = c(.025,.95),
                                 .rt.range.ms = NULL){
  
  # .data: a dataframe with at least columns for word, substitution_error, pos, production_latency, and case
  # .words: the words to subset to
  # .remove.outliers: exclude trials with RTs greater than or less than specified values (either .rt.quantile.range or .rt.range.ms)
  # .rt.quantile.range: the RT quantiles above and below which to exclude trials; provide either this or the range in RTs (.rt.range.ms; milliseconds)
  # .rt.range.ms: the RTs above and below which to exclude trials; provide either this or the range in quantiles (.rt.quantile.range; proportions)
  
  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good picture naming trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(.data$case == 'none') & 
                              is.na(.data$substitution_error) &
                              (.data$pos == 'n') &
                              (.data$word %in% .words) &
                              (! is.na(.data$case)),])
  
  if(.remove.over.3000ms){
    .data <- droplevels(.data[.data$production_latency < time.convert(3000, 'times', 'samples'),])
  }
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
}


###
### Sentence data -- first word only
###

just.first.word.sentence.data <- function(.data,
                                          .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                          .remove.over.3000ms = TRUE,
                                          .remove.outliers = FALSE,
                                          .rt.quantile.range = c(.025,.95),
                                          .rt.range.ms = NULL){
  
  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good sentence production trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(.data$case == 's') & 
                              is.na(.data$substitution_error) &
                              (((.data$dp_or_np == 'np') & (.data$pos == 'n') & (.data$word %in% .words)) | 
                                 ((.data$dp_or_np == 'dp') & (.data$pos == 'd') & (.data$word == 'the'))) &
                              (! is.na(.data$case)),])
  
  if(.remove.over.3000ms){
    .data <- droplevels(.data[.data$production_latency < time.convert(3000, 'times', 'samples'),])
  }
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
} # just.first.word.sentence.data()


###
### List data -- first word only
###

just.first.word.list.data <- function(.data,
                                      .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                      .remove.over.3000ms = TRUE,
                                      .remove.outliers = FALSE,
                                      .rt.quantile.range = c(.025,.95),
                                      .rt.range.ms = NULL){
  
  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good list production trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(.data$case == '1') & 
                              is.na(.data$substitution_error) &
                              (.data$pos == 'n') &
                              (.data$word %in% .words) &
                              (! is.na(.data$case)),])
  
  if(.remove.over.3000ms){
    .data <- droplevels(.data[.data$production_latency < time.convert(3000, 'times', 'samples'),])
  }
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
} # just.first.word.list.data()


###
### Sentence data -- subjects, verbs, and objects
###

just.sentence.noun.verb.data <- function(.data,
                                         .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                         .remove.outliers = FALSE,
                                         .rt.quantile.range = c(.025,.95),
                                         .rt.range.ms = NULL){
  
  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good sentence production trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(
    # Current word is either a subject, direct object, or by-object
    ((
      # More specifically -- either an active subject/direct-object...
      ((.data$case %in% c("do","s")) & (.data$verb_voice == "active")) | 
        # ...or a passive subject/by-object
         ((.data$case %in% c("bo","s")) & (.data$verb_voice == "passive"))
      ) &
        # in which case it's one of the 6 words
        (.data$pos == 'n') &
        (.data$word %in% .words)) |
       # or a it's verb
       (.data$pos == 'v')) &
      # and it's not an error
      is.na(.data$substitution_error),])
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
} # just.sentence.noun.verb.data()


###
### List data -- just nouns
###

just.list.noun.data <- function(.data,
                                .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                .remove.outliers = FALSE,
                                .rt.quantile.range = c(.025,.95),
                                .rt.range.ms = NULL){
  
  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good list production trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(.data$case %in% c("1","2")) &
                              (.data$pos == 'n') &
                              (.data$word %in% .words) &
                              (! is.na(.data$case)) &
                              is.na(.data$substitution_error),])
  
  if(.remove.outliers){
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
} # just.first.word.list.data()



###
### Sentence data -- dets, subjects, aux's, passive be's, verbs, by's, objects
###

just.good.sentence.data <- function(.data,
                                         .words = c('chicken','dog','dracula','frankenstein','ninja','nurse'),
                                         .remove.outliers = FALSE,
                                         .rt.quantile.range = c(.025,.95),
                                         .rt.range.ms = NULL){
  
  source(paste0(path, 'analysis/R/functions/time_convert.R'))
  
  # Subset to just good sentence production trials, removing RT outliers (>3000ms)
  .data <- droplevels(.data[(
    # Current word is either a subject, direct object, or by-object
    ((
      # More specifically -- either an active subject/direct-object...
      ((.data$case %in% c("do","s")) & (.data$verb_voice == "active")) | 
        # ...or a passive subject/by-object
        ((.data$case %in% c("bo","s")) & (.data$verb_voice == "passive"))
    ) &
      # in which case it's 'the' or one of the 6 words
      (.data$pos %in% c('n','d')) &
      (.data$word %in% c('the', .words))) |
      # or a it's verb/aux/passive be/passive by
      (.data$pos %in% c('v','a','b','p'))) &
      # and it's not an error
      is.na(.data$substitution_error),])
  
  if(.remove.outliers){
    message("Warning (just.good.sentence.data()): Removing outliers in RT across all sentence roles -- will probably just remove passive by-objects since RT is measured from stimulus onset. Should probably edit funciton or set to FALSE.")
    # Get all RTs (in samples)
    rts.samples <- .data$production_latency
    
    if(is.null(.rt.range.ms)){
      # Get cutoff RTs (in samples)
      rts.samples.lower.bound <- quantile(rts.samples, .rt.quantile.range[1])
      rts.samples.upper.bound <- quantile(rts.samples, .rt.quantile.range[2])
    }else{
      rts.samples.lower.bound <- time.convert(.rt.range.ms[1], "times", "samples")
      rts.samples.upper.bound <- time.convert(.rt.range.ms[2], "times", "samples")
    }
    
    # Subset to just trials within cutoff RTs
    .data <- droplevels(.data[(.data$production_latency >= rts.samples.lower.bound) &
                                (.data$production_latency <= rts.samples.upper.bound),])
  } # if remove outliers
  
  # Return trimmed data
  return(.data)
} # just.good.sentence.data()
