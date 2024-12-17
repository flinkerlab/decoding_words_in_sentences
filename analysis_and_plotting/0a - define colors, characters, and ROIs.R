### Plot electrode high gamma and encoding importance over time to complement brain plots
### Feb 2023
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
library('viridis') # for viridis color schemes
library('pals') # for ocean.balance, ocean.delta, ocean.curl colors
#pal.bands(ocean.balance, ocean.delta, ocean.curl)
library('scico') # for palettes cork, vik, lisbon, tofino, berlin, roma
library('reticulate') # for Python
library('stringr') # for str_split()
library('lme4') # for lmer()
library('lmerTest') # for lmer() that gives p-value
library('NMF') # for non-negative matrix factorization
library('colorspace')
library('shades')

# Clean up
rm(list=ls())
cat("\014")

### Set path
if(Sys.info()['sysname'] == 'Darwin'){ # Mac
  path = '/Users/adam/Dropbox/Research/ChickenSyntax/'
  if(Sys.info()['nodename'] == 'FLINKERLABMBP06'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
  }
  if(Sys.info()['nodename'] == 'FLINKERLABMS01'){
    path = '/Users/am4611/Dropbox/Research/ChickenSyntax/'
  }
}
if(Sys.info()['sysname'] == 'Linux'){ # Ubuntu
  path = '/home/adam/Dropbox/Research/ChickenSyntax/'
}

### Set up elec info
elec.info <- 
  read.csv(paste0(path,
                  'analysis/R/brain plots/ecog/output/data/elec info/patients - combined/row_labels_and_localizations.csv'))
elec.info$use.these.elecs <- 
  as.numeric((! elec.info$region_clinical %in% c('Unknown',
                                                 '',
                                                 'NaN',
                                                 'Left-Cerebral-White-Matter',
                                                 'Left-Inf-Lat-Vent',
                                                 'Right-Cerebral-White-Matter')) &
               (elec.info$bad_elec == 0) &
               (elec.info$visual_elec == 0))
rownames(elec.info) <- elec.info$patient_elec
elec.info$region_clinical_for_plots <- gsub("ctx-lh-", "", elec.info$region_clinical, fixed = TRUE)
use.these.elecs <- elec.info[elec.info$use.these.elecs == 1,]$patient_elec


for(band.loop in c('high_gamma','beta')){
  # band.loop = c('high_gamma','beta')[1]
  
  output.path <- paste0(path, 'analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/',band.loop,'/')
  
  
  ###
  ### Elecs
  ###
  
  # Save
  elec.info.path <- paste0(output.path, 'data/0a - definitions/elec info/')
  dir.create(elec.info.path, showWarnings = FALSE, recursive = TRUE)
  save(list = c('elec.info','use.these.elecs'), 
       file = paste0(elec.info.path, 'elec info.RData'))
  write.csv(elec.info,
            file = paste0(elec.info.path, 'elec info.csv'),
            row.names = FALSE, quote = FALSE)
  
  
  ###
  ### Words
  ###
  
  words <- c('chicken','dog','dracula','frankenstein','ninja','nurse')
  
  # Save
  words.path <- paste0(output.path, 'data/0a - definitions/characters/')
  dir.create(words.path, showWarnings = FALSE, recursive = TRUE)
  save(words, 
       file = paste0(words.path, 'words.RData'))
  write.csv(data.frame('words' = words),
            file = paste0(words.path, 'words.csv'),
            row.names = FALSE, quote = FALSE)
  
  
  ###
  ### Define colors for characters
  ###
  
  character.colors.og <- 
    character.colors <- 
    data.frame('word' = words,
               'hex' = setNames(palette.colors(12, palette = "Set3")[c(12, 6, 5, 7, 4, 1)], words))
  # Make frankenstein darker
  character.colors[character.colors$word == 'frankenstein','hex'] <- rgb(35, 204, 75, maxColorValue = 255)
  character.colors$hex <- as.character(shades::saturation(character.colors$hex, values = delta(.8)))
  character.colors <- cbind(character.colors, t(col2rgb(character.colors$hex)))
  
  # Save
  character.colors.path <- paste0(output.path, 'data/0a - definitions/character colors/')
  dir.create(character.colors.path, showWarnings = FALSE, recursive = TRUE)
  save(list = c("character.colors",
                "character.colors.og"), 
       file = paste0(character.colors.path, 'character colors.RData'))
  write.csv(character.colors,
            file = paste0(character.colors.path, 'character_colors.csv'),
            row.names = FALSE, quote = FALSE)
  
  
  ###
  ### Define colors for sentences and lists
  ###
  
    ## Active/passive noun colors
  noun.colors <- list()
  noun.colors[['passive']] <- 
    noun.colors[['active']] <- 
    data.frame('noun' = c('noun1','noun2'),
               'hex' = character.colors[c('dracula','frankenstein'),'hex'])
  noun.colors[['passive']]$hex <- 
    colorspace::darken(as.character(shades::saturation(noun.colors[['passive']]$hex, values = delta(-.2))), 
                       amount = .4, method = "relative")
  noun.colors[['list']] <- noun.colors[['sentence']] <- noun.colors[['active']]
  # Add RGB
  for(syntax.loop in names(noun.colors)){
    noun.colors[[syntax.loop]] <- cbind(noun.colors[[syntax.loop]], t(col2rgb(noun.colors[[syntax.loop]]$hex)))
    rownames(noun.colors[[syntax.loop]]) <- noun.colors[[syntax.loop]]$noun
  }; rm(syntax.loop)
  
  
  ## Save
  sentence.color.path <- paste0(output.path, 'data/0a - definitions/noun colors/')
  dir.create(sentence.color.path, showWarnings = FALSE, recursive = TRUE)
  save(noun.colors, 
       file = paste0(sentence.color.path, 'noun colors.RData'))
  write.csv(noun.colors[['active']],
            file = paste0(sentence.color.path, 'active_noun_colors.csv'),
            row.names = FALSE, quote = FALSE)
  write.csv(noun.colors[['passive']],
            file = paste0(sentence.color.path, 'passive_noun_colors.csv'),
            row.names = FALSE, quote = FALSE)
  
  
  ###
  ### Define regions and colors
  ###
  
  ### Define ROIs
  roi.categories <- list(
    'IFG' = c('parsopercularis','parstriangularis','ctx-lh-parsopercularis', 'parsorbitalis'),
    'SMC' = c('precentral','postcentral'),
    'ATL' = c('temporalpole','rMTG','rSTG'),
    'MTL' = c('mMTG','mSTG'),
    'PTL' = c('cMTG','cSTG'),
    'IPL' = c('supramarginal','inferiorparietal'),
    'MFG' = c('rostralmiddlefrontal','caudalmiddlefrontal'),
    'STG' = c('cSTG','mSTG','rSTG','ctx-lh-superiortemporal'),
    'MTG' = c('cMTG','mMTG','rMTG','ctx-lh-middletemporal')
    # 'TPL' = c('temporalpole') # insanely high but only 2 elecs
  )
  rois.all <- names(roi.categories)
  
  # ROI elecs
  roi.elecs <- list()
  for(roi.loop in rois.all){ # roi.loop = rois[1]
    roi.elecs[[roi.loop]] <- elec.info[elec.info$region_clinical %in% roi.categories[[roi.loop]], 'patient_elec']
  }; rm(roi.loop)
  
  # Temporal lobe splits: either along the rostral-caudal axis or dorsal-ventral:
  temporal.lobe.splits <- list(
    'rostral.caudal' = c('IFG','MFG','SMC','ATL','MTL','PTL','IPL'),
    'dorsal.ventral' = c('IFG','MFG','SMC','STG','MTG','IPL'))
  roi.colors <- lapply(temporal.lobe.splits, function(x){
    setNames(cubicl(length(x)), x)})
  roi.colors.dfs <- lapply(roi.colors, function(x){
    x <- data.frame('roi' = names(x),
                    'hex' = x)
    x <- cbind(x, t(col2rgb(x$hex)))
  })
  
  # Save
  save.roi.stuff <- c('rois.all','roi.categories','roi.elecs','temporal.lobe.splits','roi.colors','roi.colors.dfs')
  save.roi.dir <- paste0(output.path, 'data/0a - definitions/roi metadata/')
  dir.create(save.roi.dir, showWarnings = FALSE, recursive = TRUE)
  save(list = save.roi.stuff,
       file = paste0(save.roi.dir, 'ROIs and temporal lobe splits.RData'))
  
  for(tl.split.loop in names(roi.colors.dfs)){
    write.csv(roi.colors.dfs[[tl.split.loop]],
              paste0(save.roi.dir, 'roi colors with temporal lobe split along ',tl.split.loop,' axis.csv'),
              row.names = FALSE, quote = FALSE)
  }; rm(tl.split.loop)
  
} # band.loop





message('Script completed successfully. ',Sys.time())
