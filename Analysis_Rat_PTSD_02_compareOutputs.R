### Open .chp files from TAC
library(affxparser)
rmaSMA_norm = list.files('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/TAC/Summarization_RMA/apt_output/cc-chp',full.names=TRUE,pattern =".chp")
lapply( 
  dat  = affxparser::readChp(rmaSMA_norm[1], withQuant = TRUE)
