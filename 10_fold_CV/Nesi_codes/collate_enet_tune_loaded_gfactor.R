## library path
.libPaths(new = "/scale_wlg_persistent/filesets/home/wanyu368/R/foss-2023a/4.3" )
#for parallel processing

library(tidyverse)

#dataset <- "DTI_data.RDS"
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
data_name <- gsub(dataset, pattern = "\\.RDS", replacement = "" )

datafolder <- "/nesi/nobackup/uoo03493/Yue/stacking_gfactor_10foldcv/"
#datafolder <- "/Volumes/sci-psy-narun/Nesi/Yue/stacking_gfactor_10foldcv/"
data_input <- readRDS(file = paste0(datafolder, "/data/", dataset))

data_input <- data_input %>% drop_na(SITE_ID_L)
fold_info <- readRDS(file = paste0(datafolder, "data/part_fold.RDS"))


batch_ids <- unique(fold_info$fold)

results <- list()

for(batch_idx in 1:length(batch_ids)){
tmp <- readRDS(file = paste0(datafolder, "enet_results_loaded_gfactor/",data_name, ".",batch_idx,".RDS" ))
  results[[batch_ids[batch_idx]]] <- tmp
}


saveRDS(results, file = paste0(datafolder, "collect_enet_results_loaded_gfactor/",data_name,".enet_results.RDS"))
