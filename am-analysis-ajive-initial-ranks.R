library(dplyr)
library(magrittr)
library("ggplot2")
library(pajive)

paths <- list(
  project = "~/projects/linalg-exam/"
)

setwd(paths$project)

source("am-helper-functions.R")

filtered_data <- readRDS("data-filtered-mv5-CtM-Cm40-Cx3-Imp.RDS")

filtered_data$blocks %>% data.info

# Manually selected from viewing PL plots

ranklist <- list(
      low_init_ranks = c(miRNA= 3,mRNA= 6, DNAm= 3),
     lowB_init_ranks = c(miRNA= 7,mRNA= 9, DNAm= 8),
     high_init_ranks = c(miRNA=11,mRNA=13, DNAm=12),
    highB_init_ranks = c(miRNA=17,mRNA=20, DNAm=18),
very_high_init_ranks = c(miRNA=22,mRNA=26, DNAm=24)
)

future::plan(future::multisession(), workers=16)
future::plan(future::sequential())

blocklist = c("miRNA", "mRNA", "DNAm")

for(rank_name in names(ranklist)){
  message("ranks", rank_name)
  ajive(blocks = filtered_data$blocks[blocklist],
      initial_signal_ranks = ranklist[[rank_name]][blocklist],
      n_wedin_samples = 300) %>%
    saveRDS(paste0("ajr-", rank_name, "-300.RDS"))
}

data_files <- list(
  "ajr-low_init_ranks-300.RDS",
  "ajr-lowB_init_ranks-300.RDS",
  "ajr-high_init_ranks-300.RDS",
  "ajr-highB_init_ranks-300.RDS",
  "ajr-very_high_init_ranks-300.RDS"
)

tags <- c("low","medium", "high", "higher", "veryhigh")

library(pROC)


d <- readRDS("data-filtered-mv5-CtM-Cm40-Cx3-Imp.RDS")
for(i in seq_along(data_files)){
  ajr <- readRDS(data_files[[i]])
  names(ajr$block_decomps) <- c("miRNA", "mRNA", "DNAm")
  am_rocplot(ajr, d, paste0(tags[i],"-300-7c"))
}

d <- readRDS("data-filtered-mv5-CtM-Cm40-Cx5-Imp.RDS")
ajr <- readRDS("data-ajr-mRNA-mv5-CtM-Cm40-Cx5-Imp.RDS")
  ajr <- readRDS(data_files[[i]])
  names(ajr$block_decomps) <- c("miRNA", "mRNA", "DNAm")
  am_idv_rocplot(ajr, d, "idv-Cx5-nocovs")
  
