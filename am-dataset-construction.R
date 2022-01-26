library(dplyr)
library(magrittr)
library("ggplot2")
library(pajive)
library(stringr)
library(purrr)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")

#=#=#=#=#=#=#=#=#=#=#=#
#### CONFIGURATION ####
#=#=#=#=#=#=#=#=#=#=#=#

# ~ Paths ####
paths <- list(
  project = "~/projects/linalg-exam/",
  section = list(loading = list( save_file = "data-full-complete.RDS" ),
                   setup = list( save_file = "data-setup.RDS"),
                  #filter = list( save_file = "data-filtered-STD.RDS") ),
                  filter = list( save_file = function(tags){paste0("data-filtered-", paste0(tags,collapse="-"),".RDS")}) ),
   blocks = list(
    miRNA = "miRNA.LC.RDS",
     DNAm = "DNAm.beta.LC.RDS",
     mRNA = "mRNA-exprs.LC.RDS",
     covs = "covs.LC.RDS"
   )  
)

# ~ Set working directory ####

setwd(paths$project)

# ~ Load helper function
source("am-helper-functions.R")

#=#=#=#=#=#=#=#=#=#=#
#####  LOADING  #####
#=#=#=#=#=#=#=#=#=#=#

# ~ Load raw data ####
raw_data <- list(
  blocks = map(paths$blocks, readRDS),
  log = character()
)

raw_data$blocks %>% data.info

# ~ Correct orientation ####
blocks_to_transpose = c("DNAm", "mRNA")
raw_data$blocks %<>% map_at(blocks_to_transpose, t)

raw_data$mods %<>% c("transpose")

# ~ remove casepairs column  ####
raw_data$blocks$miRNA %<>% select(-casepairs) 
raw_data$mods %<>% c("miRNA-casepairs")

# ~ data.frames -> matrix ####
raw_data$blocks$miRNA %<>% as.matrix
raw_data$mods %<>% c("matrix(miRNA)")
raw_data$blocks$covs  %<>% as.matrix
raw_data$mods %<>% c("matrix(covs)")

raw_data$blocks %>% data.info

# ~ Identify case/ctrl pairs ####
complete_pairs <-  raw_data$blocks %>%
  map(rownames) %>% reduce(intersect) %>% 
  {  num <- function(x) str_remove(x, "case |ctrl ")
     subset(.,  num(.) %in% num(.)[duplicated(num(.))] )  }

length(complete_pairs)

# ~ Select ordered pairs ####
raw_data$blocks %<>% map(~.[complete_pairs,]) 
raw_data$mods %<>% c("complete pairs", "same row order")

raw_data$blocks %>% data.info
  
# ~ [save loaded data] ####
saveRDS(raw_data, paths$section$loading$save_file)


#=#=#=#=#=#=#=#=#=#=#
##### SETUP #####
#=#=#=#=#=#=#=#=#=#=#

# ~ [load saved data] ####
full_data <- readRDS(paths$section$loading$save_file)

# ~ miRNA - log2 conversion  #####
full_data$blocks$miRNA <- log2(full_data$blocks$miRNA + 1) 
full_data$mods %<>% c("log2(miRNA)")

# ~ DNAm - M-value conversion  #####
full_data$blocks$DNAm <- log2((full_data$blocks$DNAm)/(1-full_data$blocks$DNAm)) # TODO: Does this need a constant added?
full_data$mods %<>% c("m-value(DNAm)")

# ~ [save loaded data] ####
saveRDS(full_data, paths$section$setup$save_file)

#=#=#=#=#=#=#=#=#=#=#
##### FILTERING #####
#=#=#=#=#=#=#=#=#=#=#

# ~ [load saved data] ####
filtered_data <- readRDS(paths$section$setup$save_file)

filtered_data$blocks %>% data.info

# ~ mRNA - High Variance  #####

# Identify highest variance mRNAs
source('functions/filtering.R')
filtered_data$blocks$mRNA %<>% filter_top_cols(var, keep = 5000) 
filtered_data$mods %<>% c("top5000varmRNA")
filtered_data$nametags %<>% c("mv5")

filtered_data$blocks %>% data.info

# ~  DNAm - mRNA intersection  #####
# pick CpGs that are on the same genes as mRNA 

# load CpG-gene associations
#source("build-cpg-gene-name-dataset.R")
cpg450k_annot <- read.csv("cpg-450k-tidy-annot.csv", stringsAsFactors=FALSE)

# load mRNA-gene associations
#source("build-mrna-gene-mapping-dataset.R")
mRNAmapping <- readRDS("mRNA-gene-mapping.RDS")

# TODO: encode additional columns indicating cpg region type
# Build CpG <-> mRNA mapping
DNAm.mRNA <- cpg450k_annot[c("cgName", "geneName")] %>%  unique %>%
  left_join( mRNAmapping, by="geneName") %>%
  filter( mRNAname %in% colnames(filtered_data$blocks$mRNA) ) %>% # limit to present mRNAs
  dplyr::select(-mRNAname) %>% 
  unique

# Build CpG <-> mRNA mapping : old incorrect way - geneName is not split on ';'s.
# DNAm.mRNA.orig <- annot450k[c("cgName", "geneName")] %>%  unique %>%
#   left_join( mRNAmapping, by="geneName") %>%
#   filter( mRNAname %in% colnames(filtered_data$blocks$mRNA) ) %>% # limit to present mRNAs
#   dplyr::select(-mRNAname) %>% 
#   unique
  
# only keep mRNA-associated CpGs in DNAm
filtered_data$blocks$DNAm %<>% extract(, colnames(filtered_data$blocks$DNAm) %in% DNAm.mRNA$cgName)

filtered_data$blocks %>% data.info
filtered_data$mods %<>% c("CpG<->mRNA")
filtered_data$nametags %<>% c("CtM")

rm(cpg450k_annot, mRNAmapping, DNAm.mRNA)

# ~ DNAm - High Variance  #####

# Identify highest variance CpG sites
source('functions/filtering.R')
filtered_data$blocks$DNAm %<>% filter_top_cols(var, keep = 98000)
filtered_data$mods %<>% c("top98000varDNAm")
filtered_data$nametags %<>% c("CtV98")

filtered_data$blocks %>% data.info



# ~ DNAm - excess missing  #####

# eliminate CpGs with more than 40% missing
missingCpG_threshhold = 0.4
naprop <- apply(is.na(filtered_data$blocks$DNAm), 2, mean)

filtered_data$blocks$DNAm %<>% extract(,naprop <= missingCpG_threshhold) # TODO: Sensitivity testing

filtered_data$mods %<>% c(sprintf("LT%0.2fMissDNAm", missingCpG_threshhold))
filtered_data$nametags %<>% c("Cm40")

rm(missingCpG_threshhold, naprop)
filtered_data$blocks %>% data.info

# ~ DNAm - M-value filtering  #####

# eliminate extreme M values - absolute threshhold or percentile threshhold
  absThreshold =  5
upperThreshold =  absThreshold
lowerThreshold = -absThreshold
meanM <- apply(filtered_data$blocks$DNAm, 2, function(x) mean(x, na.rm =TRUE)) # TODO: Sensitivity testing, should this be 5?
filtered_data$blocks$DNAm %<>% extract(,between(abs(meanM),lowerThreshold, upperThreshold))

filtered_data$mods %<>% c(paste0("btw(",paste0(c(upperThreshold, lowerThreshold), collapse=","), ")DNAm"))
filtered_data$nametags %<>% c("Cx5")
filtered_data$blocks %>% data.info
rm(meanM, absThreshold, upperThreshold, lowerThreshold)

# ~ ALL - matrix completion #####


completion_list <- c("miRNA", "DNAm", "mRNA")
system.time(
  filtered_data$blocks[completion_list] %<>% map(~SVDmiss(.)$Xfill)
) # 71 seconds
filtered_data$blocks$mRNA  %<>% `rownames<-`(rownames(filtered_data$blocks$covs))
filtered_data$blocks$miRNA %<>% `rownames<-`(rownames(filtered_data$blocks$covs))
filtered_data$blocks$DNAm  %<>% `rownames<-`(rownames(filtered_data$blocks$covs))

filtered_data$blocks %>% data.info

rm(ajive.dataprep)
filtered_data$mods %<>% c("SVDmiss")
filtered_data$nametags %<>% c("Imp")

# Center columns (features)
dimnames.save <- lapply(filtered_data$blocks, dimnames)

filtered_data$blocks[completion_list] %<>% map(~t(apply(.,1,scale, center=TRUE, scale=FALSE)))

dimnames(filtered_data$blocks$DNAm) <- dimnames.save$DNAm
dimnames(filtered_data$blocks$mRNA) <- dimnames.save$mRNA
dimnames(filtered_data$blocks$miRNA) <- dimnames.save$miRNA

filtered_data$mods %<>% c("row-centering")
filtered_data$blocks %>% data.info

saveRDS(filtered_data, paths$section$filter$save_file(filtered_data$nametags))
