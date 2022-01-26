library(dplyr)
library(magrittr)
library("ggplot2")
library(pajive)

setwd("~/projects/linalg-exam/")

source("am-helper-functions.R")

filtered_data <- readRDS("filtered_data-mRNAvar.5000-DNAmXmRNA-NoMissDNAm-DNAmLT3.RDS")

filtered_data$blocks %>% data.info


#=-.-=-.-=-.-=-.-=-.-=-.-=-.-=#
#####  ALL - scree plots  #####
#=-.-=-.-=-.-=-.-=-.-=-.-=-.-=#

#source("functions/Screeplots.R")
#screeplot(filtered_data$blocks) # 10, 7, 12

#=-.-=-.-=-.-=-.-=-.-=-.-=-.-=-.-=-.-=#
#####  ALL - profile likelihood   #####
#=-.-=-.-=-.-=-.-=-.-=-.-=-.-=-.-=-.-=#

# use profile likelihood
# SVD on each source
system.time(singular.values1 <- svd(filtered_data$blocks$miRNA)[['d']]) # 0.042 sec
system.time(singular.values2 <- svd(filtered_data$blocks$mRNA )[['d']]) # 1.719 sec
system.time(singular.values3 <- svd(filtered_data$blocks$DNAm )[['d']]) # 20.22 sec
# check first if we want to eliminate first 1 or 2 
# if too extremes
singular.delta1 <- singular.values1 %>% round() %>% add(1) %>% log() %>% diff() %>% abs 
singular.delta2 <- singular.values2 %>% round() %>% add(1) %>% log() %>% diff() %>% abs 
singular.delta3 <- singular.values3 %>% round() %>% add(1) %>% log() %>% diff() %>% abs 
remove1 = max(which((singular.delta1 > 1)))
remove2 = max(which((singular.delta2 > 1)))
remove3 = max(which((singular.delta3 > 1)))


singular.val <- list(singular.values1, 
                     singular.values2, 
                     singular.values3)
pdf("init_ranks.pdf")
#r = estimate_init_ranks(singular.val, c(remove1, remove2, remove3))
r = estimate_init_ranks(singular.val, c(remove1, remove2, remove3))
r = list()
for(i in 0:5) #for(j in 1:4) for(k in 1:4)
#r[[paste0(i,j,k)]] <- map2_int(singular.val, c(i,j,k), ~ estimate_init_rank(.x[seq_along(.x) > .y]))
r[[paste0(i,i,i)]] <- map2_int(singular.val, c(i,i,i), ~ estimate_init_rank(.x[seq_along(.x) > .y]))
dev.off()

future::plan(future::multisession(), workers=16)

# 2774 seconds serial, 589 seconds parallel - old ajive data
# 198/5000/5000 - s=534,550;p=207 seconds
#for(num_samples in c(100, 200, 500, 600, 800, 1000)){
timing  <- list(pajive = list(), ajive = list())
results <- list(pajive = list(), ajive = list())
for(num_samples in c(100, 200, 1000)){
  cat("initial_signal_ranks:", r, "\n")
  timing$pajive[[as.character(num_samples)]] <- 
    system.time(
      results$pajive[[as.character(num_samples)]] <- pajive::ajive(filtered_data$blocks, initial_signal_ranks = r, n_wedin_samples = num_samples)) 
  timing$ajive[[ as.character(num_samples)]] <- 
    system.time( 
      results$ajive[[as.character(num_samples)]] <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r, n_wedin_samples = num_samples)) 
#  saveRDS(ajiveResults, paste0("ajiveResults-rsvd-wedin-samples-", num_samples, ".RDS"))
}

      results$ajive32 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`222`, joint_rank = 3)
      results$ajive33 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`333`, joint_rank = 3)
      results$ajive34 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`444`, joint_rank = 3)
      results$ajive42 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`222`, joint_rank = 4)
      results$ajive43 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`333`, joint_rank = 4)
      results$ajive44 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`444`, joint_rank = 4)
      results$ajive52 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`222`, joint_rank = 5)
      results$ajive53 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`333`, joint_rank = 5)
      results$ajive54 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r$`444`, joint_rank = 5)

      results$ajive5 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r, joint_rank = 5)
      results$ajive56 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r, joint_rank = 5)
      results$ajive57 <-  ajive::ajive(filtered_data$blocks, initial_signal_ranks = r, joint_rank = 5)

ajr100 <- readRDS("ajiveResults-wedin-samples-100.RDS")
ajr100r <- readRDS("ajiveResults-rsvd-wedin-samples-100.RDS")
ajr200 <- readRDS("ajiveResults-wedin-samples-200.RDS")
ajr200r <- readRDS("ajiveResults-rsvd-wedin-samples-200.RDS")
ajr500 <- readRDS("ajiveResults-wedin-samples-500.RDS")
ajr500r <- readRDS("ajiveResults-rsvd-wedin-samples-500.RDS")
ajr600 <- readRDS("ajiveResults-wedin-samples-600.RDS")
ajr600r <- readRDS("ajiveResults-rsvd-wedin-samples-600.RDS")
ajr800 <- readRDS("ajiveResults-wedin-samples-800.RDS")
ajr800r <- readRDS("ajiveResults-rsvd-wedin-samples-800.RDS")
ajr1000 <- readRDS("ajiveResults-wedin-samples-1000.RDS")
ajr1000r <- readRDS("ajiveResults-rsvd-wedin-samples-1000.RDS")
pdf("comp.pdf")
plot(density(ajr500$joint_rank_sel$wedin$wedin_samples %>% log), type="n")
lines(density(ajr500$joint_rank_sel$wedin$wedin_samples %>% log ), col="green")
lines(density(ajr500r$joint_rank_sel$wedin$wedin_samples %>% log ), col="green", lty=2)
lines(density(ajr100$joint_rank_sel$wedin$wedin_samples %>% log), col="red",)
lines(density(ajr100r$joint_rank_sel$wedin$wedin_samples %>% log), col="red", lty=2)
lines(density(ajr200$joint_rank_sel$wedin$wedin_samples %>% log), col="orange",)
lines(density(ajr200r$joint_rank_sel$wedin$wedin_samples %>% log), col="orange", lty=2)
#lines(density(ajr1000$joint_rank_sel$wedin$wedin_samples), col="blue",)
lines(density(ajr1000r$joint_rank_sel$wedin$wedin_samples %>% log), col="blue", lty=2)
lines(density(ajr800$joint_rank_sel$wedin$wedin_samples %>% log), col="magenta",)
lines(density(ajr800r$joint_rank_sel$wedin$wedin_samples %>% log), col="magenta", lty=2)
dev.off()
