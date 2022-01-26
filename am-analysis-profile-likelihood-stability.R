library(dplyr)
library(magrittr)
library("ggplot2")
library(pajive)

paths <- list(
  project = "~/projects/linalg-exam/"
)

setwd(paths$project)

source("am-helper-functions.R")

filtered_data <- readRDS("data-filtered-STD.RDS")

filtered_data$blocks %>% data.info

# Singular values - bootstrapping

nsteps = 5
nsamples = 100
nSVs = 50
SVs = matrix(NA, ncol=nSVs, nrow=nsamples)
d = filtered_data$blocks$miRNA
rm(all.plot.data)
for(j in 1:nsteps){ 
  for(i in 1:nsamples){
      if(i %% 10 == 1) message(i, appendLF = F) else message(".", appendLF = F)
      #s = svd( d[ sample(nrow(d), replace = T), 1:(100+(j-1)*20)])
      s = svd( d[ sample(nrow(d), replace = T), ])
      #message(length(s$d))
      SVs[i,] <-s[['d']][1:50] 
  }
  message("")
  plot.data <- data.frame(
    g = 100+(j-1)*20,
    i = map(1:nsamples, ~rep(., nSVs)) %>% unlist,
    p = rep(1:nSVs, nsamples),
    sv = SVs %>% t %>% as.vector )
  if("all.plot.data" %in% ls()) { all.plot.data %<>% rbind(plot.data)}
  else { all.plot.data = plot.data}
}

all.plot.data$g %<>% factor

pdf("miRNA-sv-bootstrap-box-plot.pdf")
ggplot(data=all.plot.data %>% filter(p < 15)) +
  #geom_line(aes(x = p, y = (sv), group = i,col=g),alpha = 0.01)
  geom_boxplot(aes(colour = factor(g), x = factor(p), y = (sv)))
dev.off()


  for(i in 1:nsamples){
      if(i %% 10 == 1) message(i, appendLF = F) else message(".", appendLF = F)
      #s = svd( d[ sample(nrow(d), replace = T), 1:(100+(j-1)*20)])
      s = svd( d[ sample(nrow(d), replace = T), ])
      #message(length(s$d))
      SVs[i,] <-s[['d']][1:50] 
  }
  message("")
  plot.data <- data.frame(
    g = 100+(j-1)*20,
    i = map(1:nsamples, ~rep(., nSVs)) %>% unlist,
    p = rep(1:nSVs, nsamples),
    sv = SVs %>% t %>% as.vector )


pdf("miRNA-sv-bootstrap-line-plot.pdf")
ggplot(data=plot.data %>% filter(p < 15)) +
  geom_line(aes(x = p, y = (sv), group = i),alpha = 0.05) + 
  geom_boxplot(aes(x = p, group = p, y = (sv)))
dev.off()

############
# Profile Likelihood

d <- filtered_data$blocks$DNAm

all.r = list()
all.r$DNAm = r



all.r$mRNA = get_profile_likelihoods(filtered_data$blocks$mRNA)
all.r$miRNA = get_profile_likelihoods(filtered_data$blocks$miRNA)
all.r$miRNA.r = get_profile_likelihoods(filtered_data$blocks$miRNA, use.rsvd = T)
r = all.r$miRNA.r

make_pl_plot.data <- function(r){
  data.frame(
    i = map(r, ~rep(.$removed, length(.$pl))) %>% unlist %>% factor,
    p = map(r, ~1:length(.$pl)) %>% unlist ,
    pl = map(r, ~log(-.$pl)) %>% { do.call(rbind,.)} %>% t %>% as.vector,
    #est = map(remove_points, ~rep(r[[.+1]]$rank.est, min(dim(d)) )) %>% unlist,
    est = map(r, ~rep(.$rank.est, length(.$pl))) %>% unlist %>% factor,
    is.est = map(r, ~if_else((1:length(.$pl)) == .$rank.est,1,0,0) )  %>% unlist
#    is.est = map(remove_points, ~ if_else((1:min(dim(d))) == r[[.+1]]$rank.est,1,0,0) )  %>% unlist
    )
}


pdf("miRNA-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data) +
  geom_line(aes(x = p, y = pl, group = i, colour=i)) + 
  geom_point(aes(x = p, y = pl, group = i, fill = i, colour = i,  alpha = is.est)) + 
  geom_label(aes(x = p, y = pl, label=est), data=plot.data %>% filter(is.est == 1)) +
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()

##### svd vs rvsd

r.T = get_profile_likelihoods(filtered_data$blocks$mRNA, use.rsvd = T)
r.F = get_profile_likelihoods(filtered_data$blocks$mRNA, use.rsvd = F)
plot.data <- rbind(
  data.frame(rsvd = T, make_pl_plot.data(r.T)) ,
  data.frame(rsvd = F, make_pl_plot.data(r.F)) 
)

pdf("miRNA-rsvd-comp-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data) +
  geom_line(aes(x = p, y = pl, group = paste0(rsvd,i), colour=paste0(rsvd,i), lty=rsvd)) + 
  geom_point(aes(x = p, y = pl, group = i, fill = i, colour = i,  alpha = is.est)) + 
  geom_label(aes(x = p, y = pl, label=est), data=plot.data %>% filter(is.est == 1)) +
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()


#### mRNA profile-likelihood.pdf ####
plot.data <- filtered_data$blocks$mRNA %>%
  get_profile_likelihoods( use.rsvd = T) %>%
  make_pl_plot.data

pdf("mRNA-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data) +
  geom_line(aes(x = p, y = pl, group = i, colour=i)) + 
  geom_point(aes(x = p, y = pl, group = i, fill = i, colour = i,  alpha = is.est)) + 
  geom_label(aes(x = p, y = pl, label=est), data=plot.data %>% filter(is.est == 1)) +
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()


#### DNAm profile-likelihood.pdf ####
plot.data <- filtered_data$blocks$DNAm %>%
  get_profile_likelihoods( use.rsvd = T) %>%
  make_pl_plot.data

pdf("DNAm-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data) +
  geom_line(aes(x = p, y = pl, group = i, colour=i)) + 
  geom_point(aes(x = p, y = pl, group = i, fill = i, colour = i,  alpha = is.est)) + 
  geom_label(aes(x = p, y = pl, label=est), data=plot.data %>% filter(is.est == 1)) +
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()


#### bootstrap subject ####
future::plan(future::multisession(), workers=16)
plot.data.b <- future.apply::future_lapply(1:100, function(x){
  set.seed(x)
  filtered_data$blocks$mRNA[sample(nrow(filtered_data$blocks$mRNA), 
                                   replace=T),] %>%
    get_profile_likelihoods(use.rsvd = T) %>%
    make_pl_plot.data %>%
    mutate(round = x)
}, future.packages = c("rsvd", "dplyr")) 

plot.data <- do.call(rbind, plot.data.b)

pdf("mRNA-bootstrap-subject-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data %>% mutate(gr = paste0(round,",",i))) +
  geom_line(aes(x = p, y = pl, group = gr, colour=i), alpha = 0.2) + 
  geom_point(aes(x = p, y = pl, group = gr, colour=i), size=0.5, fill = "black",
             position = position_jitter(width=0.1, height=0.1),
            data=plot.data %>% mutate(gr = paste0(round,",",i)) %>% filter(is.est==T)) + 
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()


#### random subject subset ####
future::plan(future::multisession(), workers=16)
plot.data.b <- future.apply::future_lapply(1:100, function(x){
  set.seed(x)
  filtered_data$blocks$mRNA[sample(nrow(filtered_data$blocks$mRNA), 
                                   round(0.63*nrow(filtered_data$blocks$mRNA)), 
                                   replace=F),] %>%
    get_profile_likelihoods(use.rsvd = T) %>%
    make_pl_plot.data %>%
    mutate(round = x)
}, future.packages = c("rsvd", "dplyr"))

plot.data <- do.call(rbind, plot.data.b)

pdf("mRNA-random-subject-subset-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data %>% mutate(gr = paste0(round,",",i))) +
  geom_line(aes(x = p, y = pl, group = gr, colour=i), alpha = 0.2) + 
  geom_point(aes(x = p, y = pl, group = gr, colour=i), size=0.5, fill = "black",
             position = position_jitter(width=0.1, height=0.1),
            data=plot.data %>% mutate(gr = paste0(round,",",i)) %>% filter(is.est==T)) + 
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()


#### random feature subsets ####
future::plan(future::multisession(), workers=16)
plot.data.r <- future.apply::future_lapply(1:100, function(x){
  set.seed(x)
  filtered_data$blocks$mRNA[,sample(ncol(filtered_data$blocks$mRNA), 
                                    round(0.9*ncol(filtered_data$blocks$mRNA)), 
                                    replace=F)] %>%
    get_profile_likelihoods(use.rsvd = T) %>%
    make_pl_plot.data %>%
    mutate(round = x)
}, future.packages = c("rsvd", "dplyr"))

plot.data <- do.call(rbind, plot.data.r)

pdf("mRNA-random-feature-subset-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data %>% mutate(gr = paste0(round,",",i))) +
  geom_line(aes(x = p, y = pl, group = gr, colour=i), alpha = 0.2) + 
  geom_point(aes(x = p, y = pl, group = gr, colour=i), size=0.5, fill = "black",
             position = position_jitter(width=0.1, height=0.1),
            data=plot.data %>% mutate(gr = paste0(round,",",i)) %>% filter(is.est==T)) + 
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()

#### decreasing subject subsets ####
future::plan(future::multisession(), workers=16)
plot.data.r <- future.apply::future_lapply(seq(30,60,2), function(x){
  set.seed(5)
  filtered_data$blocks$mRNA[sample(nrow(filtered_data$blocks$mRNA), 
                                    replace=F)[1:ceiling(x/60*nrow(filtered_data$blocks$mRNA))],] %>%
    get_profile_likelihoods(use.rsvd = T) %>%
    make_pl_plot.data %>%
    mutate(round = x)
}, future.packages = c("rsvd", "dplyr"))

plot.data <- do.call(rbind, plot.data.r)



pdf("mRNA-decreasing-subject-subsets-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data %>% mutate(gr = paste0(round,",",i))) +
  geom_line(aes(x = p, y = pl, group = gr, colour=i), alpha = 0.2) + 
  geom_point(aes(x = p, y = pl, group = gr, colour=i), size=0.5, fill = "black",
#             position = position_jitter(width=0.1, height=0.1),
            data=plot.data %>% mutate(gr = paste0(round,",",i)) %>% filter(is.est==T)) + 
  geom_line(aes(x = p, y = pl, group = i), alpha = 0.5,
            data=plot.data %>% filter(is.est==T)) + 
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()


#### decreasing feature subsets ####
future::plan(future::multisession(), workers=16)
plot.data.r <- future.apply::future_lapply(seq(30,60,2), function(x){
  set.seed(11)
  filtered_data$blocks$mRNA[,sample(ncol(filtered_data$blocks$mRNA), 
                                    replace=F)[1:ceiling(x/60*nrow(filtered_data$blocks$mRNA))]] %>%
    get_profile_likelihoods(use.rsvd = T) %>%
    make_pl_plot.data %>%
    mutate(round = x)
}, future.packages = c("rsvd", "dplyr"))

plot.data <- do.call(rbind, plot.data.r)



pdf("mRNA-decreasing-feature-subsets-profile-likelihood.pdf")
maxpts = 50
ggplot(data=plot.data %>% mutate(gr = paste0(round,",",i))) +
  geom_line(aes(x = p, y = pl, group = gr, colour=i), alpha = 0.2) + 
  geom_point(aes(x = p, y = pl, group = gr, colour=i), size=0.5, fill = "black",
#             position = position_jitter(width=0.1, height=0.1),
            data=plot.data %>% mutate(gr = paste0(round,",",i)) %>% filter(is.est==T)) + 
  geom_line(aes(x = p, y = pl, group = i), alpha = 0.5,
            data=plot.data %>% filter(is.est==T)) + 
  coord_cartesian(xlim=c(0,maxpts+1))
dev.off()
