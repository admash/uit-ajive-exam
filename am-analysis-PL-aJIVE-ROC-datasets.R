

library(dplyr)
library(magrittr)
library("ggplot2")
library(pajive)
library(pROC)

paths <- list(
  project = "~/projects/linalg-exam/"
)

setwd(paths$project)

source("am-helper-functions.R")


data_files <- list(
 # Comparison: 5000 vs 10000 mRNA 
 "data-filtered-mv5-CtM-Cm40-Cx3-Imp.RDS",
 "data-filtered-mv10-CtM-Cm40-Cx3-Imp.RDS",

 # Comparison: 40% vs 0% missing cut-off threshold
#"data-filtered-mv5-CtM-Cm40-Cx3-Imp.RDS",
 "data-filtered-mv5-CtM-Cm00-Cx3-Imp.RDS",

 # Comparison: 3 vs 5 for extreme M-value cutoff
#"data-filtered-mv5-CtM-Cm40-Cx3-Imp.RDS",
 "data-filtered-mv5-CtM-Cm40-Cx5-Imp.RDS",
 
 # Comparison: Gene matching vs high variance criterion for DNAm
#"data-filtered-mv5-CtM-Cm40-Cx5-Imp.RDS",
 "data-filtered-mv5-CtV98-Cm40-Cx5-Imp.RDS"#,

 # Comparison: Decreasing subjects
 # Comparison: Decreasing features
 #"data-filtered-*",
)

maxpts = 50
make_filename = function(prefix, tags, suffix){paste0(prefix, paste0(tags,collapse="-"),suffix)}
future::plan(future::multisession(), workers=8)

for(i in seq_along(data_files)){
  d <- readRDS(data_files[[i]])
  blocklist = setdiff(names(d$blocks), "covs")
  set.seed(11)
  pdf(make_filename("plot-PL-",d$nametags, ".pdf"))
  ranks = list()
  for(j in blocklist){
    plot.data <- d$blocks[[j]] %>% get_profile_likelihoods(use.rsvd = T) %>% make_pl_plot.data 
    saveRDS(plot.data, make_filename("data-PL-", c(j,d$nametags), ".RDS"))
    ranks[[j]] <- plot.data %>% filter(i == 2, is.est==1) %>% pull(est) %>% as.character %>% as.numeric()
    g <- ggplot(data=plot.data) +
      geom_line(aes(x = p, y = pl, group = i, colour=i)) + 
      geom_point(aes(x = p, y = pl, group = i, fill = i, colour = i,  alpha = is.est)) + 
      geom_label(aes(x = p, y = pl, label=est), data=plot.data %>% filter(is.est == 1)) +
      labs(title = paste0("Profile Likelihood Plot - ", j),
              subtitle = paste0(d$nametags, collapse=" + "),
           colour = "Removed SVs"
           ) +
      guides(group = "none", alpha = "none", fill = "none") +
      xlab("Rank") + ylab("Profile Likelihood") + 
      theme_minimal() + 
      coord_cartesian(xlim=c(0,maxpts+1))
    print(g)
  }
  dev.off()
  ajr <- pajive::ajive(blocks = d$blocks[blocklist], initial_signal_ranks = unlist(ranks), n_wedin_samples = 100)
  names(ajr$block_decomps) <- blocklist
  saveRDS(ajr, make_filename("data-ajr-", c(j,d$nametags), ".RDS"))
  am_rocplot(ajr, d, d$nametags %>% paste0(collapse="-"))
}



