library(rsvd)
library(purrr)


as.data.matrix <- function(d){
  model.matrix(
    paste0( "~",  paste(  names(d),   collapse = " + ")) %>% as.formula,  
    data = d) %>% as.matrix
}


data.info <- function(blocks){
  map_df(
    blocks,
    ~ c( 
        case_when(                          # Block Types:
          is.matrix(.) ~ "matrix",          #  matrix
          is.data.frame(.) ~ "data.frame",  #  data.frame
          TRUE ~ "other"),                  #  other(?!)
        dim(.) ,                            # Dimensions
        dim(.)[1] * dim(.)[2],              # Number of values
        sum(is.na(.)),                      # Number of missing values
        sum(colSums(is.na(.)) > 0),         # Number of columns with missing values
        sum(rowSums(is.na(.)) > 0),         # Number of rows with missing values
        dimnames(.)[[1]][1],                # First row name
        dimnames(.)[[2]][1]                 # First column name
       ) %>% 
      set_names("type", "rows", "cols", "values", "missValues", "missCols", "missRows", "rownames", "colnames"), .id="block")
}


estimate_init_rank  <- function(singular.values, plot=FALSE, lines.only=FALSE){
  
  # calculate profile loglikelihood utility function
  l <- length(singular.values)
  proflik <- c()
  for (i in 1:l){
    mu1 <- mean(singular.values[1:i])
    s1 <- sum((singular.values[1:i]- mu1)^2)
    mu2 <- mean(singular.values[(i+1):l])
    s2 <- sum((singular.values[(i+1):l]-mu2)^2)
    if (i == l) s2 <- 0
    proflik[i] <- s1+s2
  }
  # visualize results
  if(plot==TRUE){
    if(lines.only == TRUE){
      print(lines(-proflik))
    } else {
      print(plot(-proflik,))
      print(lines(-proflik))
    }
  }

  # rank will be arg min 
  # +1 because we left one out
  return( list(rank.est=which.min(proflik), pl = -proflik))
} 
  

estimate_init_ranks <- function(singular.values.list, remove){
  
  if(length(singular.values.list) != length(remove)) { stop("Length of singular values list must match number of singular values to remove")}
  
  r <- map2_int(singular.values.list, remove, ~ estimate_init_rank(.x[seq_along(.x) > .y]))
    
  # add eliminated ones
  #r <- r + c(2,2,1)
  cat("  ", paste0(r, collapse=", "), "Profile likelihood-based ranks sub-estimate", "\n")
  r <- r + remove
  cat("+ ", paste0(remove, collapse=", "), "Removed singular values", "\n")
  cat("==========\n")
  cat("  ", paste0(r, collapse=", "), "Initial ranks estimate", "\n")
  return(r)
  
}


SVDmiss <- function(X, niter=25, ncomp=min(4,dim(X)[2]), conv.reldiff=0.001, use.rsvd=FALSE, p=10, q=2)
{
  ##ensure at least one iteration
  niter <- max(niter,1)
  ##and ensure sane number of components
  if( ncomp<1 ){
    stop("ncomp should be >0, is: ", ncomp)
  }
  if(use.rsvd == TRUE) { svdfunc = rsvd } else { svdfunc = svd }
  
  ##First find missing values
  Ina <- is.na(X)
  if( all(!Ina) ){
    ##if X has no missing data this is simple
    svd0 <- svd(X)
    XF <- X
    i <- diff <- reldiff <- 0
  }else{
    ##X has missing data, use iterative method
    ##Iterative svd calculation with missing data.
    ##Initial first element of U matrix is average curve.
    message("Calculating means ... ", appendLF = FALSE)
    u1 <- rowMeans(X, na.rm = TRUE)
    message("done.")
    XM <- matrix(1, nrow(X), ncol(X))
    XM[Ina] <- 0
    XZ <- X
    # v1 is proportional to X'u1/(u1'u1), but with calculations
    # filling in zeros for missing values in the sums.
    XZ[Ina] <- 0.
    # Fill in missing values for initial complete SVD calculation.
    # Then iterate  using only complete data.
    message("Building v1 ... ", appendLF = FALSE)
      # n = 7000
      # m = 230
      # X = matrix(rnorm(n*m), ncol=n)
      # system.time(x1 <- (diag(t(X) %*% X))) #149 seconds
      # system.time(x2 <- (colSums(X*X)))  #1.09 seconds
      # all.equal(x1, x2)
    #v1 <- diag(t(XZ) %*% (XM * u1))/diag(t(XM * u1) %*% (XM * u1))
    v1 <- colSums((XZ) * (XM * u1))/colSums((XM * u1) * (XM * u1))
    message("done.")
    XF <- X
    XF[Ina] <- (matrix(u1, ncol = 1) %*% matrix(v1, nrow = 1))[Ina]
    if( any(is.na(XF)) )
      stop("Unable to complete matrix, too much missing data")
    reldiff <- conv.reldiff+1
    i <- 0
    message("iteration: reldiff")
    while(i<niter && reldiff>conv.reldiff){
      if(use.rsvd == TRUE)
        svd0 <- rsvd(XF, k=ncomp, p=p, q=q)
      else 
        svd0 <- svd(XF)
      Xnew <- X
      Xnew[Ina] <- (svd0$u[, 1:ncomp] %*%
                      diag(svd0$d[1:ncomp],nrow=length(svd0$d[1:ncomp])) %*%
                      t(svd0$v[,1:ncomp]))[Ina]
      diff <- max(abs(Xnew - XF))
      reldiff <- diff/max(abs(XF[Ina]))
      XF <- Xnew
      message(sprintf(" %02d: %0.6f", i, reldiff))
      i <- i+1
    }
    message("done.")
  }#if( all(!is.na(X)) ) ... else ...
  final.diff <- c(diff,reldiff,i,niter)
  names(final.diff) <- c("diff","rel.diff","n.iter","max.iter")
  return(list(Xfill=XF, status=final.diff))
}##function SVDmiss


get_profile_likelihoods <- function(d, remove_points=0:10, use.rsvd = FALSE, q = 10){
  r <- lapply(remove_points,function(s){
    message(s)
    if(use.rsvd) { svs =  svd(d, nu=0, nv=0      )[["d"]] }
    else         { svs = rsvd(d, nu=0, nv=0, q=10)[["d"]] }
    if(s == 0){
      est = estimate_init_rank(svs, plot=FALSE)
    } else {
      est = estimate_init_rank(svs[-(1:s)], plot=FALSE)
    }
    est$rank.est = est$rank.est + s
    est$pl = c(rep(NA, s), est$pl)
    est$removed = s
    return(est)
  } )
  return(r)
}

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

am_rocplot <- function(ajRes, d, output_suffix){

 blocks = names(ajRes$block_decomps)

# Get joint and individual ranks and scores
# We only want 5 individual each
num_indiv_ranks = 5

     ranks = list()
idv.scores = list()
       pcs = list()

 ranks$joint <- ajRes$joint_rank
  jnt.scores <- ajRes$joint_scores %>% data.frame

for(b in blocks){
       ranks[[b]] <- min(ajRes$block_decomps[[b]]$individual$rank, num_indiv_ranks)
  idv.scores[[b]] <- 
    if(ranks[[b]]==1){ 
      ajRes$block_decomps[[b]]$individual$u[,1:ranks[[b]]] %>% data.frame %>% `colnames<-`(paste0(b,".X1"))
    } else {
      ajRes$block_decomps[[b]]$individual$u[,1:ranks[[b]]] %>% data.frame
    }
         pcs[[b]] <- prcomp(d$blocks[[b]], center = TRUE, scale = TRUE, rank. = num_indiv_ranks)$x
}  

# Covariate data frame with convenient names
covs <- d$blocks$covs %>%  data.frame %>% select(age=age.sample, bmi=BMI, smoking=smoking_status_cat)
covs$smoking %<>% factor %>% relevel(ref="Never")
covs$bmi %<>% as.numeric
covs$age %<>% as.integer


# Construct model data frame
model.data <- data.frame(y = d$blocks$covs[,"case_ctrl"] %>% factor) %>% 
  cbind(covs, jnt.scores, idv.scores, pcs)


# Define models
models <- list(
  list( 
    name = "Joint and Individual Components",
    covars = c(names(covs),
               paste0("X", 1:ranks$joint),
               map(blocks, ~ paste0(., ".X", 1:ranks[[.]])) %>% unlist),
    color = "black"), 
  
  list(
    name = paste0("Joint Components", " [", ranks$joint,"]"),
    covars = c(names(covs), 
               paste0("X", 1:ranks$joint)), 
    color = "blue"),  
  
  list(
    name = paste0("Individual Components"," [",paste0(ranks[blocks],"/",
                                                      map_int(ajRes$block_decomps[blocks], ~ .$individual$rank)
                                                      , collapse=","),"]"),
    covars = c(names(covs), 
               map(blocks, ~ paste0(., ".X", 1:ranks[[.]])) %>% unlist),
    color = "green"),
  
  list(
    name = paste0("Non-integrative"," [",paste0(rep(num_indiv_ranks, length(blocks)), collapse=","),"]"),
    covars = c(names(covs), 
               map(blocks, ~ paste0(., ".PC", 1:num_indiv_ranks)) %>% unlist),
    color = "magenta"),
  
  list(
    name = "Patient Covariates",
    covars = c(names(covs)), 
    color = "red") 
)

names(models) <- map_chr(models, ~ .$name)

# Build formulas
models %<>% 
  lapply(function(m){ 
    m$exposure = "y" 
    m$formula = paste0(m$exposure, "~", paste(m$covars, collapse="+"))
    m$formula %<>% as.formula 
    return(m) })

# in-sample ROC
models %<>% 
  lapply(function(m){ 
    m$fit  = glm(m$formula, family="binomial", data=model.data)
    m$pred = predict(m$fit, model.data, type = 'response')
    m$roc  = roc(model.data$y, m$pred)
    m$auc  = auc(m$roc)
    m$label = sprintf("(%0.3f) %s", m$auc, m$name)
    return(m)
      })

roc.plot <- models %>% 
  lapply(function(m){ m$roc }) %>% 
  ggroc + 
  scale_color_discrete("(AUC) Model", labels=map(models, ~ .$label)) +
  theme_minimal() + 
  theme(legend.position = c(0.65, 0.25))

pdf(paste0("plot-ROC-",output_suffix,".pdf"))
print(roc.plot)
dev.off()

}

am_idv_rocplot <- function(ajRes, d, output_suffix){

 blocks = names(ajRes$block_decomps)

# Get joint and individual ranks and scores
# We only want 5 individual each
num_indiv_ranks = 5

     ranks = list()
idv.scores = list()
       pcs = list()

 ranks$joint <- ajRes$joint_rank
  jnt.scores <- ajRes$joint_scores %>% data.frame

for(b in blocks){
       ranks[[b]] <- min(ajRes$block_decomps[[b]]$individual$rank, num_indiv_ranks)
  idv.scores[[b]] <- 
    if(ranks[[b]]==1){ 
      ajRes$block_decomps[[b]]$individual$u[,1:ranks[[b]]] %>% data.frame %>% `colnames<-`(paste0(b,".X1"))
    } else {
      ajRes$block_decomps[[b]]$individual$u[,1:ranks[[b]]] %>% data.frame
    }
         pcs[[b]] <- prcomp(d$blocks[[b]], center = TRUE, scale = TRUE, rank. = num_indiv_ranks)$x
}  

# Covariate data frame with convenient names
covs <- d$blocks$covs %>%  data.frame %>% select(age=age.sample, bmi=BMI, smoking=smoking_status_cat)
covs$smoking %<>% factor %>% relevel(ref="Never")
covs$bmi %<>% as.numeric
covs$age %<>% as.integer


# Construct model data frame
model.data <- data.frame(y = d$blocks$covs[,"case_ctrl"] %>% factor) %>% 
  cbind(covs, jnt.scores, idv.scores, pcs)


# Define models
models <- list(
  list( 
    name = "Joint and Individual Components",
    covars = c(names(covs),
               paste0("X", 1:ranks$joint),
               map(blocks, ~ paste0(., ".X", 1:ranks[[.]])) %>% unlist),
    color = "black"), 
  
  list(
    name = paste0("Joint Components", " [", ranks$joint,"]"),
    covars = c(names(covs), 
               paste0("X", 1:ranks$joint)), 
    color = "blue"),  
  
  list(
    name = paste0("Individual Components"," [",paste0(ranks[blocks],"/",
                                                      map_int(ajRes$block_decomps[blocks], ~ .$individual$rank)
                                                      , collapse=","),"]"),
    covars = c(names(covs), 
               map(blocks, ~ paste0(., ".X", 1:ranks[[.]])) %>% unlist),
    color = "green"),
  
  list(
    name = paste0("Non-integrative"," [",paste0(rep(num_indiv_ranks, length(blocks)), collapse=","),"]"),
    covars = c(names(covs), 
               map(blocks, ~ paste0(., ".PC", 1:num_indiv_ranks)) %>% unlist),
    color = "magenta"),
  
  list(
    name = "Patient Covariates",
    covars = c(names(covs)), 
    color = "red") 
)

models <- c(
  lapply(c("miRNA", "mRNA", "DNAm"), function(blockname){
    list(
      name = paste0(blockname, " (PCA:", num_indiv_ranks, ")"),
      covars = c(names(covs),  paste0(blockname, ".PC", 1:num_indiv_ranks)),
      type="pca"
    )
  }),
  lapply(c("miRNA", "mRNA", "DNAm"), function(blockname){
    list(
      name = paste0(blockname, " (aJIVE:", num_indiv_ranks, ")"),
      covars = c(names(covs),  paste0(blockname, ".X", 1:num_indiv_ranks)),
      type="ajive"
    )
  })
)

names(models) <- map_chr(models, ~ .$name)

# Build formulas
models %<>% 
  lapply(function(m){ 
    m$exposure = "y" 
    m$formula = paste0(m$exposure, "~", paste(m$covars, collapse="+"))
    m$formula %<>% as.formula 
    return(m) })

# in-sample ROC
models %<>% 
  lapply(function(m){ 
    m$fit  = glm(m$formula, family="binomial", data=model.data)
    m$pred = predict(m$fit, model.data, type = 'response')
    m$roc  = roc(model.data$y, m$pred)
    m$auc  = auc(m$roc)
    m$label = sprintf("(%0.3f) %s", m$auc, m$name)
    return(m)
      })

roc.plot <- models %>% 
  lapply(function(m){ m$roc }) %>% 
  ggroc + 
  scale_color_discrete("(AUC) Model", labels=map(models, ~ .$label)) +
  theme_minimal() + 
  theme(legend.position = c(0.65, 0.25))

pdf(paste0("plot-ROC-",output_suffix,".pdf"))
print(roc.plot)
dev.off()

}
