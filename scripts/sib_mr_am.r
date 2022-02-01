library(here)
library(R6)
library(tidyverse)
library(simulateGP)
library(MASS)

set.seed(12345)

Family <- R6Class("Family", list(
  
  # Parameters
  nsnp=NA,
  nfam=NA,
  af=NA,
  bxy=NA,
  vg=NA,
  vc=NA,
  bux=NA,
  buy=NA,
  vu=NA,
  mux=NA,
  muy=NA,
  vex=NA,
  vey=NA,
  S=NA,
  spouse_cor=0.8, 
  b = NA,

  # Genotypes
  dads_g = NA,
  mums_g = NA,
  sib1_g = NA,
  sib2_g = NA,

  # Phenotypes
  dads_p = NA,
  mums_p = NA,
  sib1_p = NA,
  sib2_p = NA,

  # Estimates
  pop = NA,
  pop_fe = NA,
  sib_dod = NA,
  sib_fe = NA,
  trio = NA,
  summary = NA,
  snpvar = NA,

  initialize = function(nsnp=100, nfam=10000, af=runif(100, 0.01, 0.99), bxy=0.1, vg=0.4, bux=0.1, buy=0.1, vu=1, vc=1, mux=0, muy=0, vex=1, vey=1, S=1, spouse_cor=0.8)
  {
    ar <- as.list(environment())
    lapply(names(ar), function(arg)
    {
      self[[arg]] <- ar[[arg]]
    })
    self$nsnp <- length(self$af)
    self$choose_effects()
  },

  choose_effects = function(af=self$af, vg=self$vg, S=self$S)
  {
    b <- rnorm(length(af), 0, (2*af*(1-af))^(S/2))
    r <- sqrt(sum(b^2 * 2 * af * (1-af)))
    b <- b * sqrt(vg) / r
    stopifnot(all.equal(sum(b^2 * 2 * af * (1-af)), vg))
    self$b <- b
  },

  hap_to_gen = function(g)
  {
    do.call(cbind, lapply(g, function(x)
    {
      r <- x[,1] + x[,2]
      return(r - mean(r))
    }))
  },

  calc_af = function(g)
  {
    sapply(g, function(x)
    {
      sum(x[,1] + x[,2])/(2*nrow(x))
    }) %>% as.numeric()
  },

  make_gen = function()
  {
    out <- lapply(self$af, function(af)
    {
      cbind(
        rbinom(self$nfam, 1, af),
        rbinom(self$nfam, 1, af)
      )
    })
    return(out)
  },

  make_phen = function(g)
  {
    G <- self$hap_to_gen(g)
    phen <- tibble(
      gbv = G %*% self$b,
      u = rnorm(self$nfam, 0, sqrt(self$vu)),
      c = rnorm(self$nfam, 0, sqrt(self$vc)),
      ex = rnorm(self$nfam, 0, sqrt(self$vex)),
      x = self$mux + gbv + u * self$bux + c + ex,
      ey = rnorm(self$nfam, 0, sqrt(self$vey)),
      y = self$muy + x * self$bxy + u * self$buy + ey
    )
    return(phen)
  },

  make_phen_child = function(g)
  {

    G <- self$hap_to_gen(g)
    phen <- tibble(
      gbv = G %*% self$b,
      u = rnorm(self$nfam, 0, sqrt(self$vu)),
      c = (self$dads_p$c + self$mums_p$c) / sqrt(2),
      ex = rnorm(self$nfam, 0, sqrt(self$vex)),
      x = self$mux + gbv + u * self$bux + c + ex,
      ey = rnorm(self$nfam, 0, sqrt(self$vey)),
      y = self$muy + x * self$bxy + u * self$buy + ey
    )
    return(phen)
  },


  reorder_for_cor = function(x, y)
  {
    d <- tibble(
      x,
      y,
      sx = order(x),
      sy = order(y)
    )
    
    corr_mat <- matrix(c(1, self$spouse_cor, self$spouse_cor, 1), 2, 2)
    mvdat <- mvrnorm(n = nrow(d), mu=c(0,0), Sigma=corr_mat, empirical=TRUE)
    rx <- rank(mvdat[,1], ties.method = "first")
    ry <- rank(mvdat[,2], ties.method = "first")

    return(tibble(id1=d$sx[rx], id2=d$sy[ry]))
  },

  choose_spouses = function(phen1="x", phen2="x")
  {
    o <- self$reorder_for_cor(self$dads_p[[phen1]], self$mums_p[[phen2]])
    self$dads_p <- self$dads_p[o$id1,]
    self$mums_p <- self$mums_p[o$id2,]
    self$dads_g <- lapply(self$dads_g, function(x) x[o$id1,])
    self$mums_g <- lapply(self$mums_g, function(x) x[o$id2,])   
  },

  generate_child = function()
  {
    s1 <- sample(c(TRUE, FALSE), self$nfam, replace=TRUE)
    s2 <- sample(c(TRUE, FALSE), self$nfam, replace=TRUE)
    o <- lapply(1:self$nsnp, function(i)
    {
      x <- matrix(0, self$nfam, 2)
      x[s1, 1] <- self$dads_g[[i]][s1, 1]
      x[!s1, 1] <- self$dads_g[[i]][!s1, 2]
      x[s2, 2] <- self$mums_g[[i]][s2, 1]
      x[!s2, 2] <- self$mums_g[[i]][!s2, 2]
      return(x)
    })
    return(o)
  },

  fe_model = function(g1, g2, y1, y2)
  {
    res <- lapply(1:ncol(g1), function(i){
      g <- c(g1[,i]-mean(g1[,i]), g2[,i]-mean(g2[,i]))
      f <- rep(g1[,i] + g2[,i], times=2)
      y <- c(y1, y2)
      summary(lm(y ~ g+f))
    })

    self$sib_fe <- tibble(
      ahat = sapply(res, function(x) x$coef[1,1]),
      bhat = sapply(res, function(x) x$coef[2,1]),
      se = sapply(res, function(x) x$coef[2,2]),
      pval = sapply(res, function(x) x$coef[2,4])
    )

    self$pop_fe <- tibble(
      ahat = sapply(res, function(x) x$coef[1,1]),
      bhat = sapply(res, function(x) x$coef[3,1]),
      se = sapply(res, function(x) x$coef[3,2]),
      pval = sapply(res, function(x) x$coef[3,4])
    )
  },

  trio_model = function(g1, gm, gd, y1)
  {
    res <- lapply(1:ncol(g1), function(i){
      o <- summary(lm(y1 ~ g1[,i] + gm[,i] + gd[,i]))
      tibble(
        ahat=o$coef[1,1],
        bhat=o$coef[2,1],
        se=o$coef[2,2],
        pval=o$coef[2,4]
      )
    }) %>% bind_rows()
    self$trio <- res
  },

  analysis = function()
  {
    # Perform population estimate
    # Perform sibling estimate using fixed effects model
    g1 <- self$hap_to_gen(self$sib1_g)
    g2 <- self$hap_to_gen(self$sib2_g)
    gd <- self$hap_to_gen(self$dads_g)
    gm <- self$hap_to_gen(self$mums_g)
    self$pop <- simulateGP::gwas(self$sib1_p$x, g1)
    self$sib_dod <- simulateGP::gwas(self$sib1_p$x - self$sib2_p$x, self$hap_to_gen(self$sib1_g) - self$hap_to_gen(self$sib2_g))
    self$fe_model(g1, g2, self$sib1_p$x, self$sib2_p$x)
    self$trio_model(g1, gm, gd, self$sib1_p$x)
    self$af=tibble(af=self$af, af1=self$calc_af(self$sib1_g), af2=self$calc_af(self$sib2_g), afd=self$calc_af(self$dads_g), afm=self$calc_af(self$mums_g))
    self$snpvar = tibble(
      af=self$af$af,
      dadv=sapply(1:ncol(gd), function(i) var(gd[,i])),
      mumv=sapply(1:ncol(gd), function(i) var(gm[,i])),
      sib1v=sapply(1:ncol(gd), function(i) var(g1[,i])),
      sib2v=sapply(1:ncol(gd), function(i) var(g2[,i]))
    )
  },

  summarise = function()
  {
    p <- list()
  p$spouse_cor_hat <- drop(cor(self$mums_p$x, self$dads_p$x))
    p$pop <- lm(self$pop$bhat ~ self$b)$coef[2]
    p$sib_dod <- lm(self$sib_dod$bhat ~ self$b)$coef[2]
    p$sib_fe <- lm(self$sib_fe$bhat ~ self$b)$coef[2]
    p$pop_fe <- lm(self$pop_fe$bhat ~ self$b)$coef[2]
    p$trio <- lm(self$trio$bhat ~ self$b)$coef[2]
    p$sib_vargbv <- var(self$sib1_p$gbv - self$sib2_p$gbv)
    p$sib_covgbv <- cov(self$sib1_p$gbv, (self$sib1_p$gbv + self$sib2_p$gbv))
    p$po_covgbv <- cov(self$sib1_p$gbv, self$dads_p$gbv)
    p$vg_sib <- sum(self$sib_dod$bhat^2 * 2 * self$af$af1 * (1-self$af$af1))
    p$h2_sib <- p$vg_sib / var(self$sib1_p$x)
    p$vg_pop <- sum(self$pop$bhat^2 * 2 * self$af$afd * (1-self$af$afd))
    p$h2_pop <- p$vg_pop / var(self$sib1_p$x)
    p$vg_true <- sum(self$b^2 * 2 * self$af$af * (1-self$af$af))
    p$h2_true <- p$vg_true / var(self$dads_p$x)
    self$summary <- p
  },

  runsim = function()
  {
    self$dads_g <- self$make_gen()
    self$mums_g <- self$make_gen()
    self$dads_p <- self$make_phen(self$dads_g)
    self$mums_p <- self$make_phen(self$mums_g)
    self$choose_spouses("x", "x")
    self$sib1_g <- self$generate_child()
    self$sib2_g <- self$generate_child()
    self$sib1_p <- self$make_phen_child(self$sib1_g)
    self$sib2_p <- self$make_phen_child(self$sib2_g)
    self$analysis()
    self$summarise()
  }
))

##

params <- expand.grid(
  vc=c(0, 1),
  vg=c(0.5, 2),
  spouse_cor=seq(-0.9, 0.9, 0.1),
  sim=1:5
)
dim(params)

res <- lapply(1:nrow(params), function(i)
{
  message(i, " of ", nrow(params))
  p <- params[i, ]
  a <- do.call(Family$new, dplyr::select(p, -sim))
  a$runsim()
  bind_cols(p, as_tibble(a$summary))
}) %>% bind_rows() 
save(res, file=file.path(here(), "results", "am_gwas.rdata"))
