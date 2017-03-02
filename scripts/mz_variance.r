library(dplyr)
library(ggplot2)

fastAssoc <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p, n=n
	))
}


make_sibs <- function(n, maf, eff, vA, vC, model)
{
	g <- rbinom(n, 2, maf)
	A <- rnorm(n, sd=sqrt(vA))
	C <- rnorm(n, sd=sqrt(vC))
	E1 <- rnorm(n, sd=sqrt(1 - vA - vC))
	E2 <- rnorm(n, sd=sqrt(1 - vA - vC))
	if(model == 1)
	{
		z1 <- rnorm(n, mean=0, sd=sqrt(eff * g))
		z2 <- rnorm(n, mean=0, sd=sqrt(eff * g))
		y1 <- scale(A + z1 + C + E1)
		y2 <- scale(A + z2 + C + E2)
	}
	if(model == 2)
	{
		y1 <- scale(A + scale(g) * sqrt(eff) + C + E1)
		y2 <- scale(A + scale(g) * sqrt(eff) + C + E2)
	}
	if(model == 3)
	{
		z1 <- rnorm(n, mean=0, sd=sqrt(eff * g))
		z2 <- rnorm(n, mean=0, sd=sqrt(eff * g))
		y1 <- scale(A + scale(g) * sqrt(eff) + z1 + C + E1)
		y2 <- scale(A + scale(g) * sqrt(eff) + z2 + C + E2)
	}

	dat <- data.frame(
		y1 = y1,
		y2 = y2,
		A = A,
		C = C,
		g = g,
		D2 = (y1 - y2)^2,
		S2 = (y1 + y2)^2
	)
	return(dat)
}


do_test <- function(dat)
{
	res <- rbind(
		as.data.frame(fastAssoc(dat$D2, dat$g)),
		as.data.frame(fastAssoc(dat$S2, dat$g)),
		as.data.frame(fastAssoc(dat$y1, dat$g)),
		as.data.frame(fastAssoc(dat$y1^2, dat$g))
	)
	res$test <- c("D2", "S2", "y", "ysq")
	return(res)
}


###


arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
chunks <- as.numeric(arguments[2])
out <- arguments[3]

param <- expand.grid(
	nsim = 1:100,
	n = c(10000, 20000, 30000, 40000, 50000),
	eff = seq(0, 0.035, by=0.001),
	maf = c(0.05, 0.15, 0.4),
	A = 0.2,
	C = c(0.25, 0.5, 0.75),
	model = 1:3
)

chunksize <- ceiling(nrow(param) / chunks)
t1 <- (jid - 1) * chunksize + 1
t2 <- min(jid * chunksize, nrow(param))

message("total size: ", nrow(param))
message("running: ", t1, " to ", t2)

param <- param[t1:t2, ]

res <- list()
for(i in 1:nrow(param))
{
	message(i, " of ", nrow(param))
	dat <- make_sibs(
		param$n[i], 
		param$maf[i], 
		param$eff[i], 
		param$A[i], 
		param$C[i], 
		param$model[i]
	)
	right <- do_test(dat)
	left <- param[rep(i, nrow(right)), ]
	res[[i]] <- cbind(right, left)
}

res <- bind_rows(res)
save(res, file=out)

