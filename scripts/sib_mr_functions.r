library(ggplot2)
library(reshape2)
library(TwoSampleMR)
library(dplyr)


###########

###########

make_families <- function(af, nfam)
{
	nsnp <- length(af)
	dads <- matrix(0, nfam, nsnp)
	mums <- matrix(0, nfam, nsnp)
	sibs1 <- matrix(0, nfam, nsnp)
	sibs2 <- matrix(0, nfam, nsnp)
	ibd <- matrix(0, nfam, nsnp)
	ibs <- matrix(0, nfam, nsnp)
	for(i in 1:nsnp)
	{
		dad1 <- rbinom(nfam, 1, af[i]) + 1
		dad2 <- (rbinom(nfam, 1, af[i]) + 1) * -1
		mum1 <- rbinom(nfam, 1, af[i]) + 1
		mum2 <- (rbinom(nfam, 1, af[i]) + 1) * -1

		dadindex <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
		dadh <- rep(NA, nfam)
		dadh[dadindex] <- dad1[dadindex]
		dadh[!dadindex] <- dad2[!dadindex]

		mumindex <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
		mumh <- rep(NA, nfam)
		mumh[mumindex] <- mum1[mumindex]
		mumh[!mumindex] <- mum2[!mumindex]

		sib1 <- cbind(dadh, mumh)

		dadindex <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
		dadh <- rep(NA, nfam)
		dadh[dadindex] <- dad1[dadindex]
		dadh[!dadindex] <- dad2[!dadindex]

		mumindex <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
		mumh <- rep(NA, nfam)
		mumh[mumindex] <- mum1[mumindex]
		mumh[!mumindex] <- mum2[!mumindex]

		sib2 <- cbind(dadh, mumh)

		ibd[,i] <- (as.numeric(sib1[,1] == sib2[,1]) + as.numeric(sib1[,2] == sib2[,2])) / 2


		sibs1[,i] <- rowSums(abs(sib1) - 1)
		sibs2[,i] <- rowSums(abs(sib2) - 1)
		dads[,i] <- dad1 - 1 + abs(dad2) - 1
		mums[,i] <- mum1 - 1 + abs(mum2) - 1

		# l[[i]] <- (sum(sib1[,1] == sib2[,1]) / nsnp + sum(sib1[,2] == sib2[,2]) / nsnp) / 2

	}

	# This may not be correct - getting some really large values
	ibs <- scale(sibs1) * scale(sibs2)

	# Just count how many alleles are in common
	ibs_unw <- abs(abs(sibs1 - sibs2) - 2) / 2
	return(list(dads=dads, mums=mums, sibs1=sibs1, sibs2=sibs2, ibd=ibd, ibs=ibs, ibs_unw=ibs_unw))
}



makePhen <- function(effs, indep, vy=1, vx=rep(1, length(effs)), my=0)
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	stopifnot(sum(cors^2) <= 1)
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy)) + my
	return(y)
}

chooseEffects <- function(nsnp, totvar, sqrt=TRUE)
{
	eff <- rnorm(nsnp)
	aeff <- abs(eff)
	sc <- sum(aeff) / totvar
	out <- eff / sc
	if(sqrt)
	{
		out <- sqrt(abs(out)) * sign(out)
	}
	return(out)
}

make_phenotypes <- function(fam, eff_gx, eff_xy, vx, vy, mx, my)
{
	lapply(fam, function(g)
	{
		u <- rnorm(nrow(g))
		x <- makePhen(c(eff_gx), cbind(g), vy=vx, my=mx)
		y <- makePhen(c(eff_xy), cbind(x), vy=vy, my=my)
		return(data.frame(x=x, y=y))
	})
}

join_populations <- function(l)
{
	dads <- do.call(rbind, lapply(l, function(x) x$dads))
	mums <- do.call(rbind, lapply(l, function(x) x$mums))
	sibs1 <- do.call(rbind, lapply(l, function(x) x$sibs1))
	sibs2 <- do.call(rbind, lapply(l, function(x) x$sibs2))
	ibd <- do.call(rbind, lapply(l, function(x) x$ibd))
	ibs <- do.call(rbind, lapply(l, function(x) x$ibs))
	return(list(dads=dads, mums=mums, sibs1=sibs1, sibs2=sibs2, ibd=ibd, ibs=ibs))
}

sample_populations <- function(l, n)
{
	x <- nrow(l$dads)
	index <- sort(sample(1:x, n, replace=FALSE))
	l$dads <- l$dads[index,]
	l$mums <- l$mums[index,]
	l$sibs1 <- l$sibs1[index,]
	l$sibs2 <- l$sibs2[index,]
	l$ibd <- l$ibd[index,]
	l$ibs <- l$ibs[index,]
	l$ibs_unw <- l$ibs_unw[index,]
	return(l)
}



# Dynastic effects
dynastic_phen <- function(fam, eff_gx, eff_xy, eff_ux, eff_uy, eff_xu)
{
	n <- nrow(fam$sibs1)
	
	# parents x

	umums <- rnorm(n)
	udads <- rnorm(n)
	xmums <- makePhen(c(eff_gx, eff_ux), cbind(fam$mums, umums))
	xdads <- makePhen(c(eff_gx, eff_ux), cbind(fam$dads, udads))

	u <- makePhen(c(eff_xu, eff_xu), cbind(xmums, xdads))
	x1 <- makePhen(c(eff_gx, eff_ux), cbind(fam$sibs1, u))
	x2 <- makePhen(c(eff_gx, eff_ux), cbind(fam$sibs2, u))
	y1 <- makePhen(c(eff_xy, eff_uy), cbind(x1, u))
	y2 <- makePhen(c(eff_xy, eff_uy), cbind(x2, u))
	l <- list()
	l$dads <- data.frame(x=xdads)
	l$mums <- data.frame(x=xmums)
	l$sibs1 <- data.frame(x=x1, y=y1)
	l$sibs2 <- data.frame(x=x2, y=y2)
	return(l)
}


popstrat_phen <- function(fam1, fam2, eff_gx, eff_xy, eff_ux, eff_uy)
{
	index1 <- rep(1, nrow(fam1$dads))
	index2 <- rep(2, nrow(fam2$dads))
	u <- scale(c(index1, index2))
	fam <- join_populations(list(fam1, fam2))
	xdads <- makePhen(c(eff_gx, eff_ux), cbind(fam$dads, u))
	xmums <- makePhen(c(eff_gx, eff_ux), cbind(fam$mums, u))
	xsibs1 <- makePhen(c(eff_gx, eff_ux), cbind(fam$sibs1, u))
	xsibs2 <- makePhen(c(eff_gx, eff_ux), cbind(fam$sibs2, u))
	ydads <- makePhen(c(eff_xy, eff_uy), cbind(xdads, u))
	ymums <- makePhen(c(eff_xy, eff_uy), cbind(xmums, u))
	ysibs1 <- makePhen(c(eff_xy, eff_uy), cbind(xsibs1, u))
	ysibs2 <- makePhen(c(eff_xy, eff_uy), cbind(xsibs2, u))
	l <- list()
	l$dads <- data.frame(x=xdads, y=ydads)
	l$mums <- data.frame(x=xmums, y=ymums)
	l$sibs1 <- data.frame(x=xsibs1, y=ysibs1)
	l$sibs2 <- data.frame(x=xsibs2, y=ysibs2)
	return(l)
}



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
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p
	))
}

gwas <- function(y, g)
{
	out <- matrix(0, ncol(g), 5)
	for(i in 1:ncol(g))
	{
		o <- fastAssoc(y, g[,i])
		out[i, ] <- unlist(o)
	}
	out <- as.data.frame(out)
	names(out) <- names(o)
	return(out)
}

do_mr_standard <- function(x, y, g)
{
	gwasx <- gwas(x, g)
	gwasy <- gwas(y, g)
	out <- mr_ivw(gwasx$bhat, gwasy$bhat, gwasx$se, gwasy$se)
	return(out)
}

vh_regression <- function(bS, bD, seS, seD, n)
{
	vS <- seS^2
	vD <- seD^2
	w <- (1 / vD) / (1 / vD + 1 / vS)
	bDS <- 0.5 * ((1 - w) * bS - w * bD)
	vDS <- 0.25 * w * vD
	tDS <- bDS / sqrt(vDS)
	pDS <- pt(abs(tDS), df=n-1, lower.tail=FALSE) * 2
	d <- data.frame(bDS = bDS, vDS = vDS, tDS = tDS, pDS = pDS)
	return(d)
}

do_mr_wf <- function(x1, x2, y1, y2, ibd)
{
	dmodx <- gwas((x1-x2)^2, ibd)
	dmody <- gwas((y1-y2)^2, ibd)
	smodx <- gwas((x1+x2)^2, ibd)
	smody <- gwas((y1+y2)^2, ibd)
	vhx <- vh_regression(smodx$bhat, dmodx$bhat, smodx$se, dmodx$se, length(x1))
	vhy <- vh_regression(smody$bhat, dmody$bhat, smody$se, dmody$se, length(y1))
	return(mr_ivw(vhx[,1], vhy[,1], sqrt(vhx[,2]), sqrt(vhy[,2])))
}


do_mr_pop_wf <- function(fam, phen)
{
	bx <- gwas(phen$dads$x, fam$dads)$bhat
	nsnp <- ncol(fam$dads)
	sbx <- rep(NA, nsnp)
	sby <- rep(NA, nsnp)
	ssex <- rep(NA, nsnp)
	ssey <- rep(NA, nsnp)
	for(i in 1:nsnp)
	{
		sdiffgx <- fam$sibs1[,i] * bx[i] - fam$sibs2[,i] * bx[i]
		sdiffx <- phen$sibs1$x - phen$sibs2$x
		# mod <- summary(lm(sdiffx ~ sdiffgx))
		mod <- fastAssoc(sdiffx, sdiffgx)
		sbx[i] <- mod$bhat
		ssex[i] <- mod$se
		sdiffy <- phen$sibs1$y - phen$sibs2$y
		mod <- fastAssoc(sdiffy, sdiffgx)
		sby[i] <- mod$bhat
		ssey[i] <- mod$se
	}
	return(mr_ivw(sbx, sby, ssex, ssey))
}


do_mr_trio <- function(fam, phen)
{
	nsnp <- ncol(fam$dads)
	bx <- rep(NA, nsnp)
	by <- rep(NA, nsnp)
	sex <- rep(NA, nsnp)
	sey <- rep(NA, nsnp)
	for(i in 1:nsnp)
	{
		sdiffgx <- fam$sibs1[,i] * bx[i] - fam$sibs2[,i] * bx[i]
		sdiffx <- phen$sibs1$x - phen$sibs2$x
		# mod <- summary(lm(sdiffx ~ sdiffgx))

		mod <- summary(lm(phen$sibs1$x ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i]))$coefficients
		bx[i] <- mod[2,1]
		sex[i] <- mod[2,2]
		mod <- summary(lm(phen$sibs1$y ~ fam$sibs1[,i] + fam$dads[,i] + fam$mums[,i]))$coefficients
		by[i] <- mod[2,1]
		sey[i] <- mod[2,2]
	}
	return(mr_ivw(bx, by, sex, sey))
}
