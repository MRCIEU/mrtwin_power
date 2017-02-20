library(ggplot2)
library(reshape2)

makePhen <- function(effs, indep, vy=1, vx=rep(1, length(effs)))
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	stopifnot(sum(cors^2) <= 1)
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy))
	return(y)
}

base_population <- function(nparents, allele_frequency, effect_size)
{

	nsnp <- length(allele_frequency)
	males1 <- matrix(FALSE, nparents, nsnp)
	males2 <- matrix(FALSE, nparents, nsnp)
	females1 <- matrix(FALSE, nparents, nsnp)
	females2 <- matrix(FALSE, nparents, nsnp)
	for(i in 1:nsnp)
	{
		males1[,i] <- rbinom(nparents, 1, allele_frequency[i])
		males2[,i] <- rbinom(nparents, 1, allele_frequency[i])
		females1[,i] <- rbinom(nparents, 1, allele_frequency[i])
		females2[,i] <- rbinom(nparents, 1, allele_frequency[i])
	}
	population <- list(males1=males1, males2=males2, females1=females1, females2=females2)
	return(population)
}


reproduce <- function(population)
{

	# Make two children per family
	nparents <- nrow(population$males1)
	nchildren <- nparents * 2
	population_size <- nchildren
	nsnp <- ncol(population$males1)

	male_partner <- 1:nparents
	female_partner <- 1:nparents

	x <- population$males1[male_partner, , drop=FALSE]
	locus <- sample(1:ncol(x), ncol(x)/2, replace=FALSE)
	x[, locus] <- population$males2[male_partner, locus, drop=FALSE]
	child_paternal1 <- x

	x <- population$males1[male_partner, , drop=FALSE]
	locus <- sample(1:ncol(x), ncol(x)/2, replace=FALSE)
	x[, locus] <- population$males2[male_partner, locus, drop=FALSE]
	child_paternal2 <- x

	x <- population$females1[female_partner, , drop=FALSE]
	locus <- sample(1:ncol(x), ncol(x)/2, replace=FALSE)
	x[, locus] <- population$females2[female_partner, locus, drop=FALSE]
	child_maternal1 <- x

	x <- population$females1[female_partner, , drop=FALSE]
	locus <- sample(1:ncol(x), ncol(x)/2, replace=FALSE)
	x[, locus] <- population$females2[female_partner, locus, drop=FALSE]
	child_maternal2 <- x

	fam <- data.frame(
		fid = rep(1:nparents, times=4),
		iid = c(1:(nparents*4)),
		pat = c(rep(0, nparents*2), rep(1:nparents, times=2)),
		mat = c(rep(0, nparents*2), rep(1:nparents+nparents, times=2)),
		sex = c(rep(1:2, each=nparents), rep(1:2, each=nparents)),
		phen = -9
	)

	paternal <- rbind(
		population$males1,
		population$females1,
		child_paternal1,
		child_paternal2
	)

	maternal <- rbind(
		population$males2,
		population$females2,
		child_maternal1,
		child_maternal2
	)

	return(list(pat = paternal, mat = maternal, id=fam))
}


make_families <- function(nparents, allele_frequency, effect_size)
{
	population <- base_population(nparents, allele_frequency, effect_size)
	return(reproduce(population))
}








dat <- make_population(20000, runif(100), NULL)


get_sibling_similarities <- function(dat)
{
	nparents <- nrow(dat$id) / 4
	# sib1 <- 1:nparents + nparents*2
	# sib2 <- 1:nparents + nparents*3

	sib1 <- 1:nparents
	sib2 <- 1:nparents + nparents


	ibd <- (dat$pat[sib1,] == dat$pat[sib2,]) + (dat$mat[sib1,] == dat$mat[sib2,])
	pihat <- rowMeans(ibd) / 2
	hist(pihat)


}




get_freqs <- function(generations, fitness_effects)
{
	require(reshape2)
	require(lubridate)
	a <- matrix(sapply(generations, function(x)
	{
		
		a <- matrix(sapply(x, function(y)
		{
			colMeans(y)
		}), ncol=4)
		rowMeans(a)
	}), ncol=length(generations))
	a <- melt(a, c("mutation", "generation"), value.name="freq")
	b <- data.frame(mutation=1:length(fitness_effects), fitness_effects=fitness_effects)
	a <- merge(a, b, by="mutation")
	a <- a[order(a$generation, a$mutation), ]
	# a$date <- Sys.Date()
	# a$date <- a$date + years(25 * (a$generation-1))
	a$year <- 2016 + (a$generation - 1) * 25
	names(a) <- c("mutation", "generation", "allele_frequency", "fitness_effects", "year")
	a$mutation <- paste0("position_", formatC(a$mutation, width=2, format="d", flag="0"))
	return(a)
}


run_generations <- function(generations, population_size, fitness_effects, number_of_generations, mutation_rate_per_generation_per_locus)
{
	l <- length(generations)
	for(i in 1:number_of_generations)
	{
		generations[[i + l]] <- reproduce(generations[[l + i - 1]], population_size, fitness_effects, mutation_rate_per_generation_per_locus)
	}
	return(generations)
}



plot_people <- function(population)
{
	dat1 <- data.frame(val=rowSums(population$males1) + rowSums(population$males2), sex="Male")
	dat1 <- dat1[order(dat1$val, decreasing=TRUE), ]
	dat1$id <- 1:nrow(dat1)
	d <- ceiling(sqrt(nrow(dat1)))
	dat1$x <- ceiling(dat1$id / d)
	dat1$y <- (dat1$id-1) %% d

	dat2 <- data.frame(val=rowSums(population$females1) + rowSums(population$females2), sex="Female")
	dat2 <- dat2[order(dat2$val, decreasing=TRUE), ]
	dat2$id <- 1:nrow(dat2)
	d <- ceiling(sqrt(nrow(dat2)))
	dat2$x <- ceiling(dat2$id / d)
	dat2$y <- (dat2$id-1) %% d

	dat <- rbind(dat1, dat2)

	p <- ggplot(dat, aes(x=x, y=y)) +
	geom_point(aes(colour=val), size=5) +
	scale_colour_gradient(low="red", high="blue") +
	facet_grid(. ~ sex) +
	theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), legend.position="none")
	return(p)
}


plot_frequencies <- function(f)
{
	p <- ggplot(f, aes(x=generation, y=freq)) +
		geom_line(aes(group=mutation, col=fitness_effects)) +
		ylim(c(0,1)) +
		scale_colour_gradientn(colours=c("blue", "black", "red"))	
	return(p)
}



library(dplyr)


make_families <- function(af, nfam)
{

	nsnp <- length(af)
	dads <- list()
	mums <- list()
	sibs1 <- list()
	sibs2 <- list()
	for(i in 1:nfam)
	{
		dad1 <- rbinom(nsnp, 1, af1)
		dad2 <- rbinom(nsnp, 1, af1)
		mum1 <- rbinom(nsnp, 1, af1)
		mum2 <- rbinom(nsnp, 1, af1)

		dadindex <- sample(c(TRUE, FALSE), nsnp, replace=TRUE)
		dadh <- rep(NA, nsnp)
		dadh[dadindex] <- dad1[dadindex]
		dadh[!dadindex] <- dad2[!dadindex]

		mumindex <- sample(c(TRUE, FALSE), nsnp, replace=TRUE)
		mumh <- rep(NA, nsnp)
		mumh[mumindex] <- mum1[mumindex]
		mumh[!mumindex] <- mum2[!mumindex]

		sib1 <- cbind(dadh, mumh)

		dadindex <- sample(c(TRUE, FALSE), nsnp, replace=TRUE)
		dadh <- rep(NA, nsnp)
		dadh[dadindex] <- dad1[dadindex]
		dadh[!dadindex] <- dad2[!dadindex]

		mumindex <- sample(c(TRUE, FALSE), nsnp, replace=TRUE)
		mumh <- rep(NA, nsnp)
		mumh[mumindex] <- mum1[mumindex]
		mumh[!mumindex] <- mum2[!mumindex]

		sib2 <- cbind(dadh, mumh)

		sibs1[[i]] <- rowSums(sib1)
		sibs2[[i]] <- rowSums(sib2)
		dads[[i]] <- dad1 + dad2
		mums[[i]] <- mum1 + mum2

		# l[[i]] <- (sum(sib1[,1] == sib2[,1]) / nsnp + sum(sib1[,2] == sib2[,2]) / nsnp) / 2

	}

	sibs1 <- do.call(rbind, sibs1)
	sibs2 <- do.call(rbind, sibs2)
	dads <- do.call(rbind, dads)
	mums <- do.call(rbind, mums)
	return(list(dads=dads, mums=mums, sibs1=sibs1, sibs2=sibs2))
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
	return(list(dads=dads, mums=mums, sibs1=sibs1, sibs2=sibs2))
}



a <- make_families(runif(90), 10000)
b <- make_families(runif(90), 10000)
eff <- chooseEffects(nsnp, 0.1)
ap <- make_phenotypes(a, eff, sqrt(0), 4, 2, 40, 20)
bp <- make_phenotypes(b, eff, sqrt(0), 4, 2, 23, 10)

ab <- join_populations(list(a,b))
abp <- join_populations(list(ap,bp))


lapply(ap, cor)
lapply(bp, cor)
lapply(abp, cor)


bx <- summary(lm(ap$dad$x ~ a$dad))$coefficients[-1,1]
by <- summary(lm(ap$dad$y ~ a$dad))$coefficients[-1,1]

plot(by ~ bx)
summary(lm(by ~ bx))


bxab <- summary(lm(abp$dad$x ~ ab$dad))$coefficients[-1,1]
by <- summary(lm(abp$dad$y ~ ab$dad))$coefficients[-1,1]
plot(by ~ bx)
summary(lm(by ~ bx))

bxab <- summary(lm(abp$mum$x ~ ab$mum))$coefficients[-1,1]
by <- summary(lm(abp$mum$y ~ ab$mum))$coefficients[-1,1]
plot(by ~ bx)
summary(lm(by ~ bx))



sbx <- rep(NA, nsnp)
sby <- rep(NA, nsnp)
ssex <- rep(NA, nsnp)
ssey <- rep(NA, nsnp)
for(i in 1:nsnp)
{
	sdiffgx <- a$sibs1[,i] * bx[i] - a$sibs2[,i] * bx[i]
	sdiffx <- ap$sibs1$x - ap$sibs2$x
	mod <- summary(lm(sdiffx ~ sdiffgx))
	sbx[i] <- mod$coef[2,1]
	ssex[i] <- mod$coef[2,2]
	sdiffy <- ap$sibs1$y - ap$sibs2$y
	mod <- summary(lm(sdiffy ~ sdiffgx))
	sby[i] <- mod$coef[2,1]
	ssey[i] <- mod$coef[2,2]
}

plot(sby ~ sbx)
summary(lm(sby ~ sbx))
library(TwoSampleMR)

mr_ivw(sbx, sby, ssex, ssey)

plot(y = sby/sbx * sqrt(1/ssex), x=sqrt(1/ssex) )
summary(lm(sby/sbx * sqrt(1/ssex) ~ sqrt(1/ssex)))
abline(lm(sby/sbx * sqrt(1/ssex) ~ sqrt(1/ssex)))






summary(glm(rep(1:2, each=10000) ~ ab$dads))

get_pop_estimate <- function(ab, abp)
{

}





# simulate a confounder
n <- 100000
g1 <- rbinom(n, 2, 0.7)
g2 <- rbinom(n, 2, 0.3)
x1 <- 27 + g1 * 0.2 + rnorm(n, sd=4)
x2 <- 27 + g2 * 0.2 + rnorm(n, sd=4)

summary(lm(x1 ~ g1))
summary(lm(x2 ~ g2))

summary(lm(c(x1, x2) ~ c(g1, g2)))


