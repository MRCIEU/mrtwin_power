library(dplyr)
library(ggplot2)
library(mvtnorm)


# n <- 1000000

# g <- rbinom(n, 2, 0.5)
# pgs <- rnorm(n, sd=sqrt(0.2))

# # rsq due to variance effect of SNP is 0.5%
# b <- 0.005/2
# e1 <- rnorm(n, mean=0, sd=sqrt(0.8 + b * g))
# e2 <- rnorm(n, mean=0, sd=sqrt(0.8 + b * g))
# tapply(e1, g, var)

# c1 <- rmvnorm(n, mean=c(0,0), sigma=matrix(c((1-h2), ((1-h2)*r), ((1-h2)*r), ((1-h2))),2,2))


# y1 <- pgs + e1
# y2 <- pgs + e2

# ydif <- (y1 - y2)^2
# summary(lm(ydif ~ g))



# To estimate statistical power, we need to simulate the phenotype as a function of SNPs and random noise. Then we see how often a significant association at some alpha level is found using the model for a specific set of parameters.
# The models:

# Model 1: The SNP influences the variance
# y_i1 = g_i + e_i1 + c_i1
# y_i2 = g_i + e_i2 + c_i2
# g_i ~ N(0, h2)
# e_i1 ~ e_i2 ~ N(0, b*x_i)
# x_i ~ Binom(2, maf)
# c_i ~ MVN(0, sigma)
# sigma = [1-h2, C
#             C, 1-h2]
#
# Basically what this means is that for a MZ pair, each individual has the same genetic risk score (g). This has variance of h2. They also have a common environmental effect and a specific environmental effect. This is shown by sigma - the variance of the environmental effect is 1 - h2, but MZ1 and MZ2 have correlated environments, denoted by C (i.e. a mixture of specific and common). Finally, the SNP influences an extra environmental component. 


arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
chunks <- as.numeric(arguments[2])
out <- arguments[3]

param <- expand.grid(
	nsim = 1:100,
	n = c(10000, 20000, 30000, 40000, 50000),
	eff = seq(0.005, 0.07, by=0.005),
	maf = c(0.05, 0.15, 0.5),
	C = c(0.25, 0.5, 0.75),
	h2 = 0.2,
	model = 1:3
)

chunksize <- ceiling(nrow(param) / chunks)
t1 <- (jid - 1) * chunksize + 1
t2 <- min(jid * chunksize, nrow(param))

message("total size: ", nrow(param))
message("running: ", t1, " to ", t2)

param <- param[t1:t2, ]

for(i in 1:nrow(param))
{
	message(i, " of ", nrow(param))
	sig <- matrix(c(
		1 - param$h2[i], (1 - param$h2[i]) * param$C[i], (1 - param$h2[i]) * param$C[i], 1 - param$h2[i]), 2, 2
	)
	C <- rmvnorm(param$n[i], mean=c(0,0), sigma=sig)
	if(param$model[i] == 1)
	{
		g <- rbinom(param$n[i], 2, param$maf[i])
		pgs <- rnorm(param$n[i], sd=sqrt(0.2))
		e1 <- rnorm(param$n[i], mean=0, sd=sqrt(param$eff[i] * g))
		e2 <- rnorm(param$n[i], mean=0, sd=sqrt(param$eff[i] * g))
		y1 <- pgs + e1 + C[,1]
		y2 <- pgs + e2 + C[,2]
	}
	if(param$model[i] == 2)
	{
		g <- rbinom(param$n[i], 2, param$maf[i])
		pgs <- rnorm(param$n[i], sd=sqrt(0.2))
		# e1 <- rnorm(param$n[i], mean=0, sd=sqrt(0.8))
		# e2 <- rnorm(param$n[i], mean=0, sd=sqrt(0.8))
		y1 <- pgs + scale(g) * sqrt(param$eff[i]) + C[,1]
		y2 <- pgs + scale(g) * sqrt(param$eff[i]) + C[,2]
	}
	if(param$model[i] == 3)
	{
		g <- rbinom(param$n[i], 2, param$maf[i])
		pgs <- rnorm(param$n[i], sd=sqrt(0.2))
		e1 <- rnorm(param$n[i], mean=0, sd=sqrt(param$eff[i] * g))
		e2 <- rnorm(param$n[i], mean=0, sd=sqrt(param$eff[i] * g))
		y1 <- pgs + scale(g) * sqrt(param$eff[i]) + e1 + C[,1]
		y2 <- pgs + scale(g) * sqrt(param$eff[i]) + e2 + C[,2]
	}
	ydif <- (y1 - y2)^2
	vdiff <- (tapply(y1, g, var) + tapply(y2, g, var)) / 2
	mod <- summary(lm(ydif ~ g))
	param$g1[i] <- vdiff[1]
	param$g2[i] <- vdiff[2]
	param$g3[i] <- vdiff[3]
	param$pval[i] <- mod$coefficients[2,4]
	param$b[i] <- mod$coefficients[2,1]
}

save(param, file=out)

# a <- group_by(param, model, eff, maf, n, C) %>%
# summarise(pow = sum(pval < 5e-8) / n())

# ggplot(subset(a, maf==0.5), aes(x=eff, y=pow)) +
# geom_point(aes(colour=as.factor(n))) +
# geom_line(aes(colour=as.factor(n))) +
# facet_grid(C ~ model)
