library(mvtnorm)

i <- 1
param <- list()
param$h2 <- 0.2
param$n <- 1000000
param$eff <- 0.02
param$maf <- 0.5
param$C <- 0.3



sig <- matrix(c(
1 - param$h2[i], (1 - param$h2[i]) * param$C[i], (1 - param$h2[i]) * param$C[i], 1 - param$h2[i]), 2, 2
)
C <- rmvnorm(param$n[i], mean=c(0,0), sigma=sig)
g <- rbinom(param$n[i], 2, param$maf[i])
pgs <- rnorm(param$n[i], sd=sqrt(param$h2[i]))
e1 <- rnorm(param$n[i], mean=0, sd=sqrt(param$eff[i] * g))
e2 <- rnorm(param$n[i], mean=0, sd=sqrt(param$eff[i] * g))
y1 <- pgs + e1 + C[,1]
y2 <- pgs + e2 + C[,2]
ydiff1 <- abs(y1-y2)^2
summary(lm(ydiff1 ~ g))
ydiff2 <- abs(y1+y2)^2
summary(lm(ydiff2 ~ g))

cor(ydiff1, ydiff2)


var(y1)

var(e1)
var(C[,1])

vdiff <- (tapply(y1, g, var) + tapply(y2, g, var)) / 2
vdiff

tapply(e1, g, var)


var(e1)
sum(table(g) * c(0,1,2) * param$eff[i] / param$n[i])

sum(c(param$maf[i]^2, param$maf[i] * 2 * (1-param$maf[i]), (1-param$maf[i])^2) * c(2,1,0) * param$eff[i])


2 * param$maf[i]^2 * param$eff[i] + param$maf[i] * 2 * (1-param$maf[i]) * param$eff[i]

ydiff <- abs(y1-y2)

tapply(ydiff, g, mean)

summary(lm(ydiff ~ g))