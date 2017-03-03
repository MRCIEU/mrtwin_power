

n <- 1000
q <- sqrt(0.2)
f <- sqrt(0.2)
r <- sqrt(1-q^2-f^2)
alpha <- 0.01

vard2 <- 8 * r^2 * (1 - f^2) + (7/2) * q^4
vars2 <- 8 * (1 + f^2) * (1 + q^2 + f^2) + (7/2) * q^4
(tval <- 0.5 * (n - 2) * q^4 / vard2 - 0.5 * (n - 2) * q^4 / vars2) + 0.5

pt(tval, n, lower.tail=FALSE)
thr <- qt(1-alpha, n)
1 - pt(thr, n, ncp=tval)




F <- rnorm(n, sd=f)
q1 <- rnorm(n, sd=q)
p <- rbinom(n, 2, 0.5) / 2
q2 <- rep(NA, n)
q2[p==0] <- rnorm(sum(p==0), sd=q)
q2[p==0.5] <- rnorm(sum(p==0.5), mean=q1[p==0.5] / 2, sd = sqrt(0.75 * q^2))
q2[p==1] <- q1[p==1]

Y1 <- scale(q1 + F + rnorm(n, sd=r))
Y2 <- scale(q2 + F + rnorm(n, sd=r))


s2 <- (Y1 + Y2)^2
d2 <- (Y1 - Y2)^2





source("sib_mr_functions.r")

dat <- make_families(runif(90), 100000)

fam <- sample_populations(dat, 50000)
eff <- chooseEffects(90, 0.1)
phen <- dynastic_phen(fam, eff, 
	sqrt(0), 
	sqrt(0),
	sqrt(0),
	sqrt(0)
)

do_mr_pop_wf(fam, phen)$b^2




do_mr_pop_wf2(d, p)




