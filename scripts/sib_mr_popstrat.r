
# Population stratification
# 1. x causes y in A and B, trait distributions the same
# 2. A and B trait distributions different
# 3. x causes y in A and B and trait distributions different

source("sib_mr_functions.r")

param <- expand.grid(
	nsim = 1:20,
	n = c(20000, 60000, 100000),
	nsnp = 90,
	vargx = 0.1,
	eff_xy = c(0, 0.005, 0.01, 0.05, 0.1),
	eff_ux = c(0, 0.1, 0.4)
)
param$sim <- 1:nrow(param)

arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
chunks <- as.numeric(arguments[2])
out <- arguments[3]

chunksize <- ceiling(nrow(param) / chunks)
t1 <- (jid - 1) * chunksize + 1
t2 <- min(jid * chunksize, nrow(param))

message("total size: ", nrow(param))
message("running: ", t1, " to ", t2)

param <- param[t1:t2, ]



dat1 <- make_families(runif(90), 50000)
dat2 <- make_families(runif(90), 50000)

res <- list()
for(i in 1:3)
{
	message(i)
	n1 <- n2 <- param$n[i] / 2
	d1 <- sample_populations(dat1, n1)
	d2 <- sample_populations(dat2, n2)
	d <- join_populations(list(d1, d2))

	eff <- chooseEffects(param$nsnp[i], param$vargx[i])
	p <- popstrat_phen(d1, d2, eff, 
		sqrt(param$eff_xy[i]), 
		sqrt(param$eff_ux[i]),
		sqrt(param$eff_ux[i])
	)
	res[[i]] <- rbind(
		do_mr_standard(p$sibs1$x, p$sibs1$y, d$sibs1),
		do_mr_standard((p$sibs1$x - p$sibs2$x)^2, (p$sibs1$y - p$sibs2$y)^2, d$ibd),
		do_mr_standard((p$sibs1$x + p$sibs2$x)^2, (p$sibs1$y + p$sibs2$y)^2, d$ibd),
		do_mr_wf(p$sibs1$x, p$sibs2$x, p$sibs1$y, p$sibs2$y, d$ibd),
		do_mr_pop_wf(d, p)
	)
	res[[i]] <- as.data.frame(res[[i]])
	res[[i]]$sim <- i
	res[[i]]$model <- "popstrat"
	res[[i]]$test <- 1:nrow(res[[i]])
}

res <- bind_rows(res)
param <- merge(param, res, by="sim")
save(param, file=out)
