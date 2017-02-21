source("sib_mr_functions.r")

param <- expand.grid(
	nsim = 1:20,
	n = c(20000, 60000, 100000),
	nsnp = 90,
	vargx = 0.1,
	eff_xy = c(0, 0.005, 0.01, 0.05, 0.1),
	eff_ux = c(0, 0.1, 0.4),
	eff_xu = c(0, 0.1, 0.4)
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

# Analyses
# 1. Standard MR analysis using first sib only
# 2. MR analysis using SNP-exposure effects estimated in parents, and within family SNP-outcome effects estimated in children
# 3. MR analysis using within family effects estimated in children only
# 4. Mid-parental mean used as intrument?


# Dynastic effects
# parentg -> x <- confounder -> y
# 1. x causes y
# 2. parent's genotype influences child's confounder
# 3. x causes y and parent's genotype influences child's confounder

dat <- make_families(runif(90), 100000)

res <- list()
for(i in 1:nrow(param))
{
	message(i)
	d <- sample_populations(dat, param$n[i])
	eff <- chooseEffects(param$nsnp[i], param$vargx[i])
	p <- dynastic_phen(d, eff, 
		sqrt(param$eff_xy[i]), 
		sqrt(param$eff_ux[i]),
		sqrt(param$eff_ux[i]),
		sqrt(param$eff_xu[i])
	)
	res[[i]] <- rbind(
		do_mr_standard(p$sibs1$x, p$sibs1$y, d$sibs1),
		do_mr_standard((p$sibs1$x - p$sibs2$x)^2, (p$sibs1$y - p$sibs2$y)^2, d$ibd),
		do_mr_standard((p$sibs1$x + p$sibs2$x)^2, (p$sibs1$y + p$sibs2$y)^2, d$ibd),
		do_mr_wf(p$sibs1$x, p$sibs2$x, p$sibs1$y, p$sibs2$y, d$ibd),
		do_mr_pop_wf(d, p)
	)
	res[[i]] <- as.data.frame(res[[i]])
	res[[i]]$sim <- param$sim[i]
	res[[i]]$model <- "dynastic"
	res[[i]]$test <- 1:nrow(res[[i]])
}

res <- bind_rows(res)
param <- merge(param, res, by="sim")
save(param, file=out)


