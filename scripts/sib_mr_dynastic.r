source("sib_mr_functions.r")


arguments <- commandArgs(T)
paramfile <- arguments[1]
jid <- as.numeric(arguments[2])
chunks <- as.numeric(arguments[3])
out <- arguments[4]

param <- read.table(paramfile, stringsAsFactors=FALSE)

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
	right <- try(rbind(
		as.data.frame(do_mr_standard(p$sibs1$x, p$sibs1$y, d$sibs1)),
		as.data.frame(do_mr_standard((p$sibs1$x - p$sibs2$x)^2, (p$sibs1$y - p$sibs2$y)^2, d[[param$gen[i]]])),
		as.data.frame(do_mr_standard((p$sibs1$x + p$sibs2$x)^2, (p$sibs1$y + p$sibs2$y)^2, d[[param$gen[i]]])),
		as.data.frame(do_mr_wf(d, p, param$gen[i])),
		as.data.frame(do_mr_pop_wf(d, p)),
		as.data.frame(do_mr_trio(d, p))
	))
	if(class(right) != "try-error")
	{
		right$test <- 1:nrow(right)
		right$model <- "dynastic"
		left <- param[rep(i, nrow(right)), ]
		res[[i]] <- cbind(left, right)
	}
}

res <- bind_rows(res)
save(res, file=out)

