l <- list()
for(i in 1:16)
{
	load(paste0("../results/mz", i, ".rdata"))
	l[[i]] <- res
}

dat <- plyr::rbind.fill(l)
save(dat, file="../results/mz.rdata")

