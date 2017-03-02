l <- list()
for(i in 1:100)
{
	if(file.exists(paste0("../results/sib_popstrat_", i, ".rdata")))
	{
		load(paste0("../results/sib_popstrat_", i, ".rdata"))
		l[[i]] <- param
	}
	else { message(i)}
}

dat <- plyr::rbind.fill(l)
save(dat, file="../results/sib_popstrat.rdata")

