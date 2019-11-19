library(dplyr)
fn <- list.files("../results/scratch/") %>% 
	grep("sib_mr_dynastic", ., value=TRUE) %>% 
	file.path("..", "results", "scratch", .)

l <- list()
for(i in 1:length(fn))
{
	load(fn[i])
	l[[i]] <- res
}

dat <- plyr::rbind.fill(l)
save(dat, file="../results/sib_mr_dynastic.rdata")

