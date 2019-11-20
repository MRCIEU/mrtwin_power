library(dplyr)
fn <- list.files("../results/scratch/") %>% 
	grep("sib_mr_dynastic", ., value=TRUE) %>% 
	file.path("..", "results", "scratch", .)

l <- list()
message(length(fn))
for(i in 1:length(fn))
{
	message(i)
	load(fn[i])
	l[[i]] <- res
}

dat <- plyr::rbind.fill(l)
save(dat, file="../results/sib_mr_dynastic.rdata")

