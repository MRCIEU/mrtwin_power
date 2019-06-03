source("sib_mr_functions.r")
param <- data_frame(
	n = 1000,
	vargx = 0.1,
	eff_xy = 0.1,
	eff_ux = 0.1,
	eff_xu = 0.1,
	nsnp = 10,
	gen = "ibs_unw"
)

dat <- make_families(runif(10), 10000)

d <- sample_populations(dat, param$n)
eff <- chooseEffects(param$nsnp, param$vargx)
p <- dynastic_phen(d, eff, 
	sqrt(param$eff_xy), 
	sqrt(param$eff_ux),
	sqrt(param$eff_ux),
	sqrt(param$eff_xu)
)
right <- rbind(
	as.data.frame(do_mr_standard(p$sibs1$x, p$sibs1$y, d$sibs1)),
	as.data.frame(do_mr_standard((p$sibs1$x - p$sibs2$x)^2, (p$sibs1$y - p$sibs2$y)^2, d[[param$gen]])),
	as.data.frame(do_mr_standard((p$sibs1$x + p$sibs2$x)^2, (p$sibs1$y + p$sibs2$y)^2, d[[param$gen]])),
	as.data.frame(do_mr_wf(d, p)),
	as.data.frame(do_mr_pop_wf(d, p)),
	as.data.frame(do_mr_trio(d, p)),
	as.data.frame(do_mr_fixed(d, p)),
	as.data.frame(do_mr_random(d, p))

)
right$test <- c("pop", "diff of diff", "diff of sum", "diff of diff 2", "diff external", "trio", "fixed", "random")
