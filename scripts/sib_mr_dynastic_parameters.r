main <- function(seed=1, output)
{
	set.seed(seed)
	param <- expand.grid(
		n = c(10000, 20000, 40000, 60000, 100000),
		nsnp = 90,
		gen = c("ibs_unw"),
		vargx = c(0.05, 0.1),
		eff_xy = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05),
		eff_ux = c(0, 0.1, 0.2),
		eff_xu = c(0, 0.1, 0.2),
		nsim = 1:1000
	)
	param$gen <- as.character(param$gen)
	write.table(param, file=output)
	return(nrow(param))
}
main(as.numeric(commandArgs(T)[1]), commandArgs(T)[2])
