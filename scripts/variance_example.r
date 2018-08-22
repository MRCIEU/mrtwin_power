n <- 1000
x <- rbinom(n, 2, 0.5)

y1 <- rnorm(n, x, 1)
y2 <- rnorm(n, 0, x)
y3 <- rnorm(n, x, x)

d <- data.frame(x=x, y=c(y1,y2,y3), what=rep(c("Mean", "Var", "Mean and Var"), each=n))

library(ggplot2)

ggplot(d, aes(x=as.factor(x),y=y)) +
geom_boxplot() +
facet_grid(. ~ what) +
labs(y="Phenotype", x="Genotype") 
ggsave("../images/variance_het.pdf", height=7, width=14)
