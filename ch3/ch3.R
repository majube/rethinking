#3.1
Pr_Positive_Vampire <- 0.95
Pr_Positive_Mortal <- 0.01
Pr_Vampire <- 0.001
Pr_Positive <- Pr_Positive_Vampire * Pr_Vampire + Pr_Positive_Mortal * (1- Pr_Vampire)
(Pr_Vampire_Positive <- Pr_Positive_Vampire * Pr_Vampire / Pr_Positive)
#3.2
p_grid <- seq(from=0, to=1, length.out=1000)
prob_p <- rep(1, 1000)
prob_data <- dbinom(6, size=9, prob=p_grid)
posterior <- prob_data * prob_p
posterior <- posterior / sum(posterior)
#3.3
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
#3.4
plot(samples)
#3.5
library(rethinking)
dens(samples)
#3.6
sum(posterior[p_grid < 0.5])
#3.7
sum(samples < 0.5) / 1e4
#3.8
sum(samples > 0.5 & samples < 0.75) / 1e4
#3.9
quantile(samples, 0.8)
#3.10
quantile(samples, c(0.1, 0.9))
#3.11 Samples from skewed posterior
p_grid <- seq(from=0, to=1, length.out=1e3)
prior <- rep(1, 1e3)
likelihood <- dbinom(3, size=3, prob=p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
samples <- sample(p_grid, size=1e4, replace=TRUE, prob=posterior)
#3.12 Percentile compatibility interval
PI(samples, prob=0.5)
#3.13 Highest Posterior Density 
HPDI(samples, prob=0.5)
#3.14 Maximum a posteriori (MAP) estimate
p_grid[which.max(posterior)]
#3.15 Approximation if you only have samples from posterior
chainmode(samples, adj=0.01)
#3.16
mean(samples)
median(samples)
#3.17 calculate loss proportional for distance guess from true value using posterior
sum(posterior*abs(0.5 - p_grid))
#3.18 calculate loss for each value of p_grid
loss <- sapply(p_grid, function(d) sum(posterior*abs(d - p_grid)))
#3.19 find value with lowest loss
p_grid[which.min(loss)]
#3.20 Generate binomial dummy data
dbinom(0:2, size=2, prob=0.7)
#3.21 one dummy observation
rbinom(1, size=2, prob=0.7)
#3.22 ten dummy observations
rbinom(10, size=2, prob=0.7)
#3.23 1e5 dummy observations
dummy_w <- rbinom(1e5, size=2, prob=0.7)
table(dummy_w)/1e5 #close to 3.20
#3.24
dummy_w <- rbinom(1e5, size=9, prob=0.7)
simplehist(dummy_w, xlab="dummy water count")
#3.25
w <- rbinom(1e4, size=9, prob=0.6)
simplehist(w)
#3.26 implied observations averaged over the posterior
w <- rbinom(1e4, size=9, prob=samples)
simplehist(w)
