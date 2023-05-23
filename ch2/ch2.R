# Chapter 2
#2.1
ways <- c(0, 3, 8, 9, 0)
ways / sum(ways)

# 2.2
dbinom(6, size=9, prob=0.5)

# 2.3
# define grid
p_grid <- seq(from=0, to=1, length.out=20)
# define prior
prior <- rep(1, 20) #uniform prior
# compute likelihood at each value in grid
likelihood <- dbinom(6, size=9, prob=p_grid)
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
#standardize posterior
posterior <- unstd.posterior / sum(unstd.posterior)
#2.4
plot(p_grid, posterior, type="b", xlab="probability of water", ylab="posterior probability")
mtext("20 points, uniform prior")
#2.5
# uniform prior from 0.5 to 1 (more water than land)
prior <- ifelse(p_grid < 0.5, 0, 1)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)
plot(p_grid, posterior, type="b", xlab="probability of water", ylab="posterior probability")
mtext("20 points, uniform prior for 0.5 > p > 1")

# prior that exponentially decays either side of 0.5
prior <- exp(-5*abs(p_grid-0.5))
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)
plot(p_grid, posterior, type="b", xlab="probability of water", ylab="posterior probability")
mtext("20 points, exponential prior centered on 0.5")

#2.6
library(rethinking)
globe.qa <- quap(
  alist(
    W ~ dbinom(W+L, p) #binomial likelihood
    p ~ dunif (0, 1) #uniform prior
  ),
  data=list(W=6, L=3))
precis(globe.qa)

#2.7
W <- 6
L <- 3
curve(dbeta(x, W+1, L+1), from=0, to=1)
#quadratic approximation
curve(dnorm(x, 0.67, 0.16, lty=2, add=TRUE))

#2.8
#Metropolis algorithm
n_samples <- 1000
p <- rep(NA, n_samples)
p[1] <- 0.5
W <- 6
L <- 3
for (i in 2:n_samples) {
  p_new <- rnorm(1, p[i-1], 0.1) #get new p as random normal value; previous p value mean, sd=0.1
  if (p_new < 0) p_new <- abs(p_new) # take abs value if new p is < 0
  if (p_new > 1) p_new <- 2 - p_new # take 2 - p_new as p_new if p_new > 1
  q0 <- dbinom(W, W+L, p[i-1]) # probability with old p
  q1 <- dbinom(W, W+L, p_new) # probability with new p
  #
  p[i] <- ifelse(
    runif(1) < q1/q0,
    p_new,
    p[i-1]
    )
}

#2.9
dens(p, xlim=c(0,1))
curve(dbeta(x, W+1, L+1, lty=2, add=TRUE))