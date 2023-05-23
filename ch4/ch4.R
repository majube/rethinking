#4.1 normal through addition
pos <- replicate(1e4, sum(runif(16, min=-1, max=1)))
hist(pos)
plot(density(pos)) #approximately normal
#4.2 multiply 1 + 12 random values from unif(0, 0.1)
prod(1 + runif(12, min=0, max=0.1))
#4.3 normal through multiplication
growth <- replicate(1e4, prod(1 + runif(12, min=0, max=0.1)))
hist(growth)
plot(density(growth))
library(rethinking)
dens(growth, norm.comp=TRUE) #approximately normal
#4.4
big <- replicate(1e4, prod(1 + runif(12, min=0, max=0.5)))
dens(big, norm.comp=TRUE) #bigger multiplications not normal
small <- replicate(1e4, prod(1 + runif(12, min=0, max=0.1)))
dens(small, norm.comp=TRUE) #smaller are normal
#4.5 log-multiplication is normal
log.big <- replicate(1e4, log(prod(1 + runif(12, min=0, max=0.5))))
dens(log.big, norm.comp=TRUE)
#4.6 using model distribution to define the posterior
w <- 6; n <- 9;
p_grid <- seq(from=0, to=1, length.out=1e2)
posterior <- dbinom(w, n, p_grid)*dunif(p_grid,0,1) #likelihood*prior
posterior <- posterior / sum(posterior) #normalize posterior
plot(posterior)
#4.7
data(Howell1)
d <- Howell1
#4.8
str(d)
#4.9
precis(d)
#4.10 access to column of df
d$height
#4.11 only adults
d2 <- d[d$age >= 18, ]
#4.12 check prior for height
curve(dnorm(x, 178, 20), from=100, to=250)
#4.13 uniform prior for sigma
curve(dunif(x, 0, 50), from=-10, to=60)
#4.14 prior predictive simulation: sample prior
sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)
#4.15 with wider prior
sample_mu <- rnorm(1e4,178, 100)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)
sum(prior_h < 0)/1e4 #4% negative heights
#4.16 grid approximation posterior
mu.list <- seq(from=150, to=160, length.out=100)
sigma.list <- seq(from=7, to=9, length.out=100)
post <- expand.grid(mu=mu.list, sigma=sigma.list) #create grid
post$LL <- sapply(1:nrow(post), function(i) sum(dnorm(d2$height, post$mu[i], post$sigma[i], log=TRUE)))
post$prod <- post$LL + dnorm(post$mu, 178, 20, TRUE) + dunif(post$sigma, 0, 50, TRUE)
post$prob <- exp(post$prod - max(post$prod))
#4.17 contour plot posterior
contour_xyz(post$mu, post$sigma, post$prob, xlab='mu', ylab='sigma')
#4.18 heatmap posterior
image_xyz(post$mu, post$sigma, post$prob)
#4.19 sample posterior
sample.rows <- sample(1:nrow(post), size=1e4, replace=TRUE, prob=post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]
#4.20 
plot(sample.mu, sample.sigma, cex=0.5, pch=16, col=col.alpha(rangi2, 0.1))
#4.21 marginal posterior densities
dens(sample.mu, adj=1)
dens(sample.sigma)
#4.22
PI(sample.mu)
PI(sample.sigma)
#4.23 sample of 20 heights
d3 <- sample(d2$height, size=20)
#4.24 repeat of previous with 20 heights
mu.list <- seq(from=150, to=170, length.out=200)
sigma.list <- seq(from=4, to=20, length.out=200)
post2 <- expand.grid(mu=mu.list, sigma=sigma.list) #create grid
post2$LL <- sapply(1:nrow(post2), function(i) sum(dnorm(d3, post2$mu[i], post2$sigma[i], log=TRUE)))
post2$prod <- post2$LL + dnorm(post2$mu, 178, 20, TRUE) + dunif(post2$sigma, 0, 50, TRUE)
post2$prob <- exp(post2$prod - max(post2$prod))
sample2.rows <- sample(1:nrow(post2), size=1e4, replace=TRUE, prob=post2$prob)
sample2.mu <- post2$mu[sample2.rows]
sample2.sigma <- post2$sigma[sample2.rows]
plot(sample2.mu, sample2.sigma, cex=0.5, col=col.alpha(rangi2, 0.1), xlab="mu", ylab="sigma", pch=16) #much longer tail of sigma
#4.25 marginal density has fatter tail
dens(sample2.sigma, norm.comp=TRUE)
#4.26 clear env, reload data
rm(list=ls())
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[d$age >= 18, ]
#4.27 define model
flist <- alist( 
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0, 50)
)
#4.28 fit model
m4.1 <- quap(flist, data=d2)
#4.29
precis(m4.1)
#4.30 specify start values for quap
start <- list(
  mu=mean(d2$height),
  sigma=sd(d2$height)
)
m4.1 <- quap(flist, data=d2, start=start)
precis(m4.1)
#4.31 narrower prior for mu
m4.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(178, 0.1),
    sigma ~ dunif(0, 50)
  ),
  data=d2,
)
#sigma much larger because mu is constrained by prior
precis(m4.2)
#4.32 var-cov matrix of fitted model
vcov(m4.2)
#4.33 elements: variances and correlations
diag(vcov(m4.1))
cov2cor(vcov(m4.1))
#4.34 convenience function in rethinking for extracting samples from posterior
post <- extract.samples(m4.1, n=1e4)
head(post)
Sys.setlocale(locale='Chinese') #fixes histogram output
precis(post)
#4.36 sampling using mvrnorm (multivariate randon normal)
library(MASS)
post <- mvrnorm(n=1e4, mu=coef(m4.1), Sigma=vcov(m4.1))
#4.37 scatterplot height ~ weight
plot(d2$height ~ d2$weight)
#4.38 prior predictive simulation for linear regression
set.seed(2971)
N <- 100
a <- rnorm(N, 178, 20)
b <- rnorm(N, 0, 10)
#4.39 plot 100 random lines from the prior
plot(NULL, xlim=range(d2$weight), ylim=c(-100, 400), xlab="weight", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1, lwd=0.5)
mtext("b ~ dnorm(0, 10)")
xbar <- mean(d2$weight)
for (i in 1:N) curve(a[i] + b[i]*(x - xbar), 
                     from=min(d2$weight), to=max(d2$weight), add=TRUE, col=col.alpha("black", 0.2))
#4.40 constrain slope to > 0: use log-normal
b <- rlnorm(1e4, 0, 1)
dens(b, xlim=c(0,5), adj=0.1)
#4.41 plotting as before
a <- rnorm(N, 178, 20)
b <- rlnorm(N, 0, 1)
plot(NULL, xlim=range(d2$weight), ylim=c(-100, 400), xlab="weight", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1, lwd=0.5)
mtext("b ~ dlnorm(0, 1)")
for (i in 1:N) curve(a[i] + b[i]*(x - xbar), 
                     from=min(d2$weight), to=max(d2$weight), add=TRUE, col=col.alpha("black", 0.2))
#4.42 fit model again
rm(list=ls())
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[d$age >= 18, ]
xbar <- mean(d2$weight)
m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ),
  data=d2
)
precis(m4.3)
#4.43 avoiding the log-normal distribution
m4.3b <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + exp(log_b)*(weight - xbar),
    a ~ dnorm(178, 20),
    log_b ~ dnorm(0, 1),
    sigma ~ dunif(0, 20)
  ),
  data=d2
)
precis(m4.3b)
exp(coef(m4.3b)[2]) # same as before
#4.44  estimates of marginal posterior distrubitions of parameters 
precis(m4.3)
#4.45 var-cov matrix
round(vcov(m4.3), 3)
#4.46 plotting posterior inference against data
plot(height ~weight, data=d2, col=rangi2)
post <- extract.samples(m4.3)
a_map <- mean(post$a) #maximum a posteriori
b_map <- mean(post$b)
curve(a_map + b_map*(x- xbar), add=TRUE)
#4.47 get samples from posterior
post <- extract.samples(m4.3)
post[1:5, ]
#4.48 fit model to first 10 heights
N <- 10
dN <- d2[1:N, ]
mN <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - mean(weight)),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ),
  data=dN
)
#4.49 plot 20 lines from posterior
post <- extract.samples(mN, n=20)
plot(dN$weight, dN$height, xlim=range(d2$weight), ylim=range(d2$height), col=rangi2, xlab="weight", ylab="height")
mtext(concat("N = ", N))
for (i in 1:20)
  curve(post$a[i] + post$b[i]*(x-mean(dN$weight)), col=col.alpha("black", 0.3), add=TRUE)
#4.50 find value for specific value
post <- extract.samples(m4.3)
mu_at_50 <- post$a + post$b*(50 - xbar)
#4.51
dens(mu_at_50, col=rangi2, lwd=2, xlab="mu|weight=50", adj=1)
#4.52
PI(mu_at_50, prob=0.89)
#4.53 link computes the value of each linear model at each sample for each case in the data
mu <- link(m4.3)
str(mu)
#4.54
weight.seq <- seq(from=25, to=70, by=1)
mu <- link(m4.3, data=data.frame(weight=weight.seq))
str(mu)
#4.55 posterior distribution of height for each weight value
plot(height ~ weight, d2, type="n")
for (i in 1:100)
  points(weight.seq, mu[i,], pch=16, col=col.alpha(rangi2, 0.1))
#4.56 summarize distribution of mu
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
#4.57 plot raw data, MAP line and 89% PI
plot(height ~ weight, data=d2, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq) #doesn't work
#4.58 manually without link
post <- extract.samples(m4.3) #get posterior parameter samples
mu.link <- function(weight) post$a + post$b*(weight - xbar) #linear function
weight.seq <- seq(from=25, to=70, by=1) #predictor variable grid
mu <- sapply(weight.seq, mu.link) #apply link function to predictor variable
mu.mean <- apply(mu, 2, mean) #calculate mean
mu.CI <- apply(mu, 2, PI, prob=0.89) #calculate compatibility interval
#4.59 simulate heights based on weights
sim.height <- sim(m4.3, data=list(weight=weight.seq))
str(sim.height)
#4.60
height.PI <- apply(sim.height, 2, PI, prob=0.89)
#4.61 plot raw data, MAP, line and 89% PI
plot(height ~ weight, d2, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq) #doesn't work
shade(height.PI, weight.seq) #doesn't work
#4.62 smoother PI for simulated heights (large n)
sim.height <- sim(m4.3, data=list(weight=weight.seq), n=1e4)
height.PI <- apply(sim.height, 2, PI, prob=0.89)
#4.63 manually creating simulated values
post <- extract.samples(m4.3)
weight.seq <- 25:70
sim.height <- sapply(weight.seq, function(weight)
  rnorm(
    n=nrow(post),
    mean=post$a + post$b*(weight - xbar),
    sd=post$sigma
  )
)
#4.64 reload data
rm(list=ls())
library(rethinking)
data(Howell1)
d <- Howell1
#4.65 standardize variables
d$weight_s <- (d$weight - mean(d$weight))/sd(d$weight)
d$weight_s2 <- d$weight_s^2
m4.5 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0, 1),
    b2 ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ),
  data=d
)
#4.66
precis(m4.5)
#4.67 calculate mean relationship and intervals
weight.seq <- seq(from=-2.2, to=2, length.out=30)
pred_dat <- list(weight_s=weight.seq, weight_s2=weight.seq^2)
mu <- link(m4.5,data=pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.height <- sim(m4.5, data=pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)
#4.68 plot
plot(height ~ weight_s, d, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)
#4.69 cubic regression
d$weight_s3 <- d$weight_s^3
m4.6 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2 + b2*weight_s3,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0, 1),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ), 
  data=d
)
#4.70 converting back to natural scale (from standardized)
plot(height ~ weight_s, d, col=col.alpha(rangi2, 0.5), xaxt="n") #no x-axis
#4.71
at <- c(-2, -1, 0, 1, 2)
labels <- at*sd(d$weight) + mean(d$weight)
axis(side=1, at=at, labels=round(labels, 1))
#4.72
rm(list=ls())
library(rethinking)
data(cherry_blossoms)
d <- cherry_blossoms
precis(d)
#4.73 create knots for spline
d2 <- d[complete.cases(d$doy), ] #only complete cases
num_knots <- 15
knot_list <- quantile(d2$year, probs=seq(0, 1, length.out=num_knots))
#4.74 create basis functions
library(splines)
B <- bs(d2$year, knots=knot_list[-c(1, num_knots)], degree=3, intercept=TRUE)
#4.75 plot basis functions
plot(NULL, xlim=range(d2$year), ylim=c(0, 1), xlab="year", ylab="basis")
for (i in 1:ncol(B))
  lines(d2$year, B[, i])
#4.76 fit model
m4.7 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w, # %*% is matrix multiplication
    a ~ dnorm(100, 10),
    w ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ),
  data=list(D=d2$doy, B=B),
  start=list(w=rep(0, ncol(B)))
)
#4.77 plot weighted basis functions
precis(m4.7, depth=2)
post <- extract.samples(m4.7)
w <- apply(post$w, 2, mean)
plot(NULL,xlim=range(d2$year), ylim=c(-6, 6), xlab="year", ylab="basis*weight")
for (i in 1:ncol(B)) lines(d2$year, w[i]*B[, i])
#4.78 plot data and 97% PI for mu
mu <- link(m4.7)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, d2$year, col=col.alpha("black", 0.5))
#4.79 alternative model specification without matrix multiplication
m4.7alt <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + sapply(1:827, function(i) sum(B[i, ]*w)),
    a ~ dnorm(100, 1),
    w ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ),
  data=list(D=d2$doy, B=B),
  start=list(w=rep(0, ncol(B)))
)
precis(m4.7alt, depth=2)
mu <- link(m4.7alt)
mu_PI <- apply(mu, 2, PI, 0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, 0.3), pch=16)
shade(mu_PI, d2$year, col=col.alpha("black", 0.5))
