#5.1 load data
rm(list=ls())
library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce
#standardize
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)
#5.2
sd(d$MedianAgeMarriage)
#5.3 regress on median age at marriage
m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data = d
)
#5.4 simulate priors
set.seed(10)
prior <- extract.prior(m5.1)
mu <- link(m5.1, post=prior, data=list(A=c(-2, 2)))
plot(NULL, xlim=c(-2, 2), ylim=c(-2, 2))
for (i in 1:50)
  lines(c(-2, 2), mu[i, ], col=col.alpha("black", 0.4))
#5.5 posterior predictions
A_seq <- seq(from=-3, to=3.2, length.out=30)
mu <- link(m5.1, data=list(A=A_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(D ~ A, data=d, col=rangi2)
lines(A_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)
precis(m5.1)
#5.6 regress on marriage rate
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=d
)
precis(m5.2)
#5.7
library(dagitty)
dag5.1 <- dagitty("dag{A -> D; A -> M; M -> D}")
coordinates(dag5.1) <- list(x=c(A=0, D=1, M=2), y=c(A=0, D=1, M=0))
drawdag(dag5.1)
#5.8
DMA_dag2 <- dagitty('dag{D <- A -> M}')
impliedConditionalIndependencies(DMA_dag2)
#5.9
DMA_dag1 <- dagitty('dag{D <- A -> M -> D}')
impliedConditionalIndependencies(DMA_dag1)
#
graphics.off()
#5.10 approximate posterior of multivariate regression
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)
precis(m5.3)
#5.11
plot(coeftab(m5.1, m5.2, m5.3), par=c("bA", "bM"))
#5.12 simulated data
N <- 50 # 50 states
age <- rnorm(N) # sim A
mar <- rnorm(N, -age) # sim A -> M
div <- rnorm(N, age) # sim A -> D
#5.13 predictor residual plot
m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bAM*A,
    a ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)
#5.14 calculate residuals: subtract observed marriage rate from predicted rate by above model
mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean
plot(d$D ~ mu_resid, xlab="Marriage rate residuals", ylab="Divorce rate (std)", col="blue")
abline(v=0, lty="dashed")
#linear regression of D on the marriage rate residuals
resid_M_lm <- lm(D ~ mu_resid, data=d)
abline(resid_M_lm)
# regress age of marriage on marriage rate
m5.5 <- quap(
  alist(
    A ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)
mu <- link(m5.5)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$A - mu_mean
plot(d$D ~ mu_resid, xlab="Age at marriage residuals", ylab="Divorce rate (std)", col="blue")
abline(v=0, lty="dashed")
resid_A_lm <- lm(d$D ~ mu_resid)
abline(resid_A_lm)
#5.15 posterior predictive plot: predictions vs observations
mu <- link(m5.3)
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)
D_sim <- sim(m5.3, n=1e4)
D_PI <- apply(D_sim, 2, PI)
#5.16
plot(mu_mean ~ d$D, col=rangi2, ylim=range(mu_PI), xlab="Observed divorce", ylab="Predicted divorce")
abline(a=0, b=1, lty=2)
for (i in 1:nrow(d)) lines(rep(d$D[i], 2), mu_PI[,i], col=rangi2)
#5.17
identify(x=d$D, y=mu_mean, labels=d$Loc)
#5.18 spurious association: causal predictor influences both outcome and spurious predictor
N <- 1e4
x_real <- rnorm(N) # x_real as Gaussian with mean 0 and stddev 1
x_spur <- rnorm(N, x_real) # x_spur as Gaussian with mean=x_real
y <- rnorm(N, x_real) # y as Gaussian with mean=x_real
d <- data.frame(y, x_real, x_spur)
pairs(d)
m_spur_assoc <- quap(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a + b_real*x_real + b_spur*x_spur,
    a ~ dnorm(0, 0.2),
    b_real ~ dnorm(0, 0.1),
    b_spur ~ dnorm(0, 0.1),
    sigma ~ dexp(1)
  ), data=d)
precis(m_spur_assoc)
#5.19 counterfactual plot
data(WaffleDivorce)
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)
m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 ),
    ## A -> M
    M ~ dnorm( mu_M , sigma_M ),
    mu_M <- aM + bAM*A,
    aM ~ dnorm( 0 , 0.2 ),
    bAM ~ dnorm( 0 , 0.5 ),
    sigma_M ~ dexp( 1 )
  ) , data = d )
precis(m5.3_A)
#5.20
A_seq <- seq(from=-2, to=2, length.out=30)
#5.21
# prep data
sim_dat <- data.frame(A=A_seq)
# simulate M and then D, using A_seq
s <- sim(m5.3_A, data=sim_dat, vars=c("M", "D"))
#5.22
plot(sim_dat$A, colMeans(s$D), ylim=c(-2, 2), type="l", xlab="manipulated A", ylab="counterfactual D")
shade(apply(s$D, 2, PI), sim_dat$A)
mtext("Total counterfactual effect of A on D")
# A -> M
plot(sim_dat$A, colMeans(s$M), ylim=c(-2, 2), type="l", xlab="manipulated A", ylab="counterfactual M")
shade(apply(s$M, 2, PI), sim_dat$A)
mtext("Counterfactual effect A -> M")
#5.23: increase A from 20 to 30
# new DF, standardized to mean 26.1 and stdev 1.24
sim2_dat <- data.frame(A=(c(20,30)-26.1)/1.24)
s2 <- sim(m5.3_A, data=sim2_dat, vars=c('M', 'D'))
mean(s2$D[,2] - s2$D[,1])
#5.24
sim_dat <- data.frame(M=seq(from=-2, to=2, length.out=30), A=0)
s <- sim(m5.3_A, data=sim_dat, vars='D')
plot(sim_dat$M, colMeans(s), ylim=c(-2,2), type='l', xlab='manipulated M', ylab='counterfactual D')
shade(apply(s, 2, PI), sim_dat$M)
mtext('Total counterfactual effect of M on D')
#5.25 manually simulating counterfactuals (without sim)
A_seq <- seq(from=-2, to=2, length.out=30)
#5.26
post <- extract.samples(m5.3_A) #get 1e4 posterior values for parameters
M_sim <- with(post, sapply(1:30, function(i) rnorm(1e3, aM + bAM * bAM+A_seq[i], sigma_M))) #get simulated simulated values for M
#5.27
D_sim <- with(post, sapply(1:30, function(i) rnorm(1e3, a + bA*A_seq[i], )))
plot(A_seq, colMeans(D_sim), type='l', xlab='manipulated A', ylab='counterfactual D')
mtext("Total counterfactual effect of A on D")
#5.28
rm(list=ls())
dev.off()
library(rethinking)
data(milk)
d <- milk
str(milk)
#5.29 standardize
d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(log(d$mass))
#5.30 simple bivariate regression between K and N
m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0, 1),
    bN ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d
)
#fails due to missing values
#5.31
d$neocortex.perc
#5.32
dcc <- d[complete.cases(d$K, d$N, d$M), ]
#5.33
m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0, 1),
    bN ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=dcc
)
#5.34 plot prior predictive distributions
prior <- extract.prior(m5.5_draft)
xseq <- c(-2, 2)
mu <- link(m5.5_draft, post=prior, data=list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for (i in 1:50) lines(xseq, mu[i,], col=col.alpha('black', 0.3))
#5.35 tighter prior for intercept
m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=dcc
)
prior <- extract.prior(m5.5)
mu <- link(m5.5, post=prior, data=list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for (i in 1:50) lines(xseq, mu[i,], col=col.alpha('black', 0.3))
#5.36 model summary
precis(m5.5)
#5.37 plot
xseq <- seq(from=min(dcc$N)-0.15, to=max(dcc$N)+0.15, length.out=30)
mu <- link(m5.5, data=list(N=xseq))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)
plot(K ~ N, data=dcc, xlab='neocortex percent (std)', ylab='kcal/g (std)')
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)
#5.38 bivariate model of K as a function of M
m5.6 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bM*M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=dcc
)
precis(m5.6)
xseq <- seq(from=min(dcc$M)-0.15, to=max(dcc$M)+0.15, length.out=30)
mu <- link(m5.6, data=list(M=xseq))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)
plot(K ~ M, data=dcc, xlab='log body mass (std)', ylab='kcal/g (std)')
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)
#5.39  multivariate model with K a function of both N and M
m5.7 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N + bM*M,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=dcc
)
precis(m5.7)
#5.40 comparison parameters bivariate and multivariate models
plot(coeftab(m5.5, m5.6, m5.7), pars=c('bM', 'bN'))
pairs(~K + M + N, dcc)
#5.41 K as a function of M with N=0
xseq <- seq(from=min(dcc$M)-0.15, to=max(dcc$M)+0.15, length.out=30)
mu <- link(m5.7, data=data.frame(M=xseq, N=0))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)
plot(NULL, xlim=range(dcc$M), ylim=range(dcc$K), xlab='log body mass (std)', ylab='kcal/g (std)')
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)
mtext('Counterfactual holding N=0')
#K as a function of N holding M=0
xseq <- seq(from=min(dcc$N)-0.15, to=max(dcc$N)+0.15, length.out=30)
mu <- link(m5.7, data=data.frame(N=xseq, M=0))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)
plot(NULL, xlim=range(dcc$N), ylim=range(dcc$K), xlab='neocortex pct (std)', ylab='kcal/g (std)')
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)
mtext('Counterfactual holding M=0')
#5.42 simulating a masking relationship
# M -> K <- N
# M -> N
dag5.42 <- dagitty("dag{M -> K; N -> K; M -> N}")
drawdag(dag5.42)
n <- 100
M <- rnorm(n)
N <- rnorm(n, M)
K <- rnorm(n, N - M)
d_sim <- data.frame(K=K, N=N, M=M)
#m5.5 bivariate (K ~ N) with this simulated data
m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim
)
precis(m5.5)
#m5.6 bivariate (K ~ M) with this simulated data
m5.6 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bM*M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim
)
precis(m5.6)
#m5.7 (K ~ M + N) with this simulated data
m5.7 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N + bM*M,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim
)
precis(m5.7)
plot(coeftab(m5.5, m5.6, m5.7), pars=c("bM", "bN"))
#5.43 sim data other two DAGs
# M -> K <- N
# M -> N
dag5.43_1 <- dagitty("dag{M -> K; N -> K; N -> M}")
drawdag(dag5.43_1)
n <- 100
N <- rnorm(100)
M <- rnorm(n, N)
K <- rnorm(n, N - M)
d_sim2 <- data.frame(K=K, N=N, M=M)
#m5.5 bivariate (K ~ N) with this simulated data
m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim2
)
precis(m5.5)
#m5.6 bivariate (K ~ M) with this simulated data
m5.6 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bM*M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim2
)
precis(m5.6)
#m5.7 (K ~ M + N) with this simulated data
m5.7 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N + bM*M,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim2
)
precis(m5.7)
plot(coeftab(m5.5, m5.6, m5.7), pars=c("bM", "bN"))
# M -> K <- N
# M <- U -> N
dag5.43_2 <- dagitty("dag{M -> K; N -> K; U -> M; U -> N}")
drawdag(dag5.43_2)
n <- 100
U <- rnorm(n)
M <- rnorm(n, U)
N <- rnorm(n, U)
K <- rnorm(n, N - M)
d_sim3 <- data.frame(K=K, N=N, M=M)
#m5.5 bivariate (K ~ N) with this simulated data
m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim3
)
precis(m5.5)
#m5.6 bivariate (K ~ M) with this simulated data
m5.6 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bM*M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim3
)
precis(m5.6)
#m5.7 (K ~ M + N) with this simulated data
m5.7 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N + bM*M,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d_sim3
)
precis(m5.7)
plot(coeftab(m5.5, m5.6, m5.7), pars=c("bM", "bN"))
#5.44 calculate Markov equivalence set
dag5.7 <- dagitty("dag{
  M -> K <- N
  M -> N }")
coordinates(dag5.7) <- list(x=c(M=0, K=1, N=2), y=c(M=0.5, K=1, N=0.5))
MElist <- equivalentDAGs(dag5.7)
drawdag(MElist)
#5.45
rm(list=ls())
dev.off()
library(rethinking)
data(Howell1)
d <- Howell1
str(d)
#5.46 prior predictive simulation for avg male/female height with model with indicator variable
mu_female <- rnorm(1e4, 178, 20)
mu_male <- rnorm(1e4, 178, 20) + rnorm(1e4, 0, 10)
precis(data.frame(mu_female, mu_male))
# prior for males is wider (larger sd) because it uses two parameters
#5.47 use of index variable
d$sex <- ifelse(d$male == 1, 2, 1)
str(d$sex)
#5.48
m5.8 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a[sex],
    a[sex] ~ dnorm(178, 20),
    sigma ~ dunif(0, 50)
  ), data=d)
precis(m5.8, depth=2)
#5.49 difference between male and female (called a contrast)
post <- extract.samples(m5.8)
post$diff_fm <- post$a[, 1] - post$a[, 2]
precis(post, depth=2)
#5.50 many categories
rm(list=ls())
library(rethinking)
data(milk)
d <- milk
levels(d$clade)
#5.51 set clade as integer factors
d$clade_id <- as.integer(d$clade)
#5.52
d$K <- standardize(d$kcal.per.g)
m5.9 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d)
labels <- paste("a[", 1:4, "]:", levels(d$clade), sep="")
plot(precis(m5.9, depth=2, pars="a"), labels=labels, xlab="expected kcal (std)")
#5.53 assign other random categorical variables to data to include in model
set.seed(63)
d$house <- sample(rep(1:4, each=8), size=nrow(d))
#5.54 fit model and plot posterior
m5.9 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a[clade_id] + h[house],
    a[clade_id] ~ dnorm(0, 0.5),
    h[house] ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d)
precis(m5.9, depth=2)
labels <- paste("h[", 1:4, "]:", c("Gryffindor", "Hufflepuff", "Ravenclaw", "Slytherin"), sep="")
plot(precis(m5.9, depth=2, pars="a"), labels=labels, xlab="expected kcal (std)")
