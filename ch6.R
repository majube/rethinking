#6.1
set.seed(1914)
N <- 200 # num grant proposals
p <- 0.1 # proportion to select
# uncorrelated newsworthiness and trustworthiness
nw <- rnorm(N)
tw <- rnorm(N)
# select top 10% of combined scores
s <- nw + tw # total scores
q <- quantile(s, 1-p) # top 10% treshold
selected <- ifelse(s >= q, TRUE, FALSE)
cor(tw[selected], nw[selected])
#6.2
N <- 100
#set.seed(909)
height <- rnorm(N, 10, 2)
leg_prop <- runif(N, 0.4, 0.5)
leg_left <- leg_prop*height + rnorm(N, 0, 0.02)
leg_right <- leg_prop*height + rnorm(N, 0, 0.02)
d <- data.frame(height, leg_left, leg_right)
#6.3
library(rethinking)
m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2, 10),
    br ~ dnorm(2, 10),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.1)
#6.4
plot(precis(m6.1))
#6.5
post <- extract.samples(m6.1)
plot(bl ~ br, post, col=col.alpha(rangi2, 0.1), pch=16)
#6.6
sum_blbr <- post$bl + post$br
dens(sum_blbr, col=rangi2, lwd=2, xlab="sum of bl and br")
#6.7
m6.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2,  10),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.2)
mean(sum_blbr)
#6.8
data(milk)
d <- milk
d$K <- standardize(d$kcal.per.g)
d$F <- standardize(d$perc.fat)
d$L <- standardize(d$perc.lactose)
#6.9
# kcal.per.g regressed on perc.fat
m6.3 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bF*F,
    a ~ dnorm(0, 0.2),
    bF ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)
# kcal.per.g. regressed on perc.lactose
m6.4 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bL*L,
    a ~ dnorm(0, 0.2),
    bL ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.3)
precis(m6.4)
#m6.10
m6.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bF*F + bL*L,
    a ~ dnorm(0, 0.2),
    bF ~ dnorm(0, 0.5),
    bL ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.5)
#6.11
pairs( ~ kcal.per.g + perc.fat + perc.lactose, data=d, col=rangi2)
#6.12
rm(list=ls())
library(rethinking)
data(milk)
d <- milk
sim.coll <- function(r=0.9) {
  d$x <- rnorm(nrow(d), mean=r*d$perc.fat, sd=sqrt((1 - r^2)*var(d$perc.fat)))
  m <- lm(kcal.per.g ~ perc.fat + x, data=d)
  sqrt(diag(vcov(m)))[2] # stddev of parameter
}
rep.sim.coll <- function(r=0.9, n=100) {
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}
r.seq <- seq(from=0, to=0.99, by=0.01)
stddev <- sapply(r.seq, function(z) rep.sim.coll(r=z, n=100))
plot(stddev ~ r.seq, type="l", col=rangi2, lwd=2, xlab="correlation")
#6.13
set.seed(71)
# number of plants
N <- 100
# simulate initial heights
h0 <- rnorm(N, 10, 2)
# assign treatments and simulate fungus and growth
treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4)
h1 <- h0 + rnorm(N, 5 -3*fungus)
# compose a clean data frame
d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
precis(d)
#6.14
sim_p <- rlnorm(1e4, 0, 0.25)
precis(data.frame(sim_p))
#6.15
m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p ~ dlnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.6)
#6.16
m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    bf ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.7)
#6.17
m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bt*treatment,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.8)
#6.18
library(dagitty)
plant_dag <- dagitty("dag{H_0 -> H_1; T -> F -> H_1}")
coordinates(plant_dag) <- list(x=c(H_0=0,T=2,F=1.5,H_1=1), y=c(H_0=0,T=0,F=0,H_1=0))
drawdag(plant_dag)
#6.19
impliedConditionalIndependencies(plant_dag)
#6.20
set.seed(71)
N <- 1e3
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each=N/2)
M <- rbern(N)
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4 + 0.4*M)
h1 <- h0 + rnorm(N, 5 + 3*M)
d2 <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    bf ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d2
)
precis(m6.7)
m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bt*treatment,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d2
)
precis(m6.8)
#6.21
library(rethinking)
d <- sim_happiness(seed=1977, N_years=1e3)
precis(d)
#6.22
d2 <- d[d$age > 17, ] #only adults
d2$A <- (d2$age - 18)/(65- 18)
#6.23
d2$mid <- d2$married + 1
m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a[mid] + bA*A,
    a[mid] ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data=d2
)
precis(m6.9, depth=2)
#6.24
m6.10 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data=d2
)
precis(m6.10)
#6.25
N <- 200 # number of grandparent-parent-child triads
b_GP <- 1 # direct effect of G on P
b_GC <- 0 # direct effect of G on C
b_PC <- 1 # direct effect of P on C
b_U <- 2  # direct effect of U on P and C
#6.26
set.seed(1)
U <- 2*rbern(N, 0.5) -1
G <- rnorm(N)
P <- rnorm(N, b_GP*G + b_U*U)
C <- rnorm(N, b_PC*P + b_GC*G + b_U*U)
d <- data.frame(C=C, P=P, G=G, U=U)
#6.27
m6.11 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G,
    a ~ dnorm(0, 1),
    c(b_PC, b_GC) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.11)
#6.28
m6.12 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G + b_U*U,
    a ~ dnorm(0, 1),
    c(b_PC, b_GC, b_U) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.12) # after conditioning on U coefficients are correctly recovered
#6.29
library(dagitty)
dag_6.1 <- dagitty("dag{
  U [unobserved]
  X -> Y
  X <- U <- A -> C -> Y
  U -> B <- C
}")
adjustmentSets(dag_6.1, exposure="X", outcome="Y")
#6.30
library(dagitty)
dag_6.2 <- dagitty("dag{
  A -> D
  A -> M -> D
  A <- S -> M
  S -> W -> D
}")
adjustmentSets(dag_6.2, exposure="W", outcome="D")
#6.31
impliedConditionalIndependencies(dag_6.2)
