#9.1
num_weeks <- 1e5
positions <- rep(0,num_weeks)
current <- 10
for (i in 1:num_weeks) {
  # record current position
  positions[i] <- current
  # flip coin to generate proposal
  proposal <- current + sample(c(-1,1), size=1)
  # now make sure he loops around the archipelago
  if (proposal < 1) proposal <- 10
  if (proposal > 10) proposal <- 1
  # move?
  prob_move <- proposal/current
  current <- ifelse(runif(1) < prob_move, proposal, current)
}

#9.2
plot(1:100, positions[1:100])

#9.3
plot(table(positions))

#9.4
library(rethinking)
Ds <- c(1, 10, 1e2, 1e3)
T <- 1e3
rad_dist <- function(Y) sqrt(sum(Y^2))
for (i in 1:4) {
  Y <- rmvnorm(T, rep(0,Ds[i]), diag(Ds[i]))
  Rd <- sapply(1:T, function(i) rad_dist(Y[i,]))
  dens(Rd, xlim=c(0, 35), ylim=c(0, 1), add=ifelse(i == 1, FALSE, TRUE), xlab="Radial distance from mode", ylab="Density")
}

#9.5
# U needs to return neg-log-probability
U <- function(q, a=0, b=1, k=0, d=1) {
  muy <- q[1]
  mux <- q[2]
  U <- sum(dnorm(y,muy,1,log=TRUE)) + sum(dnorm(x,mux,1,log=TRUE)) +
    dnorm(muy,a,b,log=TRUE) + dnorm(mux,k,d,log=TRUE)
  return(-U)
}

#9.6
# gradient function
# need vector of partial derivatives of U with respect to vector q
U_gradient <- function(q, a=0, b=1, k=0, d=1) {
  muy <- q[1]
  mux <- q[2]
  G1 <- sum(y - muy) + (a - muy)/b^2 #dU/dmuy
  G2 <- sum(x - mux) + (k - mux)/d^2 #dU/dmux
  return(c(-G1, -G2)) # negative bc energy is neg-log-prob
}
# test data
set.seed(7)
y <- rnorm(50)
x <- rnorm(50)
x <- as.numeric(scale(x))
y <- as.numeric(scale(y))

#9.7
library(shape) # for fancy arrows
Q <- list()
Q$q <- c(-0.1,0.2) # starting point
pr <- 0.3
plot(NULL, ylab="muy", xlab="mux", xlim=c(-pr,pr), ylim=c(-pr,pr))
step <- 0.03
L <- 11 # 0.03/28 for U-turns --- 11 for working example
n_samples <- 4
path_col <- col.alpha("black",0.5)
points(Q$q[1], Q$q[2], pch=4, col="black")
for (i in 1:n_samples) {
  Q <- HMC2(U, U_gradient, step, L, Q$q)
  if (n_samples < 10) {
    for (j in 1:L) {
      K0 <- sum(Q$ptraj[j,]^2)/2 # kinetic energy
      lines(Q$traj[j:(j+1),1], Q$traj[j:(j+1),2], col=path_col, lwd=1+2*K0)
    }
    points(Q$traj[1:L+1,], pch=16, col="white", cex=0.35)
    Arrows(Q$traj[L,1], Q$traj[L,2], Q$traj[L+1,1], Q$traj[L+1,2],
            arr.length=0.35, arr.adj = 0.7)
    text(Q$traj[L+1,1], Q$traj[L+1,2], i, cex=0.8, pos=4, offset=0.4)
  }
  points(Q$traj[L+1,1], Q$traj[L+1,2], pch=ifelse(Q$accept==1, 16, 1),
          col=ifelse(abs(Q$dH)>0.1, "red", "black"))
}

#9.8
HMC2 <- function (U, grad_U, epsilon, L, current_q) {
  q = current_q
  p = rnorm(length(q), 0, 1) # random flick - p is momentum.
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # initialize bookkeeping - saves trajectory
  qtraj <- matrix(NA,nrow=L+1,ncol=length(q))
  ptraj <- qtraj
  qtraj[1,] <- current_q
  ptraj[1,] <- p
  
  #9.9
  # Alternate full steps for position and momentum
  for (i in 1:L) {
    q = q + epsilon * p # Full step for the position
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) {
      p = p - epsilon * grad_U(q) # decrease momentum
      ptraj[i+1, ] <- p # update momentum
    }
    qtraj[i+1, ] <- q # update position
  }
  
  #9.10
  # Make a half step for momentum at the end
  p = p - epsilon * grad_U(q) / 2
  ptraj[L+1, ] <- p
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept <- 0
  if (runif(1) < exp(current_U - proposed_U + current_K-proposed_K)) { # exp(0) = 1 so if sum is zero accept, otherwise stochastic acceptance
    new_q <- q  # accept
    accept <- 1
  } else new_q <- current_q  # reject
  return(list(q=new_q, traj=qtraj, ptraj=ptraj, accept=accept))
}

#9.11
library(rethinking)
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000), ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa==1, 1, 2)

#9.12
m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std - 0.215),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
 ), data=dd)
precis(m8.3, depth=2)

#9.13
dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer(dd$cid)
)
str(dat_slim)

#9.14
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std - 0.215),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
 ), data=dat_slim, chains=1)

#9.15
precis(m9.1, depth=2)

#9.16
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std - 0.215),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
 ), data=dat_slim, chains=4, cores=4)

#9.17
show(m9.1)

#9.18
precis(m9.1, 2)

#9.19
pairs(m9.1)

#9.20
traceplot(m9.1, chains=1)

#9.21
trankplot(m9.1)

#9.22
y <- c(-1,1)
set.seed(11)
m9.2 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha,
    alpha ~ dnorm(0, 1000),
    sigma ~ dexp(0.0001)
 ), data=list(y=y), chains=3)

#9.23
precis(m9.2)

pairs(m9.2@stanfit)
traceplot(m9.2)
trankplot(m9.2)

#9.24
set.seed(11)
m9.3 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha,
    alpha ~ dnorm(1, 10),
    sigma ~ dexp(1)
 ), data=list(y=y), chains=3)
precis(m9.3)
traceplot(m9.3)
trankplot(m9.3)

#9.25
set.seed(41)
y <- rnorm(100, mean=0, sd=1)

#9.26
set.seed(384)
m9.4 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a1 + a2,
    a1 ~ dnorm(0, 1000),
    a2 ~ dnorm(0, 1000),
    sigma ~ dexp(1)
 ), data=list(y=y), chains=3)
precis(m9.4)
traceplot(m9.4)
trankplot(m9.4)

#9.27
m9.5 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a1 + a2,
    a1 ~ dnorm(0, 10),
    a2 ~ dnorm(0, 10),
    sigma ~ dexp(1)
 ), data=list(y=y), chains=3)
precis(m9.5)
traceplot(m9.5)
trankplot(m9.5)
