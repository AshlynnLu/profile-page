## Load function
source("src/helper-functions/functions.R")

## Load data
load(file = "data/raw/CR.sex.chicks09.RData")

CH <- as.matrix(CR.sex.chicks09[, -14])
n.ind <- dim(CR.sex.chicks09)[1]
n.occasion <- dim(CR.sex.chicks09[1:13])[2]
marr.full <- marray(CH)


## Constant model
caplik.constant <- function(theta, data){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data)[1]
  nj = dim(data)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi <- array(0,(nj-1))
  p <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  q <- array(0,dim=c(ni,nj))
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p[t] <- 1/(1+exp(-theta[1]))
    phi[t] <- 1/(1+exp(-theta[2]))
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q[t,t] <- phi[t]*p[t]
  }
  
  # Off diagonal elements
  
  for (t in 1:(ni-1)){
    for (j in (t+1):(nj-1)) {
      q[t,j] <- prod(phi[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q[t,nj] <- 1 - sum(q[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  
  likhood <- 0
  
  for (t in 1:ni){
    for (j in t:nj) {
      likhood <- likhood + data[t,j]*log(q[t,j])
    }
  }
  
  # Output the negative loglikelihood value:
  likhood <- -likhood
}

# MLE
opt.constant <- optim(par = c(0, 0),
                      fn = caplik.constant,
                      data = marr.full)

opt.constant.hess <- optimHess(par = c(0, 0),
                               fn = caplik.constant,
                               data = marr.full)

MLE.constant <- array(0, c(2, 3))
MLE.constant[1:2, 1] <- plogis(opt.constant$par)[1:2]

# SD
MLE.constant[1:2, 2] <- sqrt(diag(solve(opt.constant.hess)))[1:2]

MLE.constant <- round(MLE.constant, 3)

# resample
n.resample <- 50
resample.theta.constant <- matrix(data = 0, ncol = 2, nrow = n.resample)
for (i in 1:n.resample){
  resample <- CR.sex.chicks09[sample(dim(CR.sex.chicks09)[1], size = dim(CR.sex.chicks09)[1], replace = TRUE), -14]
  marr.resample <- marray(resample)
  
  opt.resample <- optim(par = c(0, 0),
                        fn = caplik.constant,
                        data = marr.resample)
  
  resample.theta.constant[i, ] = plogis(opt.resample$par)
}

# CI
MLE.constant[1, 3] <- paste0("[", round(quantile(resample.theta.constant[, 1], 0.025), 3), ",", round(quantile(resample.theta.constant[, 1], 0.975), 3), "]")
MLE.constant[2, 3] <- paste0("[", round(quantile(resample.theta.constant[, 2], 0.025), 3), ",", round(quantile(resample.theta.constant[, 2], 0.975), 3), "]")

# table
MLE.constant <- data.frame(MLE.constant)
rownames(MLE.constant) <- c("$p$", "$\\phi$")
colnames(MLE.constant) <- c("MLE", "SD", "CI")

# Save data
saveRDS(MLE.constant, file = "data/derived/MLE.constant.rds")


## Sex-dependent model
CH.m <- as.matrix(CR.sex.chicks09[which(CR.sex.chicks09$Sex == "M"), 1:13])
CH.f <- as.matrix(CR.sex.chicks09[which(CR.sex.chicks09$Sex == "F"), 1:13])
marr.m <- marray(CH.m)
marr.f <- marray(CH.f)

# Likelihood function
caplik.sex <- function(theta, data.m, data.f){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data.m)[1]
  nj = dim(data.m)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi.m <- array(0,nj-1)
  phi.f <- array(0,nj-1)
  
  p <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  
  q.m <- array(0,dim=c(ni,nj))
  q.f <- array(0,dim=c(ni,nj))
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p[t] <- 1/(1+exp(-theta[1]))
    phi.m[t] <- 1/(1+exp(-theta[2]))
    phi.f[t] <- 1/(1+exp(-theta[3]))
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q.m[t,t] <- phi.m[t]*p[t]
    q.f[t,t] <- phi.f[t]*p[t]
  }
  
  # Off diagonal elements
  
  for (t in 1:(ni-1)){
    for (j in (t+1):(nj-1)) {
      q.m[t,j] <- prod(phi.m[t:j])*prod(pbar[t:(j-1)])*p[j]
      q.f[t,j] <- prod(phi.f[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q.m[t,nj] <- 1 - sum(q.m[t,t:(nj-1)])
    q.f[t,nj] <- 1 - sum(q.f[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  likhood <- 0
  likhood.m <- 0
  likhood.f <- 0
  
  for (t in 1:ni){
    for (j in t:nj) {      
      likhood.m <- likhood.m + data.m[t,j]*log(q.m[t,j])
      likhood.f <- likhood.f + data.f[t,j]*log(q.f[t,j])
    }
  }
  
  likhood <- likhood.m + likhood.f
  
  # Output the negative loglikelihood value:
  likhood <- -likhood
}

# Find the MLE
opt.sex <- optim(par = c(0, 0, 0),
                 fn = caplik.sex,
                 data.m = marr.m, 
                 data.f = marr.f)

opt.sex.hess <- optimHess(par = c(0, 0, 0),
                          fn = caplik.sex,
                          data.m = marr.m, 
                          data.f = marr.f)
# MLE
MLE.sex <- array(0, c(3, 3))
MLE.sex[1:3, 1] <- plogis(opt.sex$par)[1:3]
# SD
MLE.sex[1:3, 2] <- sqrt(diag(solve(opt.sex.hess)))[1:3]

MLE.sex <- round(MLE.sex, 3)

# resample
resample.theta.sex <- matrix(data = 0, ncol = 3, nrow = n.resample)
for (i in 1:n.resample){
  resample <- CR.sex.chicks09[sample(dim(CR.sex.chicks09)[1], size = dim(CR.sex.chicks09)[1], replace = TRUE), ]
  resample.m <- as.matrix(resample[which(resample$Sex == "M"), -14])
  resample.f <- as.matrix(resample[which(resample$Sex == "F"), -14])
  
  marr.resample.m <- marray(resample.m)
  marr.resample.f <- marray(resample.f)
  
  opt.resample <- optim(par = c(0, 0, 0),
                        fn = caplik.sex,
                        data.m = marr.resample.m, 
                        data.f = marr.resample.f)
  
  resample.theta.sex[i, ] = plogis(opt.resample$par)
}

# CI
for (i in 1:3){
  MLE.sex[i, 3] <- paste0("[", round(quantile(resample.theta.sex[, i], 0.025), 3), ",", round(quantile(resample.theta.sex[, i], 0.975), 3), "]")
}


MLE.sex <- data.frame(MLE.sex)
rownames(MLE.sex) <- c("$p$", "$\\phi_m$", "$\\phi_f$")
colnames(MLE.sex) <- c("MLE", "SD", "CI")

# Save data
saveRDS(MLE.sex, file = "data/derived/MLE.sex.rds")

