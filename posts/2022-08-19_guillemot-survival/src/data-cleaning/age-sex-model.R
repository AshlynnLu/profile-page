## Load function
source("src/helper-functions/functions.R")

## Load data
load(file = "data/raw/CR.sex.chicks09.RData")

CH.m <- as.matrix(CR.sex.chicks09[which(CR.sex.chicks09$Sex == "M"), 1:13])
CH.f <- as.matrix(CR.sex.chicks09[which(CR.sex.chicks09$Sex == "F"), 1:13])
marr.m <- marray(CH.m)
marr.f <- marray(CH.f)

# Initialise an array to hold AIC
AIC.age.sex <- array(0, 3)

## 3 age class
marr1.age3.m <- marray.age(CH.m)[,,1]
marr2.age3.m <- marray.age(CH.m)[,,2]
marr3.age3.m <- marray.age(CH.m)[,,3]
marr1.age3.f <- marray.age(CH.f)[,,1]
marr2.age3.f <- marray.age(CH.f)[,,2]
marr3.age3.f <- marray.age(CH.f)[,,3]

caplik.age3.sex <- function(theta, data1.m, data2.m, data3.m, data1.f, data2.f, data3.f){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data1.m)[1]
  nj = dim(data1.m)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi1.m <- array(0,nj-1)
  phi2.m <- array(0,nj-1)
  phi3.m <- array(0,nj-1)
  phi1.f <- array(0,nj-1)
  phi2.f <- array(0,nj-1)
  phi3.f <- array(0,nj-1)
  
  p1 <- array(0,(nj-1))
  p <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  pbar1 <- array(0,(nj-1))
  
  q1.m <- array(0,dim=c(ni,nj))
  q2.m <- array(0,dim=c(ni,nj))
  q3.m <- array(0,dim=c(ni,nj))
  q1.f <- array(0,dim=c(ni,nj))
  q2.f <- array(0,dim=c(ni,nj))
  q3.f <- array(0,dim=c(ni,nj))
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p1[t] <- 0
    p[t] <- 1/(1+exp(-theta[1]))
    phi1.m[t] <- 1
    phi2.m[t] <- 1/(1+exp(-theta[2]))
    phi3.m[t] <- 1/(1+exp(-theta[3]))
    phi1.f[t] <- 1
    phi2.f[t] <- 1/(1+exp(-theta[4]))
    phi3.f[t] <- 1/(1+exp(-theta[5]))
    pbar1[t] <- 1-p1[t]
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q1.m[t,t] <- phi1.m[t]*p1[t]
    q2.m[t,t] <- phi2.m[t]*p[t]
    q3.m[t,t] <- phi3.m[t]*p[t]
    q1.f[t,t] <- phi1.f[t]*p1[t]
    q2.f[t,t] <- phi2.f[t]*p[t]
    q3.f[t,t] <- phi3.f[t]*p[t]
  }
  
  ## t+1
  ## Ages 1 and 2
  for (t in 1:(ni-1)){
    q1.m[t,t+1] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*p[t+1]
    q2.m[t,t+1] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*p[t+1]
    q1.f[t,t+1] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*p[t+1]
    q2.f[t,t+1] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*p[t+1]
    
    # Off diagonal elements (only age 3)
    for (j in (t+1):(nj-1)) {
      q3.m[t,j] <- prod(phi3.m[t:j])*prod(pbar[t:(j-1)])*p[j]
      q3.f[t,j] <- prod(phi3.f[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  ## t+2
  ## Age 2
  for (t in 1:(ni-2)){
    ## Age 2
    for (j in (t+2):(nj-1)) {
      q2.m[t,j] <- phi2.m[t]*pbar[t]*prod(phi3.m[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
      q2.f[t,j] <- phi2.f[t]*pbar[t]*prod(phi3.f[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
      
    }
    ## Age 1
    q1.m[t,t+2] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*p[t+2]
    q1.f[t,t+2] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*p[t+2]
  }
  
  ## t+3 for Age 1
  for (t in 1:(ni-3)){
    for (j in (t+3):(nj-1)) {
      q1.m[t,j] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*prod(phi3.m[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
      q1.f[t,j] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*prod(phi3.f[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
    }
  }
  
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q1.m[t,nj] <- 1 - sum(q1.m[t,t:(nj-1)])
    q2.m[t,nj] <- 1 - sum(q2.m[t,t:(nj-1)])
    q3.m[t,nj] <- 1 - sum(q3.m[t,t:(nj-1)])
    q1.f[t,nj] <- 1 - sum(q1.f[t,t:(nj-1)])
    q2.f[t,nj] <- 1 - sum(q2.f[t,t:(nj-1)])
    q3.f[t,nj] <- 1 - sum(q3.f[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  
  likhood <- 0
  likhood1.m <- 0
  likhood2.m <- 0
  likhood3.m <- 0
  likhood1.f <- 0
  likhood2.f <- 0
  likhood3.f <- 0
  
  for (i in 1:ni){
    for (j in i:nj){
      #likhood2.m <- likhood2.m + data2.m[i,j]*log(q2.m[i,j])
      likhood3.m <- likhood3.m + data3.m[i,j]*log(q3.m[i,j])
      #likhood2.f <- likhood2.f + data2.f[i,j]*log(q2.f[i,j])
      likhood3.f <- likhood3.f + data3.f[i,j]*log(q3.f[i,j])
    }
    for (j in (i+1):nj){
      likhood1.m <- likhood1.m + data1.m[i,j]*log(q1.m[i,j])
      likhood1.f <- likhood1.f + data1.f[i,j]*log(q1.f[i,j])
    }
  }
  
  likhood <- likhood1.m + likhood2.m + likhood3.m + likhood1.f + likhood2.f + likhood3.f
  
  # Output the negative loglikelihood value:
  likhood <- -likhood
}

opt.age3.sex <- optim(par = rep(0, 5),
                      fn = caplik.age3.sex,
                      data1.m = marr1.age3.m, 
                      data2.m = marr2.age3.m,
                      data3.m = marr3.age3.m,
                      data1.f = marr1.age3.f, 
                      data2.f = marr2.age3.f,
                      data3.f = marr3.age3.f,
                      hessian = TRUE)

opt.age3.sex.hess <- optimHess(par = rep(0, 5),
                               fn = caplik.age3.sex,
                               data1.m = marr1.age3.m, 
                               data2.m = marr2.age3.m,
                               data3.m = marr3.age3.m,
                               data1.f = marr1.age3.f, 
                               data2.f = marr2.age3.f,
                               data3.f = marr3.age3.f)

# MLE
estimate3.sex <- round(plogis(opt.age3.sex$par), 3)
MLE.age3.sex <- array("-", c(6, 4))
MLE.age3.sex[2, c(1, 3)] <- estimate3.sex[1]
MLE.age3.sex[3:4, 1] <- estimate3.sex[2:3]
MLE.age3.sex[3:4, 3] <- estimate3.sex[4:5]

# SD
sd3.sex <- round(sqrt(diag(solve(opt.age3.sex.hess))), 3)
MLE.age3.sex[2, c(2, 4)] <- sd3.sex[1]
MLE.age3.sex[3:4, 2] <- sd3.sex[2:3]
MLE.age3.sex[3:4, 4] <- sd3.sex[4:5]



# AIC = -2lik + 2*no.of param
AIC.age.sex[1] <- 2*opt.age3.sex$value + 2*5


## 4 age class
marr1.age4.m <- marray.age(CH.m, mAge = 4)[,,1]
marr2.age4.m <- marray.age(CH.m, mAge = 4)[,,2]
marr3.age4.m <- marray.age(CH.m, mAge = 4)[,,3]
marr4.age4.m <- marray.age(CH.m, mAge = 4)[,,4]
marr1.age4.f <- marray.age(CH.f, mAge = 4)[,,1]
marr2.age4.f <- marray.age(CH.f, mAge = 4)[,,2]
marr3.age4.f <- marray.age(CH.f, mAge = 4)[,,3]
marr4.age4.f <- marray.age(CH.f, mAge = 4)[,,4]

## Likelihood
caplik.age4.sex <- function(theta, data1.m, data2.m, data3.m, data4.m, data1.f, data2.f, data3.f, data4.f){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data1.m)[1]
  nj = dim(data1.m)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi1.m <- array(0,nj-1)
  phi2.m <- array(0,nj-1)
  phi3.m <- array(0,nj-1)
  phi4.m <- array(0,nj-1)
  phi1.f <- array(0,nj-1)
  phi2.f <- array(0,nj-1)
  phi3.f <- array(0,nj-1)
  phi4.f <- array(0,nj-1)
  
  p1 <- array(0,(nj-1))
  p <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  pbar1 <- array(0,(nj-1))
  
  q1.m <- array(0,dim=c(ni,nj))
  q2.m <- array(0,dim=c(ni,nj))
  q3.m <- array(0,dim=c(ni,nj))
  q4.m <- array(0,dim=c(ni,nj))
  q1.f <- array(0,dim=c(ni,nj))
  q2.f <- array(0,dim=c(ni,nj))
  q3.f <- array(0,dim=c(ni,nj))
  q4.f <- array(0,dim=c(ni,nj))
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p1[t] <- 0
    p[t] <- 1/(1+exp(-theta[1]))
    phi1.m[t] <- 1
    phi2.m[t] <- 1/(1+exp(-theta[2]))
    phi3.m[t] <- 1/(1+exp(-theta[3]))
    phi4.m[t] <- 1/(1+exp(-theta[4]))
    phi1.f[t] <- 1
    phi2.f[t] <- 1/(1+exp(-theta[5]))
    phi3.f[t] <- 1/(1+exp(-theta[6]))
    phi4.f[t] <- 1/(1+exp(-theta[7]))
    pbar1[t] <- 1-p1[t]
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q1.m[t,t] <- phi1.m[t]*p1[t]
    q2.m[t,t] <- phi2.m[t]*p[t]
    q3.m[t,t] <- phi3.m[t]*p[t]
    q4.m[t,t] <- phi4.m[t]*p[t]
    q1.f[t,t] <- phi1.f[t]*p1[t]
    q2.f[t,t] <- phi2.f[t]*p[t]
    q3.f[t,t] <- phi3.f[t]*p[t]
    q4.f[t,t] <- phi4.f[t]*p[t]
  }
  
  ## t+1
  ## Ages 1, 2 and 3
  for (t in 1:(ni-1)){
    q1.m[t,t+1] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*p[t+1]
    q2.m[t,t+1] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*p[t+1]
    q3.m[t,t+1] <- phi3.m[t]*pbar[t]*phi4.m[t+1]*p[t+1]
    q1.f[t,t+1] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*p[t+1]
    q2.f[t,t+1] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*p[t+1]
    q3.f[t,t+1] <- phi3.f[t]*pbar[t]*phi4.f[t+1]*p[t+1]
    
    # Off diagonal elements (only age 4)
    for (j in (t+1):(nj-1)) {
      q4.m[t,j] <- prod(phi4.m[t:j])*prod(pbar[t:(j-1)])*p[j]
      q4.f[t,j] <- prod(phi4.f[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  ## t+2
  for (t in 1:(ni-2)){
    ## Age 3
    for (j in (t+2):(nj-1)) {
      q3.m[t,j] <- phi3.m[t]*pbar[t]*prod(phi4.m[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
      q3.f[t,j] <- phi3.f[t]*pbar[t]*prod(phi4.f[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
    }
    ## Age 1 and 2
    q1.m[t,t+2] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*p[t+2]
    q2.m[t,t+2] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*phi4.m[t+2]*p[t+2]
    q1.f[t,t+2] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*p[t+2]
    q2.f[t,t+2] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*phi4.f[t+2]*p[t+2]
  }
  
  ## t+3 for Age 1 and 2
  for (t in 1:(ni-3)){
    for (j in (t+3):(nj-1)) {
      q2.m[t,j] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*prod(phi4.m[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
      q2.f[t,j] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*prod(phi4.f[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
    }
    q1.m[t,t+3] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*pbar[t+2]*phi4.m[t+3]*p[t+3]
    q1.f[t,t+3] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*pbar[t+2]*phi4.f[t+3]*p[t+3]
  }
  
  ## t+4 for Age 1
  for (t in 1:(ni-4)){
    for (j in (t+4):(nj-1)) {
      q1.m[t,j] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*pbar[t+2]*prod(phi4.m[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
      q1.f[t,j] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*pbar[t+2]*prod(phi4.f[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
    }
  }
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q1.m[t,nj] <- 1 - sum(q1.m[t,t:(nj-1)])
    q2.m[t,nj] <- 1 - sum(q2.m[t,t:(nj-1)])
    q3.m[t,nj] <- 1 - sum(q3.m[t,t:(nj-1)])
    q4.m[t,nj] <- 1 - sum(q4.m[t,t:(nj-1)])
    q1.f[t,nj] <- 1 - sum(q1.f[t,t:(nj-1)])
    q2.f[t,nj] <- 1 - sum(q2.f[t,t:(nj-1)])
    q3.f[t,nj] <- 1 - sum(q3.f[t,t:(nj-1)])
    q4.f[t,nj] <- 1 - sum(q4.f[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  
  likhood <- 0
  likhood1.m <- 0
  likhood2.m <- 0
  likhood3.m <- 0
  likhood4.m <- 0
  likhood1.f <- 0
  likhood2.f <- 0
  likhood3.f <- 0
  likhood4.f <- 0
  
  for (i in 1:ni){
    for (j in i:nj) {
      likhood2.m <- likhood2.m + data2.m[i,j]*log(q2.m[i,j])
      likhood3.m <- likhood3.m + data3.m[i,j]*log(q3.m[i,j])
      likhood4.m <- likhood4.m + data4.m[i,j]*log(q4.m[i,j])
      likhood2.f <- likhood2.f + data2.f[i,j]*log(q2.f[i,j])
      likhood3.f <- likhood3.f + data3.f[i,j]*log(q3.f[i,j])
      likhood4.f <- likhood4.f + data4.f[i,j]*log(q4.f[i,j])
    }
    for (j in (i+1):nj){
      likhood1.m <- likhood1.m + data1.m[i,j]*log(q1.m[i,j])
      likhood1.f <- likhood1.f + data1.f[i,j]*log(q1.f[i,j])
    }
  }
  
  likhood <- likhood1.m + likhood2.m + likhood3.m +likhood4.m + likhood1.f + likhood2.f + likhood3.f + likhood4.f
  
  # Output the negative loglikelihood value:
  likhood<- -likhood
}

opt.age4.sex <- optim(par = rep(0, 7),
                      fn = caplik.age4.sex,
                      data1.m = marr1.age4.m, 
                      data2.m = marr2.age4.m,
                      data3.m = marr3.age4.m,
                      data4.m = marr4.age4.m,
                      data1.f = marr1.age4.f, 
                      data2.f = marr2.age4.f,
                      data3.f = marr3.age4.f,
                      data4.f = marr4.age4.f, 
                      hessian = TRUE)

opt.age4.sex.hess <- optimHess(par = rep(0, 7),
                               fn = caplik.age4.sex,
                               data1.m = marr1.age4.m, 
                               data2.m = marr2.age4.m,
                               data3.m = marr3.age4.m,
                               data4.m = marr4.age4.m,
                               data1.f = marr1.age4.f, 
                               data2.f = marr2.age4.f,
                               data3.f = marr3.age4.f,
                               data4.f = marr4.age4.f)

estimate4.sex <- round(plogis(opt.age4.sex$par), 3)
MLE.age4.sex <- array("-", c(6, 4))
MLE.age4.sex[2, c(1, 3)] <- estimate4.sex[1]
MLE.age4.sex[3:5, 1] <- estimate4.sex[2:4]
MLE.age4.sex[3:5, 3] <- estimate4.sex[5:7]

sd4.sex <- round(sqrt(diag(solve(opt.age4.sex.hess))), 3)
MLE.age4.sex[2, c(2, 4)] <- sd4.sex[1]
MLE.age4.sex[3:5, 2] <- sd4.sex[2:4]
MLE.age4.sex[3:5, 4] <- sd4.sex[5:7]


# AIC = -2lik + 2*no.of param
AIC.age.sex[2] <- 2*opt.age4.sex$value + 2*7



## 5 age class
marr1.age5.m <- marray.age(CH.m, mAge = 5)[,,1]
marr2.age5.m <- marray.age(CH.m, mAge = 5)[,,2]
marr3.age5.m <- marray.age(CH.m, mAge = 5)[,,3]
marr4.age5.m <- marray.age(CH.m, mAge = 5)[,,4]
marr5.age5.m <- marray.age(CH.m, mAge = 5)[,,5]
marr1.age5.f <- marray.age(CH.f, mAge = 5)[,,1]
marr2.age5.f <- marray.age(CH.f, mAge = 5)[,,2]
marr3.age5.f <- marray.age(CH.f, mAge = 5)[,,3]
marr4.age5.f <- marray.age(CH.f, mAge = 5)[,,4]
marr5.age5.f <- marray.age(CH.f, mAge = 5)[,,5]

## Likelihood
caplik.age5.sex <- function(theta, data1.m, data2.m, data3.m, data4.m, data5.m, data1.f, data2.f, data3.f, data4.f, data5.f){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data1.m)[1]
  nj = dim(data1.m)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi1.m <- array(0,nj-1)
  phi2.m <- array(0,nj-1)
  phi3.m <- array(0,nj-1)
  phi4.m <- array(0,nj-1)
  phi5.m <- array(0,nj-1)
  
  phi1.f <- array(0,nj-1)
  phi2.f <- array(0,nj-1)
  phi3.f <- array(0,nj-1)
  phi4.f <- array(0,nj-1)
  phi5.f <- array(0,nj-1)
  
  p1 <- array(0,(nj-1))
  p <- array(0,(nj-1))
  pbar1 <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  
  q1.m <- array(0,dim=c(ni,nj))
  q2.m <- array(0,dim=c(ni,nj))
  q3.m <- array(0,dim=c(ni,nj))
  q4.m <- array(0,dim=c(ni,nj))
  q5.m <- array(0,dim=c(ni,nj))
  
  q1.f <- array(0,dim=c(ni,nj))
  q2.f <- array(0,dim=c(ni,nj))
  q3.f <- array(0,dim=c(ni,nj))
  q4.f <- array(0,dim=c(ni,nj))
  q5.f <- array(0,dim=c(ni,nj))
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p1[t] <- 0
    p[t] <- 1/(1+exp(-theta[1]))
    
    phi1.m[t] <- 1
    phi2.m[t] <- 1/(1+exp(-theta[2]))
    phi3.m[t] <- 1/(1+exp(-theta[3]))
    phi4.m[t] <- 1/(1+exp(-theta[4]))
    phi5.m[t] <- 1/(1+exp(-theta[5]))
    
    phi1.f[t] <- 1
    phi2.f[t] <- 1/(1+exp(-theta[6]))
    phi3.f[t] <- 1/(1+exp(-theta[7]))
    phi4.f[t] <- 1/(1+exp(-theta[8]))
    phi5.f[t] <- 1/(1+exp(-theta[9]))
    
    pbar1[t] <- 1-p1[t]
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q1.m[t,t] <- phi1.m[t]*p1[t]
    q2.m[t,t] <- phi2.m[t]*p[t]
    q3.m[t,t] <- phi3.m[t]*p[t]
    q4.m[t,t] <- phi4.m[t]*p[t]
    q5.m[t,t] <- phi5.m[t]*p[t]
    
    q1.f[t,t] <- phi1.f[t]*p1[t]
    q2.f[t,t] <- phi2.f[t]*p[t]
    q3.f[t,t] <- phi3.f[t]*p[t]
    q4.f[t,t] <- phi4.f[t]*p[t]
    q5.f[t,t] <- phi5.f[t]*p[t]
  }
  
  ## t+1
  ## Ages 1, 2, 3 and 4
  for (t in 1:(ni-1)){
    q1.m[t,t+1] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*p[t+1]
    q2.m[t,t+1] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*p[t+1]
    q3.m[t,t+1] <- phi3.m[t]*pbar[t]*phi4.m[t+1]*p[t+1]
    q4.m[t,t+1] <- phi4.m[t]*pbar[t]*phi5.m[t+1]*p[t+1]
    
    q1.f[t,t+1] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*p[t+1]
    q2.f[t,t+1] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*p[t+1]
    q3.f[t,t+1] <- phi3.f[t]*pbar[t]*phi4.f[t+1]*p[t+1]
    q4.f[t,t+1] <- phi4.f[t]*pbar[t]*phi5.f[t+1]*p[t+1]
    
    # Off diagonal elements (only age 5)
    for (j in (t+1):(nj-1)) {
      q5.m[t,j] <- prod(phi5.m[t:j])*prod(pbar[t:(j-1)])*p[j]
      q5.f[t,j] <- prod(phi5.f[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  ## t+2
  for (t in 1:(ni-2)){
    ## Age 4
    for (j in (t+2):(nj-1)) {
      q4.m[t,j] <- phi4.m[t]*pbar[t]*prod(phi5.m[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
      q4.f[t,j] <- phi4.f[t]*pbar[t]*prod(phi5.f[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
    }
    ## Age 1, 2 and 3
    q1.m[t,t+2] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*p[t+2]
    q2.m[t,t+2] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*phi4.m[t+2]*p[t+2]
    q3.m[t,t+2] <- phi3.m[t]*pbar[t]*phi4.m[t+1]*pbar[t+1]*phi5.m[t+2]*p[t+2]
    q1.f[t,t+2] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*p[t+2]
    q2.f[t,t+2] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*phi4.f[t+2]*p[t+2]
    q3.f[t,t+2] <- phi3.f[t]*pbar[t]*phi4.f[t+1]*pbar[t+1]*phi5.f[t+2]*p[t+2]
  }
  
  ## t+3 for Age 1, 2 and 3
  for (t in 1:(ni-3)){
    for (j in (t+3):(nj-1)) {
      q3.m[t,j] <- phi3.m[t]*pbar[t]*phi4.m[t+1]*pbar[t+1]*prod(phi5.m[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
      q3.f[t,j] <- phi3.f[t]*pbar[t]*phi4.f[t+1]*pbar[t+1]*prod(phi5.f[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
    }
    q1.m[t,t+3] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*p[t+2]*phi4.m[t+3]*p[t+3]
    q1.f[t,t+3] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*p[t+2]*phi4.f[t+3]*p[t+3]
    q2.m[t,t+3] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*phi4.m[t+2]*p[t+2]*phi5.m[t+3]*p[t+3]
    q2.f[t,t+3] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*phi4.f[t+2]*p[t+2]*phi5.f[t+3]*p[t+3]
  }
  
  ## t+4 for Age 1 and 2
  for (t in 1:(ni-4)){
    for (j in (t+4):(nj-1)) {
      q2.m[t,j] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*phi4.m[t+2]*pbar[t+2]*prod(phi5.m[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
      q2.f[t,j] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*phi4.f[t+2]*pbar[t+2]*prod(phi5.f[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
    }
    q1.m[t,t+4] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*pbar[t+2]*phi4.m[t+3]*pbar[t+3]*phi5.m[t+4]*p[t+4]
    q1.f[t,t+4] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*pbar[t+2]*phi4.f[t+3]*pbar[t+3]*phi5.f[t+4]*p[t+4]
  }
  
  ## t+5 for Age 1
  for (t in 1:(ni-5)){
    for (j in (t+5):(nj-1)) {
      q1.m[t,j] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*pbar[t+2]*phi4.m[t+3]*pbar[t+3]*prod(phi5.m[(t+4):j])*prod(pbar[(t+4):(j-1)])*p[j]
      q1.f[t,j] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*pbar[t+2]*phi4.f[t+3]*pbar[t+3]*prod(phi5.f[(t+4):j])*prod(pbar[(t+4):(j-1)])*p[j]
    }
  }
  
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q1.m[t,nj] <- 1 - sum(q1.m[t,t:(nj-1)])
    q2.m[t,nj] <- 1 - sum(q2.m[t,t:(nj-1)])
    q3.m[t,nj] <- 1 - sum(q3.m[t,t:(nj-1)])
    q4.m[t,nj] <- 1 - sum(q4.m[t,t:(nj-1)])
    q5.m[t,nj] <- 1 - sum(q5.m[t,t:(nj-1)])
    
    q1.f[t,nj] <- 1 - sum(q1.f[t,t:(nj-1)])
    q2.f[t,nj] <- 1 - sum(q2.f[t,t:(nj-1)])
    q3.f[t,nj] <- 1 - sum(q3.f[t,t:(nj-1)])
    q4.f[t,nj] <- 1 - sum(q4.f[t,t:(nj-1)])
    q5.f[t,nj] <- 1 - sum(q5.f[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  
  likhood <- 0
  likhood1.m <- 0
  likhood2.m <- 0
  likhood3.m <- 0
  likhood4.m <- 0
  likhood5.m <- 0
  likhood1.f <- 0
  likhood2.f <- 0
  likhood3.f <- 0
  likhood4.f <- 0
  likhood5.f <- 0
  
  for (i in 1:ni){
    for (j in i:nj) {
      likhood2.m <- likhood2.m + data2.m[i,j]*log(q2.m[i,j])
      likhood3.m <- likhood3.m + data3.m[i,j]*log(q3.m[i,j])
      likhood4.m <- likhood4.m + data4.m[i,j]*log(q4.m[i,j])
      likhood5.m <- likhood5.m + data5.m[i,j]*log(q5.m[i,j])
      likhood2.f <- likhood2.f + data2.f[i,j]*log(q2.f[i,j])
      likhood3.f <- likhood3.f + data3.f[i,j]*log(q3.f[i,j])
      likhood4.f <- likhood4.f + data4.f[i,j]*log(q4.f[i,j])
      likhood5.f <- likhood5.f + data5.f[i,j]*log(q5.f[i,j])
    }
    for (j in (i+1):nj){
      likhood1.m <- likhood1.m + data1.m[i,j]*log(q1.m[i,j])
      likhood1.f <- likhood1.f + data1.f[i,j]*log(q1.f[i,j])
    }
  }
  
  likhood <- likhood1.m + likhood2.m + likhood3.m +likhood4.m + likhood5.m + likhood1.f + likhood2.f + likhood3.f + likhood4.f +likhood5.f
  
  # Output the negative loglikelihood value:
  likhood<- -likhood
}

opt.age5.sex <- optim(par = rep(0, 9),
                      fn = caplik.age5.sex,
                      data1.m = marr1.age5.m, 
                      data2.m = marr2.age5.m,
                      data3.m = marr3.age5.m,
                      data4.m = marr4.age5.m,
                      data5.m = marr5.age5.m,
                      data1.f = marr1.age5.f, 
                      data2.f = marr2.age5.f,
                      data3.f = marr3.age5.f,
                      data4.f = marr4.age5.f,
                      data5.f = marr5.age5.f, 
                      hessian = TRUE)

opt.age5.sex.hess <- optimHess(par = rep(0, 9),
                               fn = caplik.age5.sex,
                               data1.m = marr1.age5.m, 
                               data2.m = marr2.age5.m,
                               data3.m = marr3.age5.m,
                               data4.m = marr4.age5.m,
                               data5.m = marr5.age5.m,
                               data1.f = marr1.age5.f, 
                               data2.f = marr2.age5.f,
                               data3.f = marr3.age5.f,
                               data4.f = marr4.age5.f,
                               data5.f = marr5.age5.f)

estimate5.sex <- round(plogis(opt.age5.sex$par), 3)
MLE.age5.sex <- array("-", c(6, 4))
MLE.age5.sex[2, c(1, 3)] <- estimate5.sex[1]
MLE.age5.sex[3:6, 1] <- estimate5.sex[2:5]
MLE.age5.sex[3:6, 3] <- estimate5.sex[6:9]

sd5.sex <- round(sqrt(diag(solve(opt.age5.sex.hess))), 3)
MLE.age5.sex[2, c(2, 4)] <- sd4.sex[1]
MLE.age5.sex[3:6, 2] <- sd5.sex[2:5]
MLE.age5.sex[3:6, 4] <- sd5.sex[6:9]

# AIC = -2lik + 2*no.of param
AIC.age.sex[3] <- 2*opt.age5.sex$value + 2*9


# Table
MLE.age.sex <- array(0, c(7, 12))
MLE.age.sex[-1, 1:4] <- MLE.age3.sex
MLE.age.sex[-1, 5:8] <- MLE.age4.sex
MLE.age.sex[-1, 9:12] <- MLE.age5.sex
MLE.age.sex[1, ] <- c("male", " ", "female", " ", "male", " ", "female", " ", "male", " ", "female", " ")
MLE.age.sex[2, ] <- c("MLE", "SD", "MLE", "SD","MLE", "SD","MLE", "SD","MLE", "SD","MLE", "SD")
MLE.age.sex <- data.frame(MLE.age.sex, row.names = c("sex", " ", "$p$", "$\\phi_2$", "$\\phi_3$", "$\\phi_4$", "$\\phi_5$"))
colnames(MLE.age.sex) <- c("3 age class", " ", " ", " ", "4 age class", " ", " ", " ", "5 age class", " ", " ", " ")

# Save data
saveRDS(MLE.age.sex, file = "data/derived/MLE.age.sex.rds")



## 6 age classes
marr1.age6.m <- marray.age(CH.m, mAge = 6)[,,1]
marr2.age6.m <- marray.age(CH.m, mAge = 6)[,,2]
marr3.age6.m <- marray.age(CH.m, mAge = 6)[,,3]
marr4.age6.m <- marray.age(CH.m, mAge = 6)[,,4]
marr5.age6.m <- marray.age(CH.m, mAge = 6)[,,5]
marr6.age6.m <- marray.age(CH.m, mAge = 6)[,,6]

marr1.age6.f <- marray.age(CH.f, mAge = 6)[,,1]
marr2.age6.f <- marray.age(CH.f, mAge = 6)[,,2]
marr3.age6.f <- marray.age(CH.f, mAge = 6)[,,3]
marr4.age6.f <- marray.age(CH.f, mAge = 6)[,,4]
marr5.age6.f <- marray.age(CH.f, mAge = 6)[,,5]
marr6.age6.f <- marray.age(CH.f, mAge = 6)[,,6]

## Likelihood
caplik.age6.sex <- function(theta, data1.m, data2.m, data3.m, data4.m, data5.m, data6.m, data1.f, data2.f, data3.f, data4.f, data5.f, data6.f){
  # Number of release occasions and recovery occasions
  
  ni = dim(data1.m)[1]
  nj = dim(data1.m)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi1.m <- array(0,nj-1)
  phi2.m <- array(0,nj-1)
  phi3.m <- array(0,nj-1)
  phi4.m <- array(0,nj-1)
  phi5.m <- array(0,nj-1)
  phi6.m <- array(0,nj-1)
  
  phi1.f <- array(0,nj-1)
  phi2.f <- array(0,nj-1)
  phi3.f <- array(0,nj-1)
  phi4.f <- array(0,nj-1)
  phi5.f <- array(0,nj-1)
  phi6.f <- array(0,nj-1)
  
  p1 <- array(0,(nj-1))
  p <- array(0,(nj-1))
  pbar1 <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  
  q1.m <- array(0,dim=c(ni,nj))
  q2.m <- array(0,dim=c(ni,nj))
  q3.m <- array(0,dim=c(ni,nj))
  q4.m <- array(0,dim=c(ni,nj))
  q5.m <- array(0,dim=c(ni,nj))
  q6.m <- array(0,dim=c(ni,nj))
  
  q1.f <- array(0,dim=c(ni,nj))
  q2.f <- array(0,dim=c(ni,nj))
  q3.f <- array(0,dim=c(ni,nj))
  q4.f <- array(0,dim=c(ni,nj))
  q5.f <- array(0,dim=c(ni,nj))
  q6.f <- array(0,dim=c(ni,nj))
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p1[t] <- 0
    p[t] <- 1/(1+exp(-theta[1]))
    
    phi1.m[t] <- 1
    phi2.m[t] <- 1/(1+exp(-theta[2]))
    phi3.m[t] <- 1/(1+exp(-theta[3]))
    phi4.m[t] <- 1/(1+exp(-theta[4]))
    phi5.m[t] <- 1/(1+exp(-theta[5]))
    phi6.m[t] <- 1/(1+exp(-theta[6]))
    
    phi1.f[t] <- 1
    phi2.f[t] <- 1/(1+exp(-theta[7]))
    phi3.f[t] <- 1/(1+exp(-theta[8]))
    phi4.f[t] <- 1/(1+exp(-theta[9]))
    phi5.f[t] <- 1/(1+exp(-theta[10]))
    phi6.f[t] <- 1/(1+exp(-theta[11]))
    
    pbar1[t] <- 1-p1[t]
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q1.m[t,t] <- phi1.m[t]*p1[t]
    q2.m[t,t] <- phi2.m[t]*p[t]
    q3.m[t,t] <- phi3.m[t]*p[t]
    q4.m[t,t] <- phi4.m[t]*p[t]
    q5.m[t,t] <- phi5.m[t]*p[t]
    q6.m[t,t] <- phi6.m[t]*p[t]
    
    q1.f[t,t] <- phi1.f[t]*p1[t]
    q2.f[t,t] <- phi2.f[t]*p[t]
    q3.f[t,t] <- phi3.f[t]*p[t]
    q4.f[t,t] <- phi4.f[t]*p[t]
    q5.f[t,t] <- phi5.f[t]*p[t]
    q6.f[t,t] <- phi6.f[t]*p[t]
  }
  
  ## t+1
  ## Ages 1, 2, 3, 4 and 6
  for (t in 1:(ni-1)){
    q1.m[t,t+1] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*p[t+1]
    q2.m[t,t+1] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*p[t+1]
    q3.m[t,t+1] <- phi3.m[t]*pbar[t]*phi4.m[t+1]*p[t+1]
    q4.m[t,t+1] <- phi4.m[t]*pbar[t]*phi5.m[t+1]*p[t+1]
    q5.m[t,t+1] <- phi5.m[t]*pbar[t]*phi6.m[t+1]*p[t+1]
    
    q1.f[t,t+1] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*p[t+1]
    q2.f[t,t+1] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*p[t+1]
    q3.f[t,t+1] <- phi3.f[t]*pbar[t]*phi4.f[t+1]*p[t+1]
    q4.f[t,t+1] <- phi4.f[t]*pbar[t]*phi5.f[t+1]*p[t+1]
    q5.f[t,t+1] <- phi5.f[t]*pbar[t]*phi6.f[t+1]*p[t+1]
    
    # Off diagonal elements (only age 6)
    for (j in (t+1):(nj-1)) {
      q6.m[t,j] <- prod(phi6.m[t:j])*prod(pbar[t:(j-1)])*p[j]
      q6.f[t,j] <- prod(phi6.f[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  ## t+2
  for (t in 1:(ni-2)){
    ## Age 5
    for (j in (t+2):(nj-1)) {
      q5.m[t,j] <- phi5.m[t]*pbar[t]*prod(phi6.m[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
      q5.f[t,j] <- phi5.f[t]*pbar[t]*prod(phi6.f[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
    }
    ## Age 1, 2, 3 and 4
    q1.m[t,t+2] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*p[t+2]
    q2.m[t,t+2] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*phi4.m[t+2]*p[t+2]
    q3.m[t,t+2] <- phi3.m[t]*pbar[t]*phi4.m[t+1]*pbar[t+1]*phi5.m[t+2]*p[t+2]
    q4.m[t,t+2] <- phi4.m[t]*pbar[t]*phi5.m[t+1]*pbar[t+1]*phi6.m[t+2]*p[t+2]
    
    q1.f[t,t+2] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*p[t+2]
    q2.f[t,t+2] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*phi4.f[t+2]*p[t+2]
    q3.f[t,t+2] <- phi3.f[t]*pbar[t]*phi4.f[t+1]*pbar[t+1]*phi5.f[t+2]*p[t+2]
    q4.f[t,t+2] <- phi4.f[t]*pbar[t]*phi5.f[t+1]*pbar[t+1]*phi6.f[t+2]*p[t+2]
  }
  
  ## t+3
  for (t in 1:(ni-3)){
    ## Age 4
    for (j in (t+3):(nj-1)) {
      q4.m[t,j] <- phi4.m[t]*pbar[t]*phi5.m[t+1]*pbar[t+1]*prod(phi6.m[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
      q4.f[t,j] <- phi4.f[t]*pbar[t]*phi5.f[t+1]*pbar[t+1]*prod(phi6.f[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
    }
    ## Age 1, 2 and 3 
    q1.m[t,t+3] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*pbar[t+2]*phi4.m[t+3]*p[t+3]
    q2.m[t,t+3] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*phi4.m[t+2]*pbar[t+2]*phi5.m[t+3]*p[t+3]
    q3.m[t,t+3] <- phi3.m[t]*pbar[t]*phi4.m[t+1]*pbar[t+1]*phi5.m[t+2]*pbar[t+2]*phi6.m[t+3]*p[t+3]
    
    q1.f[t,t+3] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*pbar[t+2]*phi4.f[t+3]*p[t+3]
    q2.f[t,t+3] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*phi4.f[t+2]*pbar[t+2]*phi5.f[t+3]*p[t+3]
    q3.f[t,t+3] <- phi3.f[t]*pbar[t]*phi4.f[t+1]*pbar[t+1]*phi5.f[t+2]*pbar[t+2]*phi6.f[t+3]*p[t+3]
  }
  
  ## t+4 
  for (t in 1:(ni-4)){
    ## Age 3
    for (j in (t+4):(nj-1)) {
      q3.m[t,j] <- phi3.m[t]*pbar[t]*phi4.m[t+1]*pbar[t+1]*phi5.m[t+2]*pbar[t+2]*prod(phi6.m[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
      q3.f[t,j] <- phi3.f[t]*pbar[t]*phi4.f[t+1]*pbar[t+1]*phi5.f[t+2]*pbar[t+2]*prod(phi6.f[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
    }
    ## Age 1 and 2
    q1.m[t,t+4] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*pbar[t+2]*phi4.m[t+3]*pbar[t+3]*phi5.m[t+4]*p[t+4]
    q2.m[t,t+4] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*phi4.m[t+2]*pbar[t+2]*phi5.m[t+3]*pbar[t+3]*phi6.m[t+4]*p[t+4]
    
    q1.f[t,t+4] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*pbar[t+2]*phi4.f[t+3]*pbar[t+3]*phi5.f[t+4]*p[t+4]
    q2.f[t,t+4] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*phi4.f[t+2]*pbar[t+2]*phi5.f[t+3]*pbar[t+3]*phi6.f[t+4]*p[t+4]
  }
  
  ## t+5 
  for (t in 1:(ni-5)){
    ## Age 2
    for (j in (t+5):(nj-1)) {
      q2.m[t,j] <- phi2.m[t]*pbar[t]*phi3.m[t+1]*pbar[t+1]*phi4.m[t+2]*pbar[t+2]*phi5.m[t+3]*pbar[t+3]*prod(phi6.m[(t+4):j])*prod(pbar[(t+4):(j-1)])*p[j]
      q2.f[t,j] <- phi2.f[t]*pbar[t]*phi3.f[t+1]*pbar[t+1]*phi4.f[t+2]*pbar[t+2]*phi5.f[t+3]*pbar[t+3]*prod(phi6.f[(t+4):j])*prod(pbar[(t+4):(j-1)])*p[j]
    }
    ## Age 1
    q1.m[t,t+5] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*pbar[t+2]*phi4.m[t+3]*pbar[t+3]*phi5.m[t+4]*pbar[t+4]*phi6.m[t+5]*p[t+5]
    q1.f[t,t+5] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*pbar[t+2]*phi4.f[t+3]*pbar[t+3]*phi5.f[t+4]*pbar[t+4]*phi6.f[t+5]*p[t+5]
  }
  
  ## t+6
  ## Age 1
  for (t in 1:(ni-6)){
    for (j in (t+6):(nj-1)) {
      q1.m[t,j] <- phi1.m[t]*pbar1[t]*phi2.m[t+1]*pbar[t+1]*phi3.m[t+2]*pbar[t+2]*phi4.m[t+3]*pbar[t+3]*phi5.m[t+4]*pbar[t+4]*prod(phi6.m[(t+5):j])*prod(pbar[(t+5):(j-1)])*p[j]
      q1.f[t,j] <- phi1.f[t]*pbar1[t]*phi2.f[t+1]*pbar[t+1]*phi3.f[t+2]*pbar[t+2]*phi4.f[t+3]*pbar[t+3]*phi5.f[t+4]*pbar[t+4]*prod(phi6.f[(t+5):j])*prod(pbar[(t+5):(j-1)])*p[j]
    }
  }
  
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q1.m[t,nj] <- 1 - sum(q1.m[t,t:(nj-1)])
    q2.m[t,nj] <- 1 - sum(q2.m[t,t:(nj-1)])
    q3.m[t,nj] <- 1 - sum(q3.m[t,t:(nj-1)])
    q4.m[t,nj] <- 1 - sum(q4.m[t,t:(nj-1)])
    q5.m[t,nj] <- 1 - sum(q5.m[t,t:(nj-1)])
    q6.m[t,nj] <- 1 - sum(q6.m[t,t:(nj-1)])
    
    q1.f[t,nj] <- 1 - sum(q1.f[t,t:(nj-1)])
    q2.f[t,nj] <- 1 - sum(q2.f[t,t:(nj-1)])
    q3.f[t,nj] <- 1 - sum(q3.f[t,t:(nj-1)])
    q4.f[t,nj] <- 1 - sum(q4.f[t,t:(nj-1)])
    q5.f[t,nj] <- 1 - sum(q5.f[t,t:(nj-1)])
    q6.f[t,nj] <- 1 - sum(q6.f[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  
  likhood <- 0
  likhood1.m <- 0
  likhood2.m <- 0
  likhood3.m <- 0
  likhood4.m <- 0
  likhood5.m <- 0
  likhood6.m <- 0
  
  likhood1.f <- 0
  likhood2.f <- 0
  likhood3.f <- 0
  likhood4.f <- 0
  likhood5.f <- 0
  likhood6.f <- 0
  
  for (i in 1:ni){
    for (j in i:nj) {
      likhood2.m <- likhood2.m + data2.m[i,j]*log(q2.m[i,j])
      likhood3.m <- likhood3.m + data3.m[i,j]*log(q3.m[i,j])
      likhood4.m <- likhood4.m + data4.m[i,j]*log(q4.m[i,j])
      likhood5.m <- likhood5.m + data5.m[i,j]*log(q5.m[i,j])
      likhood6.m <- likhood6.m + data6.m[i,j]*log(q6.m[i,j])
      
      likhood2.f <- likhood2.f + data2.f[i,j]*log(q2.f[i,j])
      likhood3.f <- likhood3.f + data3.f[i,j]*log(q3.f[i,j])
      likhood4.f <- likhood4.f + data4.f[i,j]*log(q4.f[i,j])
      likhood5.f <- likhood5.f + data5.f[i,j]*log(q5.f[i,j])
      likhood6.f <- likhood6.f + data6.f[i,j]*log(q6.f[i,j])
    }
    for (j in (i+1):nj){
      likhood1.m <- likhood1.m + data1.m[i,j]*log(q1.m[i,j])
      likhood1.f <- likhood1.f + data1.f[i,j]*log(q1.f[i,j])
    }
  }
  
  likhood <- likhood1.m + likhood2.m + likhood3.m +likhood4.m + likhood5.m + likhood6.m + likhood1.f + likhood2.f + likhood3.f + likhood4.f + likhood5.f + likhood6.f
  
  # Output the negative loglikelihood value:
  likhood <- -likhood
}

opt.age6.sex <- optim(par = rep(0, 11),
                      fn = caplik.age6.sex,
                      data1.m = marr1.age6.m, 
                      data2.m = marr2.age6.m,
                      data3.m = marr3.age6.m,
                      data4.m = marr4.age6.m,
                      data5.m = marr5.age6.m,
                      data6.m = marr6.age6.m,
                      data1.f = marr1.age6.f, 
                      data2.f = marr2.age6.f,
                      data3.f = marr3.age6.f,
                      data4.f = marr4.age6.f,
                      data5.f = marr5.age6.f, 
                      data6.f = marr6.age6.f)

# p, male(phi 2, 3, 4, 5, 6), female(phi 2, 3, 4, 5, 6)
# assume p1 = 0, phi1 = 1
plogis(opt.age6.sex$par)

# Table
MLE.age6.sex <- array(0, c(7, 2))
MLE.age6.sex[1, ] <- plogis(opt.age6.sex$par)[1]
MLE.age6.sex[2, ] <- 1
MLE.age6.sex[3:7, 1] <- plogis(opt.age6.sex$par)[2:6]
MLE.age6.sex[3:7, 2] <- plogis(opt.age6.sex$par)[7:11]
MLE.age6.sex <- round(data.frame(MLE.age6.sex, row.names = c("p", "phi1", "phi2", "phi3", "phi4", "phi5", "phi6")), 3)
colnames(MLE.age6.sex) <- c("male", "female")

# AIC = -2lik + 2*no.of param
2*opt.age6.sex$value + 2*11   # 11 = dim(theta)

# Save data
AIC.age.sex <- round(data.frame(AIC.age.sex, row.names = c("3 age classes", "4 age classes", "5 age classes")), 3)

saveRDS(AIC.age.sex, file = "data/derived/AIC.age.sex.rds")

