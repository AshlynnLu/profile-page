## Load function
source("src/helper-functions/functions.R")

## Load data
load(file = "data/raw/CR.sex.chicks09.RData")

CH.m <- as.matrix(CR.sex.chicks09[which(CR.sex.chicks09$Sex == "M"), 1:13])
CH.f <- as.matrix(CR.sex.chicks09[which(CR.sex.chicks09$Sex == "F"), 1:13])
marr.m <- marray(CH.m)
marr.f <- marray(CH.f)

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


# Initialise an array to hold AIC
AIC <- array(0, 5)

#### Base model ####
## Likelihood
caplik.age5.sex1 <- function(theta, data1.m, data2.m, data3.m, data4.m, data5.m, data1.f, data2.f, data3.f, data4.f, data5.f){
  
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
  likhood <- -likhood
}

opt.age5.sex1 <- optim(par = rep(0, 9),
                       fn = caplik.age5.sex1,
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

# p, male(phi 2, 3, 4, 5), female(phi 2, 3, 4, 5), p1 = 0, phi1 = 1
plogis(opt.age5.sex1$par)

# AIC = -2lik + 2*no.of param
AIC[1] <- 2*opt.age5.sex1$value + 2*9   # 9 = dim(theta)



#### phi2 male = female ####
## Likelihood
caplik.age5.sex2 <- function(theta, data1.m, data2.m, data3.m, data4.m, data5.m, data1.f, data2.f, data3.f, data4.f, data5.f){
  
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
    phi2.f[t] <- 1/(1+exp(-theta[2]))
    phi3.f[t] <- 1/(1+exp(-theta[6]))
    phi4.f[t] <- 1/(1+exp(-theta[7]))
    phi5.f[t] <- 1/(1+exp(-theta[8]))
    
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
  likhood <- -likhood
}

opt.age5.sex2 <- optim(par = rep(0, 8),
                       fn = caplik.age5.sex2,
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

# p, phi2, male(phi 3, 4, 5), female(phi 3, 4, 5), p1 = 0, phi1 = 1
plogis(opt.age5.sex2$par)

# AIC = -2lik + 2*no.of param
AIC[2] <- 2*opt.age5.sex2$value + 2*8   # 9 = dim(theta)



#### phi 2, 3 male = female ####
## Likelihood
caplik.age5.sex3 <- function(theta, data1.m, data2.m, data3.m, data4.m, data5.m, data1.f, data2.f, data3.f, data4.f, data5.f){
  
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
    phi2.f[t] <- 1/(1+exp(-theta[2]))
    phi3.f[t] <- 1/(1+exp(-theta[3]))
    phi4.f[t] <- 1/(1+exp(-theta[6]))
    phi5.f[t] <- 1/(1+exp(-theta[7]))
    
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
  likhood <- -likhood
}

opt.age5.sex3 <- optim(par = rep(0, 7),
                       fn = caplik.age5.sex3,
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

# p, phi2, phi3, male(phi4, phi5), female(phi4, phi5), p1 = 0, phi1 = 1
plogis(opt.age5.sex3$par)

# AIC = -2lik + 2*no.of param
AIC[3] <- 2*opt.age5.sex3$value + 2*7   # 7 = dim(theta)



#### phi 2, 3, 4 male = female ####
## Likelihood
caplik.age5.sex4 <- function(theta, data1.m, data2.m, data3.m, data4.m, data5.m, data1.f, data2.f, data3.f, data4.f, data5.f){
  
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
    phi2.f[t] <- 1/(1+exp(-theta[2]))
    phi3.f[t] <- 1/(1+exp(-theta[3]))
    phi4.f[t] <- 1/(1+exp(-theta[4]))
    phi5.f[t] <- 1/(1+exp(-theta[6]))
    
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

opt.age5.sex4 <- optim(par = rep(0, 6),
                       fn = caplik.age5.sex4,
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

# p, phi2, phi3, phi4, male(phi5), female(phi5), p1 = 0, phi1 = 1
plogis(opt.age5.sex4$par)

# AIC = -2lik + 2*no.of param
AIC[4] <- 2*opt.age5.sex4$value + 2*6   # 6 = dim(theta)




#### constant phi's ####
## Likelihood
caplik.age5.sex5 <- function(theta, data1, data2, data3, data4, data5){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data1)[1]
  nj = dim(data1)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi1 <- array(0,nj-1)
  phi2 <- array(0,nj-1)
  phi3 <- array(0,nj-1)
  phi4 <- array(0,nj-1)
  phi5 <- array(0,nj-1)
  
  p1 <- array(0,(nj-1))
  p <- array(0,(nj-1))
  pbar1 <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  
  q1 <- array(0,dim=c(ni,nj))
  q2 <- array(0,dim=c(ni,nj))
  q3 <- array(0,dim=c(ni,nj))
  q4 <- array(0,dim=c(ni,nj))
  q5 <- array(0,dim=c(ni,nj))
  
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p1[t] <- 0
    p[t] <- 1/(1+exp(-theta[1]))
    
    phi1[t] <- 1
    phi2[t] <- 1/(1+exp(-theta[2]))
    phi3[t] <- 1/(1+exp(-theta[3]))
    phi4[t] <- 1/(1+exp(-theta[4]))
    phi5[t] <- 1/(1+exp(-theta[5]))
    
    pbar1[t] <- 1-p1[t]
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q1[t,t] <- phi1[t]*p1[t]
    q2[t,t] <- phi2[t]*p[t]
    q3[t,t] <- phi3[t]*p[t]
    q4[t,t] <- phi4[t]*p[t]
    q5[t,t] <- phi5[t]*p[t]
  }
  
  ## t+1
  ## Ages 1, 2, 3 and 4
  for (t in 1:(ni-1)){
    q1[t,t+1] <- phi1[t]*pbar1[t]*phi2[t+1]*p[t+1]
    q2[t,t+1] <- phi2[t]*pbar[t]*phi3[t+1]*p[t+1]
    q3[t,t+1] <- phi3[t]*pbar[t]*phi4[t+1]*p[t+1]
    q4[t,t+1] <- phi4[t]*pbar[t]*phi5[t+1]*p[t+1]
    
    # Off diagonal elements (only age 5)
    for (j in (t+1):(nj-1)) {
      q5[t,j] <- prod(phi5[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  ## t+2
  for (t in 1:(ni-2)){
    ## Age 4
    for (j in (t+2):(nj-1)) {
      q4[t,j] <- phi4[t]*pbar[t]*prod(phi5[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
    }
    ## Age 1, 2 and 3
    q1[t,t+2] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]
    q2[t,t+2] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*p[t+2]
    q3[t,t+2] <- phi3[t]*pbar[t]*phi4[t+1]*pbar[t+1]*phi5[t+2]*p[t+2]
  }
  
  ## t+3 for Age 1, 2 and 3
  for (t in 1:(ni-3)){
    for (j in (t+3):(nj-1)) {
      q3[t,j] <- phi3[t]*pbar[t]*phi4[t+1]*pbar[t+1]*prod(phi5[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
    }
    q1[t,t+3] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]*phi4[t+3]*p[t+3]
    q2[t,t+3] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*p[t+2]*phi5[t+3]*p[t+3]
  }
  
  ## t+4 for Age 1 and 2
  for (t in 1:(ni-4)){
    for (j in (t+4):(nj-1)) {
      q2[t,j] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*pbar[t+2]*prod(phi5[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
    }
    q1[t,t+4] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*pbar[t+2]*phi4[t+3]*pbar[t+3]*phi5[t+4]*p[t+4]
  }
  
  ## t+5 for Age 1
  for (t in 1:(ni-5)){
    for (j in (t+5):(nj-1)) {
      q1[t,j] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*pbar[t+2]*phi4[t+3]*pbar[t+3]*prod(phi5[(t+4):j])*prod(pbar[(t+4):(j-1)])*p[j]
    }
  }
  
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q1[t,nj] <- 1 - sum(q1[t,t:(nj-1)])
    q2[t,nj] <- 1 - sum(q2[t,t:(nj-1)])
    q3[t,nj] <- 1 - sum(q3[t,t:(nj-1)])
    q4[t,nj] <- 1 - sum(q4[t,t:(nj-1)])
    q5[t,nj] <- 1 - sum(q5[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  
  likhood <- 0
  likhood1 <- 0
  likhood2 <- 0
  likhood3 <- 0
  likhood4 <- 0
  likhood5 <- 0
  
  for (i in 1:ni){
    for (j in i:nj) {
      likhood2 <- likhood2 + data2[i,j]*log(q2[i,j])
      likhood3 <- likhood3 + data3[i,j]*log(q3[i,j])
      likhood4 <- likhood4 + data4[i,j]*log(q4[i,j])
      likhood5 <- likhood5 + data5[i,j]*log(q5[i,j])
    }
    for (j in (i+1):nj){
      likhood1 <- likhood1 + data1[i,j]*log(q1[i,j])
    }
  }
  
  likhood <- likhood1 + likhood2 + likhood3 +likhood4 + likhood5 
  
  # Output the negative loglikelihood value:
  likhood<- -likhood
}

opt.age5.sex5 <- optim(par = rep(0, 5),
                       fn = caplik.age5.sex5,
                       data1 = marr1.age5, 
                       data2 = marr2.age5,
                       data3 = marr3.age5,
                       data4 = marr4.age5,
                       data5 = marr5.age5)

# p, phi2, phi3, phi4, phi5, p1 = 0, phi1 = 1
plogis(opt.age5.sex5$par)

# AIC = -2lik + 2*no.of param
AIC[5] <- 2*opt.age5.sex5$value + 2*5   # 5 = dim(theta)









##### AIC #####
model.selection <- array(0,dim=c(5, 2))

for (i in 1:5){
  model.selection[i, 2] <- round(AIC[i], 3)
}
model.selection[, 1] <- c("1. all phi depends on sex", "2. phi2.m = phi2.f", "3. phi2,3.m = phi2,3.f", "4. phi2,3,4.m = phi2,3,4.f","5. constant phi")

model.selection <- data.frame(model.selection)
colnames(model.selection) <- c("Model", "AIC")

# Save result
saveRDS(model.selection, file = "data/derived/model.selection.rds")




## Best model
opt.age5.sex4.hess <- optimHess(par = rep(0, 6),
                                fn = caplik.age5.sex4,
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

# Lowest AIC: phi2, 3, 4 constant, phi5 depends on sex
params <- array(0, c(6, 3))
params[, 1] <- plogis(opt.age5.sex4$par)

# SD
params[, 2] <- sqrt(diag(solve(opt.age5.sex4.hess)))

params <- round(params, 3)

### resample ###
#### phi1,2,3,4.m = phi1,2,3,4.f ####
## Likelihood
caplik.model4 <- function(theta, data1.m, data2.m, data3.m, data4.m, data5.m, data1.f, data2.f, data3.f, data4.f, data5.f){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data1.m)[1]
  nj = dim(data1.m)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi1 <- array(0,nj-1)
  phi2 <- array(0,nj-1)
  phi3 <- array(0,nj-1)
  phi4 <- array(0,nj-1)
  
  phi5.m <- array(0,nj-1)
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
    
    phi1[t] <- 1
    phi2[t] <- 1/(1+exp(-theta[2]))
    phi3[t] <- 1/(1+exp(-theta[3]))
    phi4[t] <- 1/(1+exp(-theta[4]))
    
    phi5.m[t] <- 1/(1+exp(-theta[5]))
    phi5.f[t] <- 1/(1+exp(-theta[6]))
    
    pbar1[t] <- 1-p1[t]
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q1.m[t,t] <- phi1[t]*p1[t]
    q2.m[t,t] <- phi2[t]*p[t]
    q3.m[t,t] <- phi3[t]*p[t]
    q4.m[t,t] <- phi4[t]*p[t]
    q5.m[t,t] <- phi5.m[t]*p[t]
    
    q1.f[t,t] <- phi1[t]*p1[t]
    q2.f[t,t] <- phi2[t]*p[t]
    q3.f[t,t] <- phi3[t]*p[t]
    q4.f[t,t] <- phi4[t]*p[t]
    q5.f[t,t] <- phi5.f[t]*p[t]
  }
  
  ## t+1
  ## Ages 1, 2, 3 and 4
  for (t in 1:(ni-1)){
    q1.m[t,t+1] <- phi1[t]*pbar1[t]*phi2[t+1]*p[t+1]
    q2.m[t,t+1] <- phi2[t]*pbar[t]*phi3[t+1]*p[t+1]
    q3.m[t,t+1] <- phi3[t]*pbar[t]*phi4[t+1]*p[t+1]
    q4.m[t,t+1] <- phi4[t]*pbar[t]*phi5.m[t+1]*p[t+1]
    
    q1.f[t,t+1] <- phi1[t]*pbar1[t]*phi2[t+1]*p[t+1]
    q2.f[t,t+1] <- phi2[t]*pbar[t]*phi3[t+1]*p[t+1]
    q3.f[t,t+1] <- phi3[t]*pbar[t]*phi4[t+1]*p[t+1]
    q4.f[t,t+1] <- phi4[t]*pbar[t]*phi5.f[t+1]*p[t+1]
    
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
      q4.m[t,j] <- phi4[t]*pbar[t]*prod(phi5.m[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
      q4.f[t,j] <- phi4[t]*pbar[t]*prod(phi5.f[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
    }
    ## Age 1, 2 and 3
    q1.m[t,t+2] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]
    q2.m[t,t+2] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*p[t+2]
    q3.m[t,t+2] <- phi3[t]*pbar[t]*phi4[t+1]*pbar[t+1]*phi5.m[t+2]*p[t+2]
    q1.f[t,t+2] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]
    q2.f[t,t+2] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*p[t+2]
    q3.f[t,t+2] <- phi3[t]*pbar[t]*phi4[t+1]*pbar[t+1]*phi5.f[t+2]*p[t+2]
  }
  
  ## t+3 for Age 1, 2 and 3
  for (t in 1:(ni-3)){
    for (j in (t+3):(nj-1)) {
      q3.m[t,j] <- phi3[t]*pbar[t]*phi4[t+1]*pbar[t+1]*prod(phi5.m[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
      q3.f[t,j] <- phi3[t]*pbar[t]*phi4[t+1]*pbar[t+1]*prod(phi5.f[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
    }
    q1.m[t,t+3] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]*phi4[t+3]*p[t+3]
    q1.f[t,t+3] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]*phi4[t+3]*p[t+3]
    q2.m[t,t+3] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*p[t+2]*phi5.m[t+3]*p[t+3]
    q2.f[t,t+3] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*p[t+2]*phi5.f[t+3]*p[t+3]
  }
  
  ## t+4 for Age 1 and 2
  for (t in 1:(ni-4)){
    for (j in (t+4):(nj-1)) {
      q2.m[t,j] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*pbar[t+2]*prod(phi5.m[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
      q2.f[t,j] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*pbar[t+2]*prod(phi5.f[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
    }
    q1.m[t,t+4] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*pbar[t+2]*phi4[t+3]*pbar[t+3]*phi5.m[t+4]*p[t+4]
    q1.f[t,t+4] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*pbar[t+2]*phi4[t+3]*pbar[t+3]*phi5.f[t+4]*p[t+4]
  }
  
  ## t+5 for Age 1
  for (t in 1:(ni-5)){
    for (j in (t+5):(nj-1)) {
      q1.m[t,j] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*pbar[t+2]*phi4[t+3]*pbar[t+3]*prod(phi5.m[(t+4):j])*prod(pbar[(t+4):(j-1)])*p[j]
      q1.f[t,j] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*pbar[t+2]*phi4[t+3]*pbar[t+3]*prod(phi5.f[(t+4):j])*prod(pbar[(t+4):(j-1)])*p[j]
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

n.resample <- 50
resample.theta.model4 <- matrix(data = 0, ncol = 6, nrow = n.resample)
for (i in 1:n.resample){
  resample <- CR.sex.chicks09[sample(dim(CR.sex.chicks09)[1], size = dim(CR.sex.chicks09)[1], replace = TRUE), ]
  resample.m <- as.matrix(resample[which(resample$Sex == "M"), 1:13])
  resample.f <- as.matrix(resample[which(resample$Sex == "F"), 1:13])
  
  re.marr1.age5.m <- marray.age(resample.m, mAge = 5)[,,1]
  re.marr2.age5.m <- marray.age(resample.m, mAge = 5)[,,2]
  re.marr3.age5.m <- marray.age(resample.m, mAge = 5)[,,3]
  re.marr4.age5.m <- marray.age(resample.m, mAge = 5)[,,4]
  re.marr5.age5.m <- marray.age(resample.m, mAge = 5)[,,5]
  
  re.marr1.age5.f <- marray.age(resample.f, mAge = 5)[,,1]
  re.marr2.age5.f <- marray.age(resample.f, mAge = 5)[,,2]
  re.marr3.age5.f <- marray.age(resample.f, mAge = 5)[,,3]
  re.marr4.age5.f <- marray.age(resample.f, mAge = 5)[,,4]
  re.marr5.age5.f <- marray.age(resample.f, mAge = 5)[,,5]
  
  re.opt.model4 <- optim(par = rep(0, 6),
                         fn = caplik.model4,
                         data1.m = re.marr1.age5.m, 
                         data2.m = re.marr2.age5.m,
                         data3.m = re.marr3.age5.m,
                         data4.m = re.marr4.age5.m,
                         data5.m = re.marr5.age5.m,
                         data1.f = re.marr1.age5.f, 
                         data2.f = re.marr2.age5.f,
                         data3.f = re.marr3.age5.f,
                         data4.f = re.marr4.age5.f, 
                         data5.f = re.marr5.age5.f)
  
  resample.theta.model4[i, ] <- plogis(re.opt.model4$par)
}

# CI
for (i in 1:6){
  params[i, 3] <- paste0("[", round(quantile(resample.theta.model4[, i], 0.025), 3), ",", round(quantile(resample.theta.model4[, i], 0.975), 3), "]")
}

params <- data.frame(params, row.names = c("$p$", "$\\phi_2$", "$\\phi_3$", "$\\phi_4$", "$\\phi_{5,m}$", "$\\phi_{5,f}$"))
colnames(params) <- c("MLE", "SD", "CI")

# Save result
saveRDS(params, file = "data/derived/params.rds")
