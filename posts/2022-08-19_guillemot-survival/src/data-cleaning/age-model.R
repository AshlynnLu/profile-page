## Load function
source("src/helper-functions/functions.R")

## Load data
load(file = "data/raw/CR.sex.chicks09.RData")
CH <- as.matrix(CR.sex.chicks09[, -14])

## 3 age classes
marr1.age3 <- marray.age(CH)[,,1]
marr2.age3 <- marray.age(CH)[,,2]
marr3.age3 <- marray.age(CH)[,,3]

caplik.age3 <- function(theta, data1, data2, data3){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data1)[1]
  nj = dim(data1)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi1 <- array(0,nj-1)
  phi2 <- array(0,nj-1)
  phi3 <- array(0,nj-1)
  
  p1 <- array(0,(nj-1))
  p <- array(0,(nj-1))
  pbar1 <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  
  q1 <- array(0,dim=c(ni,nj))
  q2 <- array(0,dim=c(ni,nj))
  q3 <- array(0,dim=c(ni,nj))
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p1[t] <- 0
    p[t] <- 1/(1+exp(-theta[1]))
    phi1[t] <- 1
    phi2[t] <- 1/(1+exp(-theta[2]))
    phi3[t] <- 1/(1+exp(-theta[3]))
    pbar1[t] <- 1-p1[t]
    pbar[t] <- 1-p[t]
  }
  
  # Calculate the multinomial cell probabilities
  
  # Diagonal elements
  
  for (t in 1:ni){
    q1[t,t] <- phi1[t]*p1[t]
    q2[t,t] <- phi2[t]*p[t]
    q3[t,t] <- phi3[t]*p[t]
  }
  
  ## t+1
  ## Ages 1 and 2
  for (t in 1:(ni-1)){
    q1[t,t+1] <- phi1[t]*pbar1[t]*phi2[t+1]*p[t+1]
    q2[t,t+1] <- phi2[t]*pbar[t]*phi3[t+1]*p[t+1]
    
    # Off diagonal elements (only age 3)
    for (j in (t+1):(nj-1)) {
      q3[t,j] <- prod(phi3[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  ## t+2
  ## Age 2
  for (t in 1:(ni-2)){
    ## Age 2
    for (j in (t+2):(nj-1)) {
      q2[t,j] <- phi2[t]*pbar[t]*prod(phi3[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
    }
    ## Age 1
    q1[t,t+2] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]
  }
  
  ## t+3 for Age 1
  for (t in 1:(ni-3)){
    for (j in (t+3):(nj-1)) {
      q1[t,j] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*prod(phi3[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
    }
  }
  
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q1[t,nj] <- 1 - sum(q1[t,t:(nj-1)])
    q2[t,nj] <- 1 - sum(q2[t,t:(nj-1)])
    q3[t,nj] <- 1 - sum(q3[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  
  likhood <- 0
  likhood1 <- 0
  likhood2 <- 0
  likhood3 <- 0
  
  for (i in 1:ni){
    for (j in i:nj) {
      #likhood2 <- likhood2 + data2[i,j]*log(q2[i,j])
      likhood3 <- likhood3 + data3[i,j]*log(q3[i,j])
    }
    for (j in (i+1):nj){
      likhood1 <- likhood1 + data1[i,j]*log(q1[i,j])
    }
  }
  
  likhood <- likhood1 + likhood2 + likhood3
  
  # Output the negative loglikelihood value:
  likhood<- -likhood
}

opt.age3 <- optim(par = rep(0, 3),
                  fn = caplik.age3,
                  data1 = marr1.age3, 
                  data2 = marr2.age3,
                  data3 = marr3.age3)

opt.age3.hess <- optimHess(par = rep(0, 3),
                           fn = caplik.age3,
                           data1 = marr1.age3, 
                           data2 = marr2.age3,
                           data3 = marr3.age3)

# AIC
AIC.age3 <- 2*opt.age3$value + 2*3   # 3 = dim(theta)

# MLE
estimate3 <- round(plogis(opt.age3$par), 3)
MLE.age3 <- array("-", c(6, 3))
MLE.age3[2:4, 1] <- estimate3[1:3]


# SD
sd3 <- round(sqrt(diag(solve(opt.age3.hess))), 3)
MLE.age3[2:4, 2] <- sd3[1:3]

# CI
n.resample <- 50
resample.theta.age3 <- matrix(data = 0, ncol = 3, nrow = n.resample)
for (i in 1:n.resample){
  resample <- as.matrix(CR.sex.chicks09[sample(dim(CR.sex.chicks09)[1], size = dim(CR.sex.chicks09)[1], replace = TRUE), -14])

  re.marr1.age3 <- marray.age(resample, mAge = 3)[,,1]
  re.marr2.age3 <- marray.age(resample, mAge = 3)[,,2]
  re.marr3.age3 <- marray.age(resample, mAge = 3)[,,3]
  
  opt.resample <- optim(par = rep(0, 3),
                        fn = caplik.age3,
                        data1 = re.marr1.age3, 
                        data2 = re.marr2.age3,
                        data3 = re.marr3.age3)
  
  resample.theta.age3[i, ] <- plogis(opt.resample$par)
}
for (i in 2:4){
  MLE.age3[i, 3] <- paste0("[", round(quantile(resample.theta.age3[, (i-1)], 0.025), 3), ",", round(quantile(resample.theta.age3[, (i-1)], 0.975), 3), "]")
}




## 4 age classes
marr1.age4 <- marray.age(CH, mAge = 4)[,,1]
marr2.age4 <- marray.age(CH, mAge = 4)[,,2]
marr3.age4 <- marray.age(CH, mAge = 4)[,,3]
marr4.age4 <- marray.age(CH, mAge = 4)[,,4]

## Likelihood
caplik.age4 <- function(theta, data1, data2, data3, data4){
  
  # Number of release occasions and recovery occasions
  
  ni = dim(data1)[1]
  nj = dim(data1)[2]
  
  # Set up the size of the arrays containing the survival probabilities, recapture 
  # probabilities and cell probabilities and set them all initially to zero 
  
  phi1 <- array(0,nj-1)
  phi2 <- array(0,nj-1)
  phi3 <- array(0,nj-1)
  phi4 <- array(0,nj-1)
  
  p1 <- array(0,(nj-1))
  p <- array(0,(nj-1))
  pbar1 <- array(0,(nj-1))
  pbar <- array(0,(nj-1))
  
  q1 <- array(0,dim=c(ni,nj))
  q2 <- array(0,dim=c(ni,nj))
  q3 <- array(0,dim=c(ni,nj))
  q4 <- array(0,dim=c(ni,nj))
  
  # Define the parameters to be constant or time-dependent
  
  for (t in 1:nj) {
    p1[t] <- 0
    p[t] <- 1/(1+exp(-theta[1]))
    phi1[t] <- 1
    phi2[t] <- 1/(1+exp(-theta[2]))
    phi3[t] <- 1/(1+exp(-theta[3]))
    phi4[t] <- 1/(1+exp(-theta[4]))
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
  }
  
  ## t+1
  ## Ages 1, 2 and 3
  for (t in 1:(ni-1)){
    q1[t,t+1] <- phi1[t]*pbar1[t]*phi2[t+1]*p[t+1]
    q2[t,t+1] <- phi2[t]*pbar[t]*phi3[t+1]*p[t+1]
    q3[t,t+1] <- phi3[t]*pbar[t]*phi4[t+1]*p[t+1]
    
    # Off diagonal elements (only age 4)
    for (j in (t+1):(nj-1)) {
      q4[t,j] <- prod(phi4[t:j])*prod(pbar[t:(j-1)])*p[j]
    }
  }
  
  ## t+2
  for (t in 1:(ni-2)){
    ## Age 3
    for (j in (t+2):(nj-1)) {
      q3[t,j] <- phi3[t]*pbar[t]*prod(phi4[(t+1):j])*prod(pbar[(t+1):(j-1)])*p[j]
    }
    ## Age 1 and 2
    q1[t,t+2] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]
    q2[t,t+2] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*phi4[t+2]*p[t+2]
  }
  
  ## t+3 for Age 1 and 2
  for (t in 1:(ni-3)){
    for (j in (t+3):(nj-1)) {
      q2[t,j] <- phi2[t]*pbar[t]*phi3[t+1]*pbar[t+1]*prod(phi4[(t+2):j])*prod(pbar[(t+2):(j-1)])*p[j]
    }
    q1[t,t+3] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*p[t+2]*phi4[t+3]*p[t+3]
  }
  
  ## t+4 for Age 1
  for (t in 1:(ni-4)){
    for (j in (t+4):(nj-1)) {
      q1[t,j] <- phi1[t]*pbar1[t]*phi2[t+1]*pbar[t+1]*phi3[t+2]*pbar[t+2]*prod(phi4[(t+3):j])*prod(pbar[(t+3):(j-1)])*p[j]
    }
  }
  # Calculate the disappearing animal probabilities
  
  for (t in 1:ni){
    q1[t,nj] <- 1 - sum(q1[t,t:(nj-1)])
    q2[t,nj] <- 1 - sum(q2[t,t:(nj-1)])
    q3[t,nj] <- 1 - sum(q3[t,t:(nj-1)])
    q4[t,nj] <- 1 - sum(q4[t,t:(nj-1)])
  }
  
  # Calculate the likelihood function
  
  likhood <- 0
  likhood1 <- 0
  likhood2 <- 0
  likhood3 <- 0
  likhood4 <- 0
  
  for (i in 1:ni){
    for (j in i:nj) {
      #likhood2 <- likhood2 + data2[i,j]*log(q2[i,j])
      likhood3 <- likhood3 + data3[i,j]*log(q3[i,j])
      likhood4 <- likhood4 + data4[i,j]*log(q4[i,j])
    }
    for (j in (i+1):nj){
      likhood1 <- likhood1 + data1[i,j]*log(q1[i,j])
    }
  }
  
  likhood <- likhood1 + likhood2 + likhood3 + likhood4
  
  # Output the negative loglikelihood value:
  likhood <- -likhood
}

opt.age4 <- optim(par = rep(0, 4),
                  fn = caplik.age4,
                  data1 = marr1.age4, 
                  data2 = marr2.age4,
                  data3 = marr3.age4,
                  data4 = marr4.age4)

opt.age4.hess <- optimHess(par = rep(0, 4),
                           fn = caplik.age4,
                           data1 = marr1.age4, 
                           data2 = marr2.age4,
                           data3 = marr3.age4,
                           data4 = marr4.age4)

# AIC
AIC.age4 <- 2*opt.age3$value + 2*4   # 4 = dim(theta)

# MLE
estimate4 <- round(plogis(opt.age4$par), 3)
MLE.age4 <- array("-", c(6, 3))
MLE.age4[2:5, 1] <- estimate4[1:4]

# SD
sd4 <- round(sqrt(diag(solve(opt.age4.hess))), 3)
MLE.age4[2:5, 2] <- sd4[1:4]

# CI
resample.theta.age4 <- matrix(data = 0, ncol = 4, nrow = n.resample)
for (i in 1:n.resample){
  resample <- as.matrix(CR.sex.chicks09[sample(dim(CR.sex.chicks09)[1], size = dim(CR.sex.chicks09)[1], replace = TRUE), -14])

  re.marr1.age4 <- marray.age(resample, mAge = 4)[,,1]
  re.marr2.age4 <- marray.age(resample, mAge = 4)[,,2]
  re.marr3.age4 <- marray.age(resample, mAge = 4)[,,3]
  re.marr4.age4 <- marray.age(resample, mAge = 4)[,,4]
  
  opt.resample <- optim(par = rep(0, 4),
                        fn = caplik.age4,
                        data1 = re.marr1.age4, 
                        data2 = re.marr2.age4,
                        data3 = re.marr3.age4, 
                        data4 = re.marr4.age4)
  
  resample.theta.age4[i, ] <- plogis(opt.resample$par)
}

for (i in 2:5){
  MLE.age4[i, 3] <- paste0("[", round(quantile(resample.theta.age4[, (i-1)], 0.025), 3), ",", round(quantile(resample.theta.age4[, (i-1)], 0.975), 3), "]")
}


## 5 age classes
marr1.age5 <- marray.age(CH, mAge = 5)[,,1]
marr2.age5 <- marray.age(CH, mAge = 5)[,,2]
marr3.age5 <- marray.age(CH, mAge = 5)[,,3]
marr4.age5 <- marray.age(CH, mAge = 5)[,,4]
marr5.age5 <- marray.age(CH, mAge = 5)[,,5]

## Likelihood
caplik.age5 <- function(theta, data1, data2, data3, data4, data5){
  
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
  likhood <- -likhood
}

opt.age5 <- optim(par = rep(0, 5),
                  fn = caplik.age5,
                  data1 = marr1.age5, 
                  data2 = marr2.age5,
                  data3 = marr3.age5,
                  data4 = marr4.age5,
                  data5 = marr5.age5)

opt.age5.hess <- optimHess(par = rep(0, 5),
                           fn = caplik.age5,
                           data1 = marr1.age5, 
                           data2 = marr2.age5,
                           data3 = marr3.age5,
                           data4 = marr4.age5,
                           data5 = marr5.age5)
# AIC
AIC.age5 <- 2*opt.age3$value + 2*5   # 5 = dim(theta)

# MLE
estimate5 <- round(plogis(opt.age5$par), 3)
MLE.age5 <- array("-", c(6, 3))
MLE.age5[2:6, 1] <- estimate5[1:5]

# SD
sd5 <- round(sqrt(diag(solve(opt.age5.hess))), 3)
MLE.age5[2:6, 2] <- sd5[1:5]

# CI
resample.theta.age5 <- matrix(data = 0, ncol = 5, nrow = n.resample)
for (i in 1:n.resample){
  resample <- as.matrix(CR.sex.chicks09[sample(dim(CR.sex.chicks09)[1], size = dim(CR.sex.chicks09)[1], replace = TRUE), -14])

  re.marr1.age5 <- marray.age(resample, mAge = 5)[,,1]
  re.marr2.age5 <- marray.age(resample, mAge = 5)[,,2]
  re.marr3.age5 <- marray.age(resample, mAge = 5)[,,3]
  re.marr4.age5 <- marray.age(resample, mAge = 5)[,,4]
  re.marr5.age5 <- marray.age(resample, mAge = 5)[,,5]
  
  opt.resample <- optim(par = rep(0, 5),
                        fn = caplik.age5,
                        data1 = re.marr1.age5, 
                        data2 = re.marr2.age5,
                        data3 = re.marr3.age5, 
                        data4 = re.marr4.age5,
                        data5 = re.marr5.age5)
  
  resample.theta.age5[i, ] <- plogis(opt.resample$par)
}

for (i in 2:6){
  MLE.age5[i, 3] <- paste0("[", round(quantile(resample.theta.age5[, (i-1)], 0.025), 3), ",", round(quantile(resample.theta.age5[, (i-1)], 0.975), 3), "]")
}


# table
MLE.age <- array(0, c(6, 9))
MLE.age[ , 1:3] <- MLE.age3
MLE.age[ , 4:6] <- MLE.age4
MLE.age[ , 7:9] <- MLE.age5
MLE.age[1, ] <- c("MLE", "SD", "CI", "MLE", "SD", "CI", "MLE", "SD", "CI")
MLE.age <- data.frame(MLE.age, row.names = c(" ", "$p$", "$\\phi_2$", "$\\phi_3$", "$\\phi_4$", "$\\phi_5$"))
colnames(MLE.age) <- c( "3 age classes"," ", " ", "4 age classes"," ", " ", "5 age classes", " ",  " ")

saveRDS(MLE.age, file = "data/MLE.age.rds")


# AIC
AIC.age <- array(0, c(3, 1))
AIC.age[1, 1] <- AIC.age3
AIC.age[2, 1] <- AIC.age4
AIC.age[3, 1] <- AIC.age5
AIC.age <- data.frame(AIC.age)
rownames(AIC.age) <- c("3 age classes", "4 age classes", "5 age classes")
colnames(AIC.age) <- "AIC"

saveRDS(AIC.age, file = "data/AIC.age.rds")

