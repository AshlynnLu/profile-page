#' Function to transfer capture history to m-array
#'
#' @param ch capture history
#'
#' @return m-array of the input capture history
#'
#' @examples 
#' CH <- as.matrix(CR.sex.chicks09[, -14])
#' marray(CH)
marray <- function(ch){
  nind <- dim(ch)[1]    # no. of individuals 
  n.occasions <- dim(ch)[2]    # no. of occasions
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  
  # calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t, 1] <- sum(ch[, t])
  }
  for (i in 1:nind){
    pos <- which(ch[i, ] != 0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z], pos[z+1]] <- m.array[pos[z], pos[z+1]] + 1
    } #z
  } #i
  
  # calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t, n.occasions+1] <- m.array[t, 1] - sum(m.array[t, 2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1), 2:(n.occasions+1)]
  return(out)
}



#' Function to split CH into different age class
#'
#' @param CH capture history
#' @param age age of a guillemot
#' @param mAge number of age classes, default is set to 3 classes
#'
#' @return A [, , mAge] m-array matrix 
#' The m-array of corresponding age class can be extract by the third index
#'
#' @examples 
#' # 3 age classes
#' CH <- as.matrix(CR.sex.chicks09[, -14])
#' marr1 <- marray.age(CH)[,,1]
#' marr2 <- marray.age(CH)[,,2]
#' marr3 <- marray.age(CH)[,,3]
#' 
#' # 4 age classes
#' marr1 <- marray.age(CH, mAge = 4)[,,1]
#' marr2 <- marray.age(CH, mAge = 4)[,,2]
#' marr3 <- marray.age(CH, mAge = 4)[,,3]
#' marr4 <- marray.age(CH, mAge = 4)[,,4]
marray.age <- function(CH,age,mAge=3){
  
  age <- numeric()
  for (i in 1:dim(CH)[1]){
    age[i] <- 1
  }
  
  # 1.2. Function to remove histories without any capture from a capture-recapture matrix
  clean.CH <- function(CH){
    incl <- which(rowSums(CH)>=1) 
    CH <- CH[incl,]
    return(CH)
  }
  
  # 1.3. Function to remove the first capture in a capture-recapture matrix
  rm.first <- function(CH) {
    get.first <- function(x) min(which(x==1))
    first <- apply(CH, 1, get.first)
    for (i in 1:nrow(CH)){
      CH[i,first[i]] <- 0
    }
    return(CH)
  }
  
  # 1.4. Function to calculate the occasion of first capture
  get.first <- function(x) min(which(x != 0))
  
  
  # 2. Calculations   
  if (is.matrix(CH) == FALSE) CH <- matrix(CH, nrow = 1)   
  maxAge <- max(c(max(age), mAge))
  nind <- nrow(CH)
  n.occasions <- ncol(CH)
  
  first <- apply(CH, 1, get.first)
  age.matrix <- matrix(0, ncol = n.occasions, nrow = nind)
  for (i in 1:nind){
    age.matrix[i,first[i]:n.occasions] <- 1:(n.occasions-first[i]+1)+(age[i]-1)
  }
  age.matrix[age.matrix > maxAge] <- maxAge
  
  # Recode capture history
  CH.rec <- CH
  for (i in 1:nind){
    h <- which(CH.rec[i,]==1)
    for (j in 1:length(h)){
      CH.rec[i,h[j]] <- j
    } # j
  } # i
  CH.rec[CH.rec > maxAge] <- maxAge
  
  CH.split <- array(0, dim = c(nrow(CH), ncol(CH), maxAge))
  for (a in 1:maxAge){
    for (i in 1:nind){
      j <- which(CH.rec[i,]==a | CH.rec[i,]==(a+1))
      if (length(j)==0) next
      CH.split[i,j[1:2],age.matrix[i,j[1]]] <- 1
      if (length(j)>1){
        CH.split[i,j[2:length(j)],age.matrix[i,j[2]]] <- 1
      }
    } # i
  } # a
  
  marr <- array(0, dim = c(n.occasions-1, n.occasions, maxAge))
  for (a in 1:(maxAge-1)){
    for (i in 1:nind){
      u <- which(CH.split[i,,a]==1)
      if (length(u)==0) next
      if (u[1]==n.occasions) next
      if (length(u)==1) marr[u,n.occasions,a] <- marr[u,n.occasions,a] + 1
      if (length(u)==2) marr[u[1],u[2]-1,a] <- marr[u[1],u[2]-1,a] + 1
    } # i
  } # a
  a <- maxAge
  
  if (is.matrix(CH.split[,,a])==FALSE){ 
    CH.split1 <- matrix(CH.split[,,a], nrow = 1)
    marr[,,a] <- marray(CH.split1)
  } # if
  else marr[,,a] <- marray(CH.split[,,a])      
  return(marr)
}

