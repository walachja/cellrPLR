# Data simulation
gendata1_c <- function(n1=20,n2=20,a=c(10,8,4,2,1,2,1,2,2)*20,peaks=c(2,3),
                       r=c(6,7,4,5,8,9,2,2,3)*100,sigN=0.05,sigS=0.3,sigB=0.001,sigM=0.01,T=0.1){
  library(pls)
  n <- n1+n2
  conc <- a/r
  Tot <- sum(conc)
  p <- length(conc)
  N <- matrix(rnorm(n*p,0,sigN),n,p)
  S <- rnorm(n,0,sigS)%*%t(rep(1,p))
  B <- matrix(rnorm(n*p,0,sigB),n,p)
  R <- rep(1,n)%*%t(r)
  M <- matrix(rnorm(n*p,0,sigM),n,p)
  C <- rep(1,n)%*%t(conc)
  #Closure
  C[1:n1,peaks[1]] <- C[1:n1,peaks[1]]+matrix(1,n1,length(peaks[1]))*T
  C[1:n1,peaks[2]] <- C[1:n1,peaks[2]]+matrix(1,n1,length(peaks[2]))*T*(-1)
  
  X <- N+(1-S)*(C*R+B)*exp(M)
  y <- c(rep(1,n1),rep(-1,n2))
  list(X=X,y=y)
}


gendata2 <- function(n1=20,n2=20,v=500,peaks=1:20, sigN=0.05,sigS=0.3,sigB=0.8,sigM=0.021,T=1.8)
{
  #v...number of variables
  #r <- sample(1:v)*0.02 # In beata's code but in rande from 0.02-10
  r <- runif(v, min = 1, max = 10)
  conc <- runif(v, min = 0, max = 1)*5+matrix(1,v,1)*5
  
  a <- conc*r
  
  n <- n1+n2
  conc <- a/r
  Tot <- sum(conc)
  p <- length(conc)
  N <- matrix(rnorm(n*p,0,sigN),n,p)
  S <- rnorm(n,0,sigS)%*%t(rep(1,p))
  B <- matrix(rnorm(n*p,0,sigB),n,p)
  R <- rep(1,n)%*%t(r)
  M <- matrix(rnorm(n*p,0,sigM),n,p)
  C <- rep(1,n)%*%t(conc)
  
  ###################################
  # Close to closure part
  for (ii in 1:length(peaks))
    C[1:n1,peaks[ii]] <- C[1:n1,peaks[ii]]+matrix(1,n1,length(peaks[ii]))*T*sample(c(-1,1),1)
  X <- N +(1-S)*(C*R+B)*exp(M)
  y <- c(rep(1,n1),rep(-1,n2))
  
  #####################################
  # Getting rid of negative values
  zeros <- which(X<=0,arr.ind = TRUE)
  if (nrow(zeros) > 0)
  {
    for (j in 1:nrow(zeros))
    {
      X[zeros[j,1],zeros[j,2]] <- conc[zeros[j,1]]
    }
  }
  ############################################
  
  colnames(X) <- 1:ncol(X)
  
  list(X=X,y=y,S=S[,1],a=a,r=r,C=C,B=B,N=N,R=R,M=M,g1=which(y==1),g2=which(y==-1))
}