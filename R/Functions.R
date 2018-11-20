cellrPLR_biom <- function(data, type = 'biweight', g1, g2, mainGroup='max',biomarker)
{
  # Main group
  if (mainGroup == 1 | mainGroup == '1'){
    g <- g1
  } else if (mainGroup == 2 | mainGroup == '2'){
    g <- g2
  } else if (mainGroup == 'max'){
       if(length(g1)>length(g2)) {g <- g1
          }  else {g <- g2}
  } else if (mainGroup == 'none'){
    g <- c(g1,g2)
  } else (stop("mainGroup agrument must be: 1,2,'max' or 'none' \n check ?cellrPLR_biom"))
  
  
  # Cube
  vvv1 <- .varmatrixDiagNEW(data,type = type, g=g)
  
  # Aggregation
  vv5 <- apply(vvv1,c(3,1),median,na.rm=TRUE)
  
  # Sorting
  s1 <-  abs((apply(vv5[g1,],2,median,na.rm=TRUE))-(apply(vv5[g2,],2,median,na.rm=TRUE)))
  names(s1) <- colnames(data)
  
  # Ranking
  res <- rep(NA,length(biomarker))
  for (i in 1:length(biomarker))
  {
    res[i] <- which( order(s1,decreasing = TRUE) == which(colnames(data) == biomarker[i]))
  }
  
  output <- list(Difference=s1,Biomarker_results=data.frame(Ordered_position=res,
                                       Biomarker_names=biomarker,
                                       Variable_position=which(colnames(data) %in% biomarker)))
  output$name  <- "Identification of biomarkers"
  class(output) <- "cellrPLR_biom"
  
  return(output)
}



print.cellrPLR_biom <- function (x, ...){
  names <- x$Biomarker_results$Biomarker_names
  df <- data.frame(Biomarker=names,Ordered_Position=x$Biomarker_results$Ordered_position)
  rownames(df) <- NULL
  cat("Positions of biomarkers: ",
      "\n","\n")
  print(df, row.names = FALSE)
}









library(MASS)

# Weights are computed from tau, than data are normalized by weighted mean and computed by Huber,Hampel,Biweight,...
.HampelWNew <- function(x,g, c1 = 4.5, c2 = 3)
{
  n <- length(x)
  medx <- median(x[g])
  x. <- abs(x - medx)
  sigma0 <- median(x.[g])
  mu <- if (c1 > 0) {
    x. <- x./(sigma0 * c1)
    w <- 1 - x. * x.
    w <- ((abs(w) + w)/2)^2
    sum(x * w)/sum(w)
    w
  }
  
  x <- (x-weighted.mean(x[g],w[g]))/mad(x[g])
  w <- psi.hampel(x)
  w <- w-1
  w[x>0] <- -w[x>0]
  return(w)
}

.HuberWNew <- function(x,g, c1 = 4.5, c2 = 3)
{
  n <- length(x)
  medx <- median(x[g])
  x. <- abs(x - medx)
  sigma0 <- median(x.[g])
  mu <- if (c1 > 0) {
    x. <- x./(sigma0 * c1)
    w <- 1 - x. * x.
    w <- ((abs(w) + w)/2)^2
    sum(x * w)/sum(w)
    w
  }
  
  x <- (x-weighted.mean(x[g],w[g]))/mad(x[g])
  w <- psi.huber(x)
  w <- w-1
  w[x>0] <- -w[x>0]
  return(w)
}


.BisquareWNew <- function(x,g, c1 = 4.5, c2 = 3)
{
  n <- length(x)
  medx <- median(x[g])
  x. <- abs(x - medx)
  sigma0 <- median(x.[g])
  mu <- if (c1 > 0) {
    x. <- x./(sigma0 * c1)
    w <- 1 - x. * x.
    w <- ((abs(w) + w)/2)^2
    sum(x * w)/sum(w)
    w
  }
  
  x <- (x-weighted.mean(x[g],w[g]))/mad(x[g])
  w <- psi.bisquare(x)
  w <- w-1
  w[x>0] <- -w[x>0]
  return(w)
}

.psi.tau <- function (x, c1 = 4.5, c2 = 3,g=g, consistency = TRUE, mu.too = FALSE, ...) 
{
  n <- length(x)
  medx <- median(x[g])
  x. <- abs(x - medx)
  sigma0 <- median(x.[g])
  mu <- if (c1 > 0) {
    x. <- x./(sigma0 * c1)
    w <- 1 - x. * x.
    w <- ((abs(w) + w)/2)^2
    sum(x * w)/sum(w)
    w
  }
}



.TauWNew <- function(x,g, c1 = 4.5, c2 = 3)
{
  n <- length(x)
  medx <- median(x[g])
  x. <- abs(x - medx)
  sigma0 <- median(x.[g])
  mu <- if (c1 > 0) {
    x. <- x./(sigma0 * c1)
    w <- 1 - x. * x.
    w <- ((abs(w) + w)/2)^2
    sum(x * w)/sum(w)
    #w
  }
  
  x <- (x-mu)/mad(x[g])#weighted.mean(x[g],w[g])
  w <- psi.tau(x)
  w <- w-1
  w[x>0] <- -w[x>0]
  return(w)
}


.varmatrixDiagNEW <- function (x, type,typeNEW='w', g)
{
  diagnost <- array(NA,dim=c(dim(x)[2],dim(x)[2],dim(x)[1]))
  ss <- NULL
  #pb <- txtProgressBar(min=1,max=ncol(x),style=3)
  #library(robustbase)
  for (i in 1:ncol(x)) {
    for (j in 1:ncol(x)) {
      if (i > j) {
        ss <- .scaleTauDiagNEW(log(x[, i]/x[, j]),type=type, g=g)
        diagnost[i,j,] <- ss
        diagnost[j,i,] <- -ss
      } 
      #setTxtProgressBar(pb,i)
    }
  }
  return(d=diagnost)
}



.scaleTauDiagNEW <- function (x,type,g=g)
{
    if (type == 'tau'){
      w <- .TauWNew(x,g=g)
    } else if (type == 'biweight'){
      w <- .BisquareWNew(x,g=g)
    } else if (type == 'huber'){
      w <- .HuberWNew(x,g=g)
    } else if (type == 'hampel') {
      w <- .HampelWNew(x,g=g)
    }
  return(w)
}
