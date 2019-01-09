cellmap <- function(data, g1, g2, mainGroup='max',simul)
{
  library(reshape2)
  library(ggplot2)
  library(plotly)
  
  if (!is.data.frame(data) & !is.matrix(x)) {
    stop("Wrong x data format. Use data.frame or matrix")
  }
  if (min(data)<=0) {
    stop("Only positive values of data can be used")
  }
  
  n <- nrow(data)
  p <- ncol(data)
  n1=length(g1)
  n2=length(g2)
  
  if (max(c(g1,g2))>n) {
    stop("g1 or g2 bigger than dimensionality of data")
  }
  if (n1+n2<n) {
    warning("Not all samples from data are being used")
  }
  if (n1+n2>n) {
    warning("Some samples are selected in both groups")
  }
  
  
  # Main group
  if (mainGroup == 1 | mainGroup == '1'){
    g <- g1; group <- 1
  } else if (mainGroup == 2 | mainGroup == '2'){
    g <- g2; group <- 1
  } else if (mainGroup == 'max'){
    if(length(g1)>length(g2)) {g <- g1; group <- 1
    }  else {g <- g2; group <- 1}
  } else if (mainGroup == 'all'){
    g <- c(g1,g2); group <- 'all'
  } else (stop("mainGroup agrument must be: 1,2,'max' or 'all' \n check ?cellrPLR_biom"))
  
  # Cube
  if (simul == TRUE) {
    vvv1 <- .varmatrixDiagNEW_cell(data, g=g)
    vv1 <- .aggregation(vvv1[[1]],median)
    vv2 <- .aggregation(vvv1[[1]],mean)
    vv3 <- .aggregation(vvv1[[2]],median)
    vv4 <- .aggregation(vvv1[[2]],mean)
    vv5 <- .aggregation(vvv1[[3]],median)
    vv6 <- .aggregation(vvv1[[3]],mean)
  } else {
    vvv1 <- .varmatrixDiagNEW(data,type = 'biweight', g=g)
    vv1 <- .aggregation(vvv1,median)
    vv2 <- .aggregation(vvv1,mean)
    
    vvv1 <- .varmatrixDiagNEW(data,type = 'huber', g=g)
    vv3 <- .aggregation(vvv1,median)
    vv4 <- .aggregation(vvv1,mean)
    
    vvv1 <- .varmatrixDiagNEW(data,type = 'hampel', g=g)
    vv5 <- .aggregation(vvv1,median)
    vv6 <- .aggregation(vvv1,mean)
    
  }
  
  return(list(vv1=vv1,
              vv2=vv2,
              vv3=vv3,
              vv4=vv4,
              vv5=vv5,
              vv6=vv6))
}

plot_cellheatmap <- function(data, type = 'biweight', g1, g2, mainGroup='max', plotly=FALSE, grid=FALSE, title=NULL)
{
  library(reshape2)
  library(ggplot2)
  library(plotly)
  
  if (!is.data.frame(data) & !is.matrix(x)) {
    stop("Wrong x data format. Use data.frame or matrix")
  }
  if (min(data)<=0) {
    stop("Only positive values of data can be used")
  }
  
  n <- nrow(data)
  p <- ncol(data)
  n1=length(g1)
  n2=length(g2)
  
  if (max(c(g1,g2))>n) {
    stop("g1 or g2 bigger than dimensionality of data")
  }
  if (n1+n2<n) {
    warning("Not all samples from data are being used")
  }
  if (n1+n2>n) {
    warning("Some samples are selected in both groups")
  }
  
  
  # Main group
  if (mainGroup == 1 | mainGroup == '1'){
    g <- g1; group <- 1
  } else if (mainGroup == 2 | mainGroup == '2'){
    g <- g2; group <- 1
  } else if (mainGroup == 'max'){
    if(length(g1)>length(g2)) {g <- g1; group <- 1
    }  else {g <- g2; group <- 1}
  } else if (mainGroup == 'all'){
    g <- c(g1,g2); group <- 'all'
  } else (stop("mainGroup agrument must be: 1,2,'max' or 'all' \n check ?cellrPLR_biom"))

  # Cube
  vvv1 <- .varmatrixDiagNEW(data,type = type, g=g)
  
  # Aggregation
  vv5 <- t(.aggregation(vvv1))
  rownames(vv5) <- colnames(data)
  colnames(vv5) <- rownames(data)
  
  v <- reshape2::melt(vv5)
  v$Var1 <- as.factor(v$Var1)
  v$Var1 <- factor(v$Var1,levels = levels(v$Var1))
  v$Var2 <- as.factor(v$Var2)
  v$Var2 <- factor(v$Var2,levels = levels(v$Var2))
  
  
  grid <- ifelse(grid==TRUE,'gray','white')
  p1 <- ggplot(v, aes(Var1, Var2)) + 
        geom_tile(aes(fill = value),color=grid) +
        scale_fill_gradient2(low = "blue",mid='white',high = "red",midpoint = 0,limits=c(-1,1)) +
        theme_grey(base_size = 15)  + scale_y_discrete(expand = c(0, 0))+
        scale_x_discrete(expand = c(0, 0))+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.3), axis.text.x = element_text(angle = 90, hjust = 1))+
        xlab('Variables') + ylab('Samples') + ggtitle(title)
 if(plotly == FALSE){
   p1
 } else print(ggplotly(p1))
  
  return(t(vv5))
}


plot.cellrPLR_biom <- function(x, plotly=FALSE, grid=FALSE, title=NULL)
{
  library(reshape2)
  library(ggplot2)
  library(plotly)
  
  vv5 <- x$Cell_outliers
  
  # Aggregation
  vv5 <- t(vv5)
  rownames(vv5) <- colnames(data)
  colnames(vv5) <- rownames(data)
  
  v <- reshape2::melt(vv5)
  v$Var1 <- as.factor(v$Var1)
  v$Var1 <- factor(v$Var1,levels = levels(v$Var1))
  v$Var2 <- as.factor(v$Var2)
  v$Var2 <- factor(v$Var2,levels = levels(v$Var2))
  
  
  grid <- ifelse(grid==TRUE,'gray','white')
  p1 <- ggplot(v, aes(Var1, Var2)) + 
    geom_tile(aes(fill = value),color=grid) +
    scale_fill_gradient2(low = "blue",mid='white',high = "red",midpoint = 0,limits=c(-1,1)) +
    theme_grey(base_size = 15)  + scale_y_discrete(expand = c(0, 0))+
    scale_x_discrete(expand = c(0, 0))+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.3), axis.text.x = element_text(angle = 90, hjust = 1))+
    xlab('Variables') + ylab('Samples') + ggtitle(title)
  if(plotly == FALSE){
    p1
  } else print(ggplotly(p1))
  
  #return(t(vv5))
}



cellrPLR_biom <- function(data, type = 'biweight', g1, g2, mainGroup='max', biomarker, permutation=FALSE, B=1000, p.alpha=0.95)
{
  if (!is.data.frame(data) & !is.matrix(data)) {
    stop("Wrong x data format. Use data.frame or matrix")
  }
  if (min(data)<=0) {
    stop("Only positive values of data can be used")
  }
  
  n <- nrow(data)
  p <- ncol(data)
  n1=length(g1)
  n2=length(g2)
  
  if (max(c(g1,g2))>n) {
    stop("g1 or g2 bigger than dimensionality of data")
  }
  if (n1+n2<n) {
    warning("Not all samples from data are being used")
  }
  if (n1+n2>n) {
    warning("Some samples are selected in both groups")
  }
  
  
  # Main group
  if (mainGroup == 1 | mainGroup == '1'){
    g <- g1
  } else if (mainGroup == 2 | mainGroup == '2'){
    g <- g2
  } else if (mainGroup == 'max'){
       if(length(g1)>length(g2)) {g <- g1
          }  else {g <- g2}
  } else if (mainGroup == 'all'){
    g <- c(g1,g2)
  } else (stop("mainGroup agrument must be: 1,2,'max' or 'all' \n check ?cellrPLR_biom"))
  
  doBiomarker <- !missing(biomarker)
  
  # Cube
  vvv1 <- .varmatrixDiagNEW(data,type = type, g=g)
  
  # Aggregation
  vv5 <- .aggregation(vvv1)
  
  # Sorting + ranking
  s1 <-  .sorting(vv5, g1, g2, colnames(data))
  
  output <- list(Difference=s1)
  output$Ordered_list <- rank(-output$Difference)
  
  
  # Ranking
  if (!missing(biomarker))
    {
      res <- .ranking(s1, biomarker, colnames(data))
      output$Biomarker_results <- data.frame(Ordered_position=res,
                                             Biomarker_names=biomarker,
                                             Variable_position=which(colnames(data) %in% biomarker))
   }

    # Permutation tests
  if (permutation == TRUE)
    {
      perm <- .permutation(data, vvv1, vv5, s1, g1, g2, B=B, alpha = p.alpha, doBiomarker = doBiomarker, biomarker = biomarker)
      output$Permutation_tests <- perm
    }
  
  output$Cell_outliers <- vv5
  output$name  <- "cell_rPLR"
  class(output) <- "cellrPLR_biom"
  
  return(output)
}



print.cellrPLR_biom <- function (x, ...){
  names <- x$Biomarker_results$Biomarker_names
  df <- data.frame(Biomarker=names,Ordered_Position=x$Biomarker_results$Ordered_position)
  rownames(df) <- NULL
  
  if(ncol(df)>0 | nrow(df)>0)
    {
      cat("Positions of biomarkers: ")
      print(df, row.names = FALSE)
      cat("\n","\n")
  }
  
  names <- x$Permutation_tests$biomarkers_names
  df <- data.frame(Biomarker=names)
  rownames(df) <- NULL
  
  if(ncol(df)>0 | nrow(df)>0)
  {
    cat("Permutation tests - selected variables: ")
    print(df, row.names = FALSE)
  }
  
}



.aggregation <- function(vvv1,type=median)
{
  return(apply(vvv1,c(3,1),type,na.rm=TRUE))
}

.sorting <- function(vv5,g1,g2,names=colnames(data))
{
  s1 <- abs((apply(vv5[g1,],2,median,na.rm=TRUE))-(apply(vv5[g2,],2,median,na.rm=TRUE)))
  names(s1) <- colnames(data)
  return(s1)
}

.ranking <- function(s1,biomarker,names=colnames(data))
{
res <- rep(NA,length(biomarker))
for (i in 1:length(biomarker))
  {
    res[i] <- which( order(s1,decreasing = TRUE) == which(names == biomarker[i]))
  }
  return(res)
}

.permutation <- function(data,vvv1, vv5, s1, g1, g2, B, alpha=0.95, doBiomarker , biomarker)
{
  r <- rep(0,ncol(data))
 
  vvP1 <- array(NA,dim(vvv1))
  
  for (i in 1:B)
    {
      vvP1[,,g1] <- vvv1[sample(1:dim(vvv1)[1]),,g1]
      vvP1[,,g2] <- vvv1[sample(1:dim(vvv1)[1]),,g2]
    
      vvP <- apply(vvP1,c(3,1),median,na.rm=TRUE)
      s1_New <-  abs((apply(vvP[g1,],2,median))-(apply(vvP[g2,],2,median)))
    
      r[s1>s1_New] <- r[s1>s1_New] + 1
    }
  r <- r/B 
  names(r) <- colnames(data)
  b <- which(r > (alpha))
  
  if (doBiomarker == TRUE)
  {
    biom <- which(colnames(data) %in% biomarker)
    
    TP <- sum(biom %in% b)
    TPR <- TP / length(biom)
    
    FP <- (length(b)-TP) 
    FPR <- FP/ ncol(data)
    FDR <- FP/ length(biom)
    
    statistics <- c(TP, TPR, FP, FPR, FDR)
  }
  names(statistics)= c('TP','TPR','FP','FPR','FDR')

  output <- list(p_value = r,TRUE_FALSE = r > (alpha), biomarkers=b, biomarkers_names = colnames(data)[b],statistics= statistics)
  return(output)
  
}



library(MASS)

# Weights are computed from tau, than data are normalized by weighted mean and computed by Huber,Hampel,Biweight,...
.HampelWNew <- function(x,g, c1 = 4.5, c2 = 3)
{
  library(MASS)
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
  library(MASS)
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
  library(MASS)
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


.varmatrixDiagNEW <- function (x, type, g)
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

.varmatrixDiagNEW_cell <- function (x, g)
{
  diagnost1 <- array(NA,dim=c(dim(x)[2],dim(x)[2],dim(x)[1]))
  diagnost2 <- array(NA,dim=c(dim(x)[2],dim(x)[2],dim(x)[1]))
  diagnost3 <- array(NA,dim=c(dim(x)[2],dim(x)[2],dim(x)[1]))
  
  ss1 <- NULL
  ss2 <- NULL
  ss3 <- NULL
  
  #pb <- txtProgressBar(min=1,max=ncol(x),style=3)
  #library(robustbase)
  for (i in 1:ncol(x)) {
    for (j in 1:ncol(x)) {
      if (i > j) {
        lr <- log(x[, i]/x[, j])
        ss1 <- .scaleTauDiagNEW(lr,type='biweight', g=g)
        ss2 <- .scaleTauDiagNEW(lr,type='huber', g=g)
        ss3 <- .scaleTauDiagNEW(lr,type='hampel', g=g)
        diagnost1[i,j,] <- ss1
        diagnost1[j,i,] <- -ss1
        
        diagnost2[i,j,] <- ss2
        diagnost2[j,i,] <- -ss2
        
        diagnost3[i,j,] <- ss3
        diagnost3[j,i,] <- -ss3
      } 
      #setTxtProgressBar(pb,i)
    }
  }
  return(list(biweight=diagnost1,huber=diagnost2,hampel=diagnost3))
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
