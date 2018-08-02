#' lmrse
#'
#' lmrse fits a linear model with cluster robust standard errors for all markers (e.g. CpG sites of DNA methylation).
#' @docType package
#' @name lmrse
#' @import sandwich
#' @importFrom Rcpp evalCpp
#' @useDynLib lmrse
#' @aliases lmrse lmrse-package
#' @param formula containing the methylation matrix as the response and the exposure and covariates as the dependent terms.
#' @param cluster clustering variable.
#' @param data an optional data.frame which contains the covariates specified in the formula.
#' @return List of coefficients, SE and p-values matrices.
#' @return \item{coef}{a matrix of regression coefficients}
#' @return \item{se}{a matrix of standard errors}
#' @return \item{p}{a matrix of p-values}
#' @examples
#' ### Data
#' y <- rnorm(5000000)
#' y <- matrix(y, ncol=1000)
#' colnames(y) <- paste0("var",1:1000)
#' x <- rnorm(5000)
#' cluster <- rep(1:1000,5)
#' c1 <- rbinom(5000,1,0.5)
#' c2 <- rnorm(5000)
#'
#' ### Analyses
#' res <- lmrse(y ~ x + c1 + c2, cluster=cluster)
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

lmrse <- function(formula, cluster, data=NULL){
  
  # Formulae
  if(!is.null(data)){mf <- model.frame(formula, data, na.action=NULL)}else{mf <- model.frame(formula, na.action=NULL)}
  y <- model.response(mf, "numeric")
  if(!is.null(data)){x <- model.matrix(formula, data)}else{x <- model.matrix(formula, mf)}
  
  # Error messages
  if(nrow(y)!=nrow(x)) stop("the number of rows in the methylation matrix is not equal to the number of rows in the covariates")
  if(nrow(y)!=length(cluster)) stop("the number of rows of the methylation matrix is not equal to the length of the clustering variable")
  
  # Missing covariates
  rm <- apply(is.na(x),1,any)
  if(nrow(x)!=length(rm)) stop("the number of rows of the covairates is not equal to the length of the missing value variable")
  y <- y[!rm,]
  x <- x[!rm,]
  cluster <- cluster[!rm]
  
  # Error messages
  if(nrow(y)!=nrow(x)) stop("the number of rows in the methylation matrix is not equal to the number of rows in the covariates")
  if(nrow(y)!=length(cluster)) stop("the number of rows of the methylation matrix is not equal to the length of the clustering variable")

  # Missing phenotypes
  miss <- apply(is.na(y),2,any)
  y_c <- as.matrix(y[,miss==F])
  y_nc <- as.matrix(y[,miss==T])
  
  # Beta
  if(any(!miss)==T){
    b_c <- t(lm(y_c~x-1)$coef)
  }
  if(any(miss)==T){
    b_nc <- data.frame()
    for(j in 1:ncol(y_nc)){
      b_nc <- rbind(b_nc, t(lm(y_nc[,j]~x-1)$coef))
    }
    b_nc <- as.matrix(b_nc)
  }
  
  # Combine betas
  if(any(miss)==T & any(!miss)==T){
    b <- matrix(NA, nrow=ncol(y), ncol=ncol(x))
    b[miss==F] <- b_c
    b[miss==T] <- b_nc
  }
  if(any(miss)==F & any(!miss)==T){
    b <- b_c
  }
  if(any(miss)==T & any(!miss)==F){
    b <- b_nc
  }
  rownames(b) <- colnames(y)
  colnames(b) <- colnames(x)
  
  # Robust SE
  se <- robustse(y, x, cluster)
  rownames(se) <- colnames(y)
  colnames(se) <- colnames(x)
  
  # P-value
  p <- 2*pnorm(abs(b/se), lower.tail=F)
  
  # Return
  results <- list(coef=b, se=se, p=p)
  class(results) <- "lmrse"
  return(results)
  
}

#' Print lmrse
#'
#' print method for class "lmrse".
#' @param x an object of class "lmrse".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
print.lmrse <- function(x, ...){
  cat("Call: \nlmrse")
  cat("\n\nCoefficients for first marker (",rownames(x[[1]])[1],"):\n", sep="")
  cat(x$coef[1,])
  cat("\n\n")
}

#' Summarizing lmrse
#'
#' summary method for class "lmrse".
#' @param x an object of class "lmrse".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
summary.lmrse <- function(x, ...){
  summ <- x
  class(summ) <- "summary.lmrse"
  return(summ)
}

#' Print summary lmrse
#'
#' print.summary method for class "lmrse".
#' @param x an object of class "lmrse".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
print.summary.lmrse <- function(x, ...){
  cat("Call: lmrse")
  cat("\n\nNumber of markers:", nrow(x[[1]]))
  cat("\n\nCoefficients for the first 10 markers:\n")
  if(nrow(x[[1]])>10){
    cat("\n\nCoefficients for the first 10 markers:\n")
    printCoefmat(x$coef[1:10,])
    cat("\nCluster robust standard errors for the first 10 markers:\n")
    printCoefmat(x$se[1:10,])
    cat("\nP-values for the first 10 markers:\n")
    printCoefmat(x$p[1:10,])
  }else{
    cat("\n\nCoefficients for the markers:\n")
    printCoefmat(x$coef)
    cat("\nCluster robust standard errors for the markers:\n")
    printCoefmat(x$se)
    cat("\nP-values for the markers:\n")
    printCoefmat(x$p)    
  }
}

#' Combine lmrse object into a results data.frame
#'
#' coerce method for class "lmrse", where coefficients, standard errors and p-values are consecutive columns for each covariable in the model.
#' @param x an object of class "lmrse".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
coerce.lmrse <- function(x) {
  coef <- x$coef; se <- x$se; p <- x$p 
  colnames(coef) <- paste0("coef.", sub("\\(Intercept\\)", "intercept", colnames(coef)))
  colnames(se) <- paste0("se.", sub("\\(Intercept\\)", "intercept", colnames(se)))
  colnames(p) <- paste0("pval.", sub("\\(Intercept\\)", "intercept", colnames(p)))
  results <- as.data.frame(cbind(coef, se, p)[, order(c(seq(ncol(coef)), seq(ncol(se)), seq(ncol(p))))])
  results$marker <- rownames(results); rownames(results) <- NULL; results <- results[,c(ncol(results),1:(ncol(results)-1))]
  return(results)
}

#' Robust SE
#'
#' robustse fits cluster robust standard errors (with robusteCpp and/or robustseR) for longitudinal analyses across markers using standard linear regression.
#' @param y matrix of methylation.
#' @param x matrix of other covariates (including intercept).
#' @param cluster clustering variable.
#' @return A matrix of robust standard errors where the rows are the markers (e.g. CpG sites of DNA methylation) and the columns are the covariates including the intercept.
#' @examples
#' ### Data
#' y <- rnorm(5000000)
#' y <- matrix(y, ncol=1000)
#' colnames(y) <- paste0("var",1:1000)
#' x <- rnorm(5000)
#' cluster <- rep(1:1000,5)
#' covar <- cbind(rbinom(5000,1,0.5),rnorm(5000))
#' X <- cbind(1,x,covar)
#' colnames(X) <- c("intercept", "X", "covar1", "covar2")
#'
#' ### Analyses
#' beta <- t((solve(crossprod(X))%*%t(X))%*%y)
#' se <- robustse(y=y, x=X, cluster=cluster)
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

robustse <- function(y=y, x=x, cluster=cluster){
  
  # Cluster variable
  cluster <- as.numeric(factor(cluster)); cluster <- cluster - 1; cluster <- as.integer(cluster)
  
  # Missing phenotypes
  miss <- apply(is.na(y),2,sum)>0
  y_c <- as.matrix(y[,miss==F])
  y_nc <- as.matrix(y[,miss==T])
  
  # Robust SEs
  if(any(!miss)==T){
    robse_c <- robustseCpp(y_c, x, cluster)
  }
  
  # Robust SEs for CpG sites with missing values
  if(any(miss)==T){
    robse_nc <- robustseR(y_nc, x, cluster)
  }
  
  # Robust SEs
  if(any(miss)==T & any(!miss)==T){
    robse <- matrix(NA, nrow=ncol(y), ncol=ncol(x))
    robse[miss==F] <- robse_c
    robse[miss==T] <- robse_nc
  }
  if(any(miss)==F & any(!miss)==T){
    robse <- robse_c
  }
  if(any(miss)==T & any(!miss)==F){
    robse <- robse_nc
  }
  
  # Naming
  colnames(robse) <- NULL
  
  # Return
  return(robse)
  
}

#' Robust SE (C++)
#'
#' robustseCpp fits cluster robust standard errors (using C++Eigen) for longitudinal analyses across CpG sites using standard linear regression.
#' @param y matrix of methylation.
#' @param x matrix of other covariates (including intercept).
#' @param cluster clustering variable.
#' @return A matrix of robust standard errors where the rows are the CpGs and the columns are the covariates including the intercept.
#' @examples
#' ### Data
#' y <- rnorm(5000000)
#' y <- matrix(y, ncol=1000)
#' colnames(y) <- paste0("var",1:1000)
#' x <- rnorm(5000)
#' cluster <- rep(1:1000,5)
#' covar <- cbind(rbinom(5000,1,0.5),rnorm(5000))
#' X <- cbind(1,x,covar)
#' colnames(X) <- c("intercept", "X", "covar1", "covar2")
#'
#' ### Analyses
#' beta <- t((solve(crossprod(X))%*%t(X))%*%y)
#' se <- robustse(y=y, x=X, cluster=cluster)
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

robustseCpp <- function(y=y, x=x, cluster=cluster){
  
  # Cluster variable
  c <- model.matrix(~factor(cluster))[,-1]; colnames(c) <- NULL; c <- cbind(0,c); c[rowSums(c)==0,1] <- 1
  
  # Inverse of xTx
  inv_xtx <- solve(crossprod(x))
  
  # Beta
  b <- ((inv_xtx%*%t(x))%*%y)
  
  # Residuals
  e <- y - x%*%b
  
  # Calculate degree of freedom adjustment
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- ncol(x)
  dfc <- (M/(M-1))*((N-1)/(N-K))
  
  # Robust SEs
  cov <- matrix(0, nrow=ncol(y), ncol=ncol(x))
  robvar <- robustseEigen(e,x,c,inv_xtx,cov)
  robse <- sqrt(dfc*robvar)
  
  # Return
  return(robse)
  
}

#' Robust SE (R)
#'
#' robustseR fits cluster robust standard errors (using R) for longitudinal analyses across CpG sites using standard linear regression.
#' @param y matrix of methylation.
#' @param x matrix of other covariates (including intercept).
#' @param cluster clustering variable.
#' @return A matrix of robust standard errors where the rows are the CpGs and the columns are the covariates including the intercept.
#' @examples
#' ### Data
#' y <- rnorm(5000000)
#' y <- matrix(y, ncol=1000)
#' colnames(y) <- paste0("var",1:1000)
#' x <- rnorm(5000)
#' cluster <- rep(1:1000,5)
#' covar <- cbind(rbinom(5000,1,0.5),rnorm(5000))
#' X <- cbind(1,x,covar)
#' colnames(X) <- c("intercept", "X", "covar1", "covar2")
#'
#' ### Analyses
#' beta <- t((solve(crossprod(X))%*%t(X))%*%y)
#' se <- robustse(y=y, x=X, cluster=cluster)
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

robustseR <- function(y=NULL, x=NULL, cluster=cluster){
  
  # Loop linear regressions and sandwich variance covariances
  n <- ncol(y)
  robse <- data.frame()
  for(j in 1:n){
    d <- na.omit(data.frame(x,y=y[,j],cluster))
    yj <- d$y
    xj <- as.matrix(d[,-c(1,(ncol(d)-1),(ncol(d)))])
    clusterj <- d$cluster
    mod <- lm(yj~xj)
    se <- sand_se(mod, clusterj)
    robse <- rbind(robse,se)
  }
  robse <- as.matrix(robse)
  colnames(robse) <- NULL
  
  # Return
  return(robse=robse)
  
}

#' Sandwich SE
#'
#' sand_se fits cluster robust standard errors using a sandwich estimator.
#' @param model output from linear model.
#' @param cluster clustering variable.
#' @return A vector of robust standard errors for the covariates including the intercept.
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export

sand_se <- function(model, cluster){
  
  # Calculate degree of freedom adjustment
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  
  # Calculate the uj's
  uj  <- apply(estfun(model),2, function(x) tapply(x, cluster, sum))
  
  # Use sandwich to get the covariance matrix
  vcov <- dfc*sandwich(model, meat=crossprod(uj)/N)
  
  # Robust SE
  robse <- sqrt(diag(vcov))
  
  # Return results
  return(robse)
  
}
