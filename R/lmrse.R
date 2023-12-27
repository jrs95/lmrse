#' @title Linear regression model with cluster robust standard errors
#'
#' @description `lmrse` fits a linear regression model with cluster robust
#'   standard errors for all markers.
#'
#' @docType package
#'
#' @name lmrse
#'
#' @import sandwich
#'
#' @importFrom Rcpp evalCpp
#'
#' @useDynLib lmrse
#'
#' @aliases lmrse lmrse-package
#'
#' @param formula containing the marker matrix as the response and the exposure
#'   and covariates as the dependent terms
#'
#' @param cluster clustering variable
#'
#' @param data an optional data.frame which contains the covariates specified
#'   in the formula
#'
#' @return List of coefficients, SE and p-values matrices
#'
#' @return \item{coef}{a matrix of regression coefficients}
#'
#' @return \item{se}{a matrix of standard errors}
#'
#' @return \item{p}{a matrix of p-values}
#'
#' @examples
#' ## Data
#' y <- rnorm(5000000)
#' y <- matrix(y, ncol = 1000)
#' colnames(y) <- paste0("var", 1:1000)
#' x <- rnorm(5000)
#' cluster <- rep(1:1000, 5)
#' c1 <- rbinom(5000, 1, 0.5)
#' c2 <- rnorm(5000)
#'
#' ## Analyses
#' res <- lmrse(y ~ x + c1 + c2, cluster = cluster)
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
lmrse <- function(formula, cluster, data = NULL) {

  # Formulae
  if (!is.null(data)) {
    mf <- model.frame(formula, data, na.action = NULL)
  } else {
    mf <- model.frame(formula, na.action = NULL)
  }
  y <- model.response(mf, "numeric")
  x <- model.matrix(formula, mf)
  names_y <- colnames(y)
  names_x <- colnames(x)

  # Error messages
  if (class(y)[1] != "matrix" || ncol(y) < 2) {
    stop("the outcome should be a marker matrix with at least two columns")
  }
  if (class(x)[1] != "matrix" || ncol(x) < 2) {
    stop(
      "the covariates matrix should contain at least an intercept and an ",
      "exposure"
    )
  }
  if (!all(x[, 1] == 1)) {
    stop("the first column of the covariates matrix should be the intercept")
  }
  if (nrow(y) != nrow(x)) {
    stop(
      "the number of rows in the marker matrix is not equal to the number of ",
      "rows in the covariates"
    )
  }
  if (nrow(y) != length(cluster)) {
    stop(
      "the number of rows of the marker matrix is not equal to the length of ",
      "the clustering variable"
    )
  }
  if (any(is.na(cluster))) {
    stop("there are missing values in the clustering variable")
  }

  # Missing covariates
  rm <- apply(is.na(x), 1, any)
  if (nrow(x) != length(rm)) {
    stop(
      "the number of rows of the covariates is not equal to the length of the ",
      "missing value variable"
    )
  }
  y <- y[!rm, , drop = FALSE]
  x <- x[!rm, , drop = FALSE]
  cluster <- cluster[!rm]

  # Error messages
  if (nrow(y) != nrow(x)) {
    stop(
      "the number of rows in the marker matrix is not equal to the number of ",
      "rows in the covariates"
    )
  }
  if (nrow(y) != length(cluster)) {
    stop(
      "the number of rows of the marker matrix is not equal to the length of ",
      "the clustering variable"
    )
  }

  # Missing phenotypes
  miss <- apply(is.na(y), 2, any)
  if (ncol(y) != length(miss)) {
    stop(
      "the number of columns of the marker matrix is not equal to the length ",
      "of the missing value variable"
    )
  }
  y_c <- y[, miss == FALSE, drop = FALSE]
  y_nc <- y[, miss == TRUE, drop = FALSE]

  # Beta
  if (any(!miss) == TRUE) {
    b_c <- t(lm(y_c ~ x - 1)$coef)
    if (any(is.na(b_c))) {
      stop("there are missing values in the regression coefficients")
    }
  }
  if (any(miss) == TRUE) {
    b_nc <- data.frame()
    for (j in seq_len(ncol(y_nc))) {
      b_nc <- rbind(b_nc, t(lm(y_nc[, j] ~ x - 1)$coef))
    }
    b_nc <- as.matrix(b_nc)
    if (any(is.na(b_nc))) {
      stop("there are missing values in the regression coefficients")
    }
  }

  # Combine betas
  if (any(miss) == TRUE && any(!miss) == TRUE) {
    b <- matrix(NA, nrow = ncol(y), ncol = ncol(x))
    b[miss == FALSE, ] <- b_c
    b[miss == TRUE, ] <- b_nc
  }
  if (any(miss) == FALSE && any(!miss) == TRUE) {
    b <- b_c
  }
  if (any(miss) == TRUE && any(!miss) == FALSE) {
    b <- b_nc
  }
  rownames(b) <- names_y
  colnames(b) <- names_x

  # Robust SE
  se <- robustse(y, x, cluster)
  rownames(se) <- names_y
  colnames(se) <- names_x

  # P-value
  p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)

  # Return
  results <- list(coef = b, se = se, p = p)
  class(results) <- "lmrse"
  return(results)

}

#' @title Print lmrse
#'
#' @description print method for class `"lmrse"`.
#'
#' @name print.lmrse
#'
#' @param x an object of class `"lmrse"`
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
print.lmrse <- function(x, ...) {
  cat("Call: \nlmrse")
  cat(
    "\n\nCoefficients for first marker (",
    rownames(x[[1]])[1],
    "):\n",
    sep = ""
  )
  cat(x$coef[1, ])
  cat("\n\n")
}

#' @title Summary of lmrse
#'
#' @description summary method for class `"lmrse"`.
#'
#' @name summary.lmrse
#'
#' @param x an object of class `"lmrse"`
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
summary.lmrse <- function(x, ...) {
  summ <- x
  class(summ) <- "summary.lmrse"
  return(summ)
}

#' @title Print summary lmrse
#'
#' @description print.summary method for class `"lmrse"`.
#'
#' @name print.summary.lmrse
#'
#' @param x an object of class `"lmrse"`
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
print.summary.lmrse <- function(x, ...) {
  cat("Call: lmrse")
  cat("\n\nNumber of markers:", nrow(x[[1]]))
  cat("\n\nCoefficients for the first 10 markers:\n")
  if (nrow(x[[1]]) > 10) {
    cat("\n\nCoefficients for the first 10 markers:\n")
    printCoefmat(x$coef[1:10, , drop = FALSE])
    cat("\nCluster robust standard errors for the first 10 markers:\n")
    printCoefmat(x$se[1:10, , drop = FALSE])
    cat("\nP-values for the first 10 markers:\n")
    printCoefmat(x$p[1:10, , drop = FALSE])
  } else {
    cat("\n\nCoefficients for the markers:\n")
    printCoefmat(x$coef)
    cat("\nCluster robust standard errors for the markers:\n")
    printCoefmat(x$se)
    cat("\nP-values for the markers:\n")
    printCoefmat(x$p)
  }
}

#' @title Combine lmrse object into a results data.frame
#'
#' @description coerce method for class `"lmrse"`, where
#'   coefficients, standard errors and p-values are
#'   consecutive columns for each covariable in the model.
#'
#' @name coerce.lmrse
#'
#' @param x an object of class "lmrse"
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
coerce.lmrse <- function(x) {
  coef <- x$coef
  se <- x$se
  p <- x$p
  colnames(coef) <- paste0(
    "coef.", sub("\\(Intercept\\)", "intercept", colnames(coef))
  )
  colnames(se) <- paste0(
    "se.", sub("\\(Intercept\\)", "intercept", colnames(se))
  )
  colnames(p) <- paste0(
    "pval.", sub("\\(Intercept\\)", "intercept", colnames(p))
  )
  results <- as.data.frame(cbind(coef, se, p))
  results <- results[
    , order(c(seq_len(ncol(coef)), seq_len(ncol(se)), seq_len(ncol(p)))),
    drop = FALSE
  ]
  results$marker <- rownames(results)
  rownames(results) <- NULL
  results <- results[, c(ncol(results), 1:(ncol(results) - 1)), drop = FALSE]
  return(results)
}

#' @title Robust SE
#'
#' @description `robustse` fits cluster robust standard errors
#'   (with robusteCpp and/or robustseR) for longitudinal analyses
#'   across markers using standard linear regression.
#'
#' @name robustse
#'
#' @param y matrix of markers
#'
#' @param x matrix of other covariates (including intercept)
#'
#' @param cluster clustering variable
#'
#' @return `robustse` returns a matrix of robust standard errors where
#'   the rows are the markers (e.g. CpG sites of DNA methylation)
#'   and the columns are the covariates including the intercept.
#'
#' @examples
#' ## Data
#' y <- rnorm(5000000)
#' y <- matrix(y, ncol = 1000)
#' colnames(y) <- paste0("var", 1:1000)
#' x <- rnorm(5000)
#' cluster <- rep(1:1000, 5)
#' covar <- cbind(rbinom(5000, 1, 0.5), rnorm(5000))
#' X <- cbind(1, x, covar)
#' colnames(X) <- c("intercept", "X", "covar1", "covar2")
#'
#' ## Analyses
#' beta <- t((solve(crossprod(X)) %*% t(X)) %*% y)
#' se <- robustse(y = y, x = X, cluster = cluster)
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
robustse <- function(y, x, cluster) {

  # Cluster variable
  cluster <- as.numeric(factor(cluster))
  cluster <- cluster - 1
  cluster <- as.integer(cluster)

  # Missing phenotypes
  miss <- apply(is.na(y), 2, any)
  if (ncol(y) != length(miss)) {
    stop(
      "the number of columns of the marker matrix is not equal to the length ",
      "of the missing value variable"
    )
  }
  y_c <- y[, miss == FALSE, drop = FALSE]
  y_nc <- y[, miss == TRUE, drop = FALSE]

  # Robust SEs
  if (any(!miss) == TRUE) {
    robse_c <- robustseCpp(y_c, x, cluster)
  }

  # Robust SEs for CpG sites with missing values
  if (any(miss) == TRUE) {
    robse_nc <- robustseR(y_nc, x, cluster)
  }

  # Robust SEs
  if (any(miss) == TRUE && any(!miss) == TRUE) {
    robse <- matrix(NA, nrow = ncol(y), ncol = ncol(x))
    robse[miss == FALSE, ] <- robse_c
    robse[miss == TRUE, ] <- robse_nc
  }
  if (any(miss) == FALSE && any(!miss) == TRUE) {
    robse <- robse_c
  }
  if (any(miss) == TRUE && any(!miss) == FALSE) {
    robse <- robse_nc
  }

  # Naming
  colnames(robse) <- NULL

  # Return
  return(robse)

}

#' @title Robust SE (C++)
#'
#' @description `robustseCpp` fits cluster robust standard errors
#'   (using C++Eigen) for longitudinal analyses across CpG sites
#'   using standard linear regression.
#'
#' @name robustseCpp
#'
#' @param y matrix of markers
#'
#' @param x matrix of other covariates (including intercept)
#'
#' @param cluster clustering variable
#'
#' @return `robustseCpp` returns a matrix of robust standard errors
#'   where the rows are the CpGs and the columns are the covariates
#'   including the intercept.
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
robustseCpp <- function(y, x, cluster) {

  # Cluster variable
  cc <- model.matrix(~ factor(cluster))[, -1]
  colnames(cc) <- NULL
  cc <- cbind(0, cc)
  cc[rowSums(cc) == 0, 1] <- 1

  # Inverse of xTx
  inv_xtx <- solve(crossprod(x))

  # Beta
  b <- ((inv_xtx %*% t(x)) %*% y)

  # Residuals
  e <- y - x %*% b

  # Calculate degree of freedom adjustment
  m <- length(unique(cluster))
  n <- length(cluster)
  k <- ncol(x)
  dfc <- (m / (m - 1)) * ((n - 1) / (n - k))

  # Robust SEs
  cov <- matrix(0, nrow = ncol(y), ncol = ncol(x))
  robvar <- robustseEigen(e, x, cc, inv_xtx, cov)
  robse <- sqrt(dfc * robvar)

  # Return
  return(robse)

}

#' @title Robust SE (R)
#'
#' @description `robustseR` fits cluster robust standard errors
#'   (using R) for longitudinal analyses across CpG sites using
#'   standard linear regression.
#'
#' @name robustseR
#'
#' @param y matrix of markers
#'
#' @param x matrix of other covariates (including intercept)
#'
#' @param cluster clustering variable
#'
#' @return `robustseR` returns a matrix of robust standard errors
#'   where the rows are the CpGs and the columns are the covariates
#'   including the intercept.
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
robustseR <- function(y, x, cluster) {

  # Linear regression model and sandwich variance covariances for each marker
  n <- ncol(y)
  robse <- data.frame()
  for (j in seq_len(n)) {
    yj <- y[, j]
    miss <- is.na(yj)
    if (length(yj) != length(miss)) {
      stop(
        "the number of rows in the marker matrix is not equal to the ",
        "length of the missing variable"
      )
    }
    if (nrow(x) != length(miss)) {
      stop(
        "the number of rows in the covariate matrix is not equal to the ",
        "length of the missing variable"
      )
    }
    if (length(cluster) != length(miss)) {
      stop(
        "the cluster variable is not equal to the length of the missing ",
        "variable"
      )
    }
    yj <- yj[!miss]
    xj <- x[!miss, , drop = FALSE]
    clusterj <- cluster[!miss]
    mod <- lm(yj ~ xj - 1)
    se <- sandwich.se(mod, clusterj)
    robse <- rbind(robse, se)
  }
  robse <- as.matrix(robse)
  colnames(robse) <- NULL

  # Return
  return(robse)

}

#' @title Sandwich SE
#'
#' @description `sandwich.se` fits cluster robust standard errors using a
#'   sandwich estimator.
#'
#' @name sandwich.se
#'
#' @import sandwich
#'
#' @param model output from linear model
#'
#' @param cluster clustering variable
#'
#' @return `sandwich.se` returns a vector of robust standard errors for the
#'   covariates including the intercept.
#'
#' @examples
#' ## Data
#' y <- rnorm(5000)
#' x <- rnorm(5000)
#' cluster <- rep(1:1000, 5)
#' c1 <- rbinom(5000, 1, 0.5)
#' c2 <- rnorm(5000)
#'
#' ## Analyses
#' model <- lm(y ~ x + c1 + c2)
#' se <- sandwich.se(model = model, cluster = cluster)
#'
#' @author James R Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
sandwich.se <- function(model, cluster) {

  # Calculate degree of freedom adjustment
  m <- length(unique(cluster))
  n <- length(cluster)
  k <- model$rank
  dfc <- (m / (m - 1)) * ((n - 1) / (n - k))

  # Calculate the uj's
  uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum))

  # Use sandwich to get the covariance matrix
  vcov <- dfc * sandwich(model, meat = crossprod(uj) / n)

  # Robust SE
  robse <- sqrt(diag(vcov))

  # Return results
  return(robse)

}
