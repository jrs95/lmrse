# lmrse <img src='man/figures/logo.png' align="right" height="139"/>

This package is used to fit linear models with cluster robust standard errors across high-dimensional phenotypes to assess change over time.  

## Installation
```
install.packages("remotes")
remotes::install_github("jrs95/lmrse")
```

## Functions
* `lmrse`: fits a linear model with cluster robust standard errors for all markers.  

## Example
```
## Data  
y <- rnorm(5000000)
y <- matrix(y, ncol = 1000) # a matrix of phenotypes (rows = individuals, columns = markers)
colnames(y) <- paste0("pheno", 1:1000)
x <- rnorm(5000) # a vector of exposure
cluster <- rep(1:1000, 5) # cluster variable
c1 <- rbinom(5000, 1, 0.5) # covariate 1
c2 <- rnorm(5000) # covariate 2

## Analyses  
res <- lmrse(y ~ x + c1 + c2, cluster = cluster)
summary(res)
results <- coerce.lmrse(res)
```

## Citation
Staley JR *et al.* Longitudinal analysis strategies for modelling epigenetic trajectories. [Int J Epidemiol](https://pubmed.ncbi.nlm.nih.gov/29462323/) 2018;47(2):516-525.
