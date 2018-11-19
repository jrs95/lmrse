# lmrse
This package is used to fit linear models with cluster robust standard errors across high-dimensional phenotypes (e.g. DNA methylation at CpG sites) to assess change over time. 

# Functions
* lmrse - fits a linear model with cluster robust standard errors for all markers (e.g. CpG sites of DNA methylation). 

# Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("jrs95/lmrse")
4. library(lmrse)

# Example
\#\#\# Data  
y <- rnorm(5000000)  
y <- matrix(y, ncol=1000) # a matrix of phenotypes with rows of individuals and columns of phenotypes (e.g. rows of samples and coloumns of CpG sites)   
colnames(y) <- paste0("var",1:1000)  
x <- rnorm(5000) # a vector of exposure   
cluster <- rep(1:1000,5) # cluster variable   
c1 <- rbinom(5000,1,0.5) # covariate 1  
c2 <- rnorm(5000) # covariate 2    

\#\#\# Analyses  
res <- lmrse(y ~ x + c1 + c2, cluster=cluster)  
summary(res)  
results <- coerce.lmrse(res)
