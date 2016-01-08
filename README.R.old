################################################################################
# File: README.R
# Aim : A brief introduction about the usage for classo and SparCC
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 04/10/2015
#-------------------------------------------------------------------------------
# Package required: 
#                   Matrix for cclasso
#                   gtools for SparCC
# Files needed: 
#               cclasso.R for cclasso (including cclasso)
#               SparCC.R for SparCC (including SparCC.count and SparCC.frac)
# Function parameter description:
#   function: cclasso
#   input:
#          x ------ nxp data matrix, row is sample, col is variable
#     counts ------ x is counts or fraction? default is fraction
#     pseudo ------ pseudo counts if x is counts. default is 0.5
#        sig ------ initial value for covariance matrix. default is NULL     
#       lams ------ tuning parameter sequences. 
#                   default is 10^(seq(0, -8, by = -0.01))
#          K ------ folds of crossvalidation. default is 3     
#       kmax ------ max iteration for augmented lagrangian method. default is 5000
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
#        lam ------ final tuning parameter
#
#   function: SparCC.count
#   input:
#          x ------ nxp count data matrix, row is sample, col is variable
#       imax ------ resampling times from posterior distribution. default 20
#       kmax ------ max iteration steps for SparCC. default is 10
#      alpha ------ the threshold for strong correlation. default is 0.1
#       Vmin ------ minimal variance if negative variance appears. default is 1e-4
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
#
#   function: SparCC.frac
#   input:
#          x ------ nxp fraction data matrix, row is sample, col is variable
#       kmax ------ max iteration steps for SparCC. default is 10
#      alpha ------ the threshold for strong correlation. default is 0.1
#       Vmin ------ minimal variance if negative variance appears. default is 1e-4
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
#-------------------------------------------------------------------------------
# Basic example
# 
source("R/cclasso.R");
source("R/SparCC.R");
# 1. generate logistic normal variables
n <- 100;
p <- 20;
x <- matrix(rnorm(n * p), nrow = n); 
x.frac <- exp(x) / rowSums(exp((x)));
totCount <- round(runif(n = n,  min = 1000, max = 2000));
x.count <- x.frac * totCount;
# 2. run cclasso 
# using fraction
res_ccl_frac <- cclasso(x = x.frac, counts = F);
# using counts
res_ccl_count <- cclasso(x = x.count, counts = T);
# 3. run SparCC.count and SparCC.frac
res_spa_count <- SparCC.count(x = x.count);
res_spa_frac <- SparCC.frac(x = x.frac);
# 4. get the correlation matrix
{
  cat("CCLasso using fraction data:\n");
  print(round(res_ccl_frac$cor.w, 2));
  cat("CCLasso using count data:\n");
  print(round(res_ccl_count$cor.w, 2));
  cat("SparCC using fraction data:\n");
  print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(res_spa_count$cor.w, 2));
}
#-------------------------------------------------------------------------------
