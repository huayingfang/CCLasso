################################################################################
# File: README.R
# Aim : A brief introduction about the usage for classo and SparCC
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 2015-01-08
#-------------------------------------------------------------------------------
# Package required: 
#                   gtools for SparCC
# Files needed: 
#               cclasso.R for cclasso (including cclasso)
#               SparCC.R for SparCC (including SparCC.count and SparCC.frac)
#-------------------------------------------------------------------------------
# Function parameter description:
# function: cclasso
#   Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#                    n samples & p compositional variables
#      counts ------ Is the compositional data matrix a count matrix? 
#                    Default: FALSE
#      pseudo ------ pseudo count if counts = TRUE
#                    Default: 0.5
#        k_cv ------ folds of cross validation
#                    Default: 3     
#     lam_int ------ tuning parameter interval
#                    Default: [1e-4, 1]
#       k_max ------ maximum iterations for golden section method
#                    Default: 20
#      n_boot ------ Bootstrap times
#                    Default: 20
#   Output: 
#      A list structure contains:
#       var_w ------ variance estimation
#       cor_w ------ correlation estimation
#      p_vals ------ p-values for elements of cor_w equal 0 or not
#      lambda ------ final tuning parameter
#     info_cv ------ information for cross validation
#-------------------------------------------------------------------------------
# function: SparCC.count
#   input:
#          x ------ nxp count data matrix, row is sample, col is variable
#       imax ------ resampling times from posterior distribution 
#                   default: 20
#       kmax ------ max iteration steps for SparCC 
#                   default: 10
#      alpha ------ the threshold for strong correlation
#                   default: 0.1
#       Vmin ------ minimal variance if negative variance appears
#                   default: 1e-4
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
#
# function: SparCC.frac
#   input:
#          x ------ nxp fraction data matrix, row is sample, col is variable
#       kmax ------ max iteration steps for SparCC
#                   default: 10
#      alpha ------ the threshold for strong correlation
#                   default: 0.1
#       Vmin ------ minimal variance if negative variance appears
#                   default: 1e-4
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
#-------------------------------------------------------------------------------
# Basic example
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
  print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(res_ccl_count$cor_w, 2));
  cat("SparCC using fraction data:\n");
  print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(res_spa_count$cor.w, 2));
}
#-------------------------------------------------------------------------------
