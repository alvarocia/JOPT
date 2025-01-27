#' J-Optimal Subsample Selection
#'
#' This function implements the J-optimal subsample selection method, as described in
#' Cia-Mina et al. (2025). It takes a dataset of covariates \code{x}, a subsample
#' proportion \code{alpha} (between 0 and 1), and a vector defining a regression model.
#' Additional parameters can be specified to control the selection process.
#'
#' @param x A dataset (data frame) containing the covariates for the regression model.
#' @param alpha A numeric value between 0 and 1 specifying the subsample proportion.
#' @param model_vec A character vector defining the regression model. Each element should
#'   represent a term in the model, written as an expression involving \code{x}. For example,
#'   \code{"1"} for the intercept, \code{"x[1]"} for the first covariate, or \code{"x[1]*x[2]^2"}
#'   for an interaction term.
#' @param k0 An integer specifying the initial size of the subsample. Defaults to \code{5*length(model_vec)}.
#' @param q A numeric value between \code{0.5} and \code{1}. Defaults to \code{5/8}.
#' @param gamma A numeric value between \code{0} and \code{q-0.5}. Defaults to \code{1/10}.
#' @param eps1 A small positive value. Defaults to \code{0}.
#'
#' @return A list with the following components:
#' \item{x_j}{A subsample of \code{x} containing the selected observations (rows) according to J-optimality.}
#' \item{idx}{A vector of indices corresponding to the selected rows of \code{x}.}
#'
#' @details
#' The J-optimal subsample selection algorithm selects a subset of observations from the dataset \code{x}
#' that optimizes the statistical efficiency of the model defined by \code{model_vec}. For technical details,
#' refer to Cia-Mina et al. (2025).
#'
#' @importFrom Matrix Diagonal
#' @importFrom latex2exp TeX
#'
#' @examples
#' # Example 1: Bivariate regression
#' set.seed(123)
#' x1 <- runif(1e3, min = -1, max = 1)
#' x2 <- runif(1e3, min = -1, max = 1)
#' x <- data.frame(x1 = x1, x2 = x2)
#' model_vec <- c("1", "x[1]", "x[2]", "x[1]*x[2]", "x[1]^2", "x[2]^2")
#' result <- jseq(x, 0.3, model_vec)
#'
#' # Plot the full dataset and the selected subsample
#' plot(x$x1, x$x2, col = "black", pch = 16, cex = 0.7, xlab = "x1", ylab = "x2")
#' points(result$x_j$x1, result$x_j$x2, col = "red", pch = 16, cex = 0.7)
#' title(main = "J-OPT", line = 1)
#'
#' # Example 2: Univariate regression
#' set.seed(123)
#' x <- data.frame(x = rnorm(1e4))
#' model_vec <- c("1", "x[1]", "x[1]^2")
#' result <- jseq(x, 0.3, model_vec)
#'
#' # Plot the density of the selected subsample
#' plot(density(result$x_j$x, bw = 2 / 100, kernel = "epanechnikov"),
#'      ylab = "", lwd = 1.7, xlim = c(-3.5, 3.5), main = "", xlab = "")
#'
#' @export
jseq <- function(x, alpha, model_vec, k0=5*length(model_vec), q=5/8, gamma=1/10, eps1=0){

n <- nrow(x) # Total sample size (streaming data)
model_function <- create_model_function(model_vec)

#####################################
# STEP 1

my_alpha <- alpha # Subsampling proportion
m <- length(model_vec) # Number of model parameters


#####################################
# STEP 2

# Initial subsample
idx <- 1:k0
x_0_j <- x[1:k0, ,drop=FALSE]

model_matrix <- function(xx){
  return(model_function(xx)%*% t(model_function(xx)) )
}


result_matrix <- model_matrix(as.numeric(x_0_j[1,]))*0
for (i in 1:k0) {
  result_matrix <- result_matrix + model_matrix(as.numeric(x_0_j[i,]))
}
N_0 <- result_matrix/k0
N_0_j_inv <- solve(N_0)

W_app <- N_0

# Criterion
crit_0_j <- -sum(diag(N_0_j_inv %*% W_app)) # Approximated J-opt.

# Directional derivatives
dd_j_opt_app <- numeric(k0)

for (i in 1:k0){
  xx <- as.numeric(x_0_j[i,])
  f <- model_function(xx)
  fTM <- t(f) %*% N_0_j_inv
  dd_j_opt_app[i] <- + fTM %*% W_app %*% t(fTM) + crit_0_j
}

nk_j_opt_app <- k0

k0plus <- ceiling((1-my_alpha/2)*k0)
k0minus <- max(c(floor((1-3*my_alpha/2)*k0),1))

dd_j_opt_app_sorted <- sort(dd_j_opt_app)

# Initialization

# Quantile estimation: c_k_0 hat
c_k_j_opt_app <- dd_j_opt_app_sorted[ceiling((1-my_alpha)*k0)]

beta_0 <- k0/(k0plus-k0minus)

h_j_opt_app <- dd_j_opt_app_sorted[k0plus] - dd_j_opt_app_sorted[k0minus]

h_k_j_opt_app <- h_j_opt_app/(k0^gamma)

f_k_j_opt_app <- sum(abs(dd_j_opt_app_sorted-c_k_j_opt_app)<=h_k_j_opt_app)/(2*k0*h_k_j_opt_app)


#####################################
# STEP 3
# Sequential selection: from k_0+1 to n

pb = txtProgressBar(min = k0+1, max = n, initial = k0+1)

for (i in (k0+1):n){
  # i <- k0+1
  setTxtProgressBar(pb,i)

  f <- model_function(as.numeric(x[i,]))
  N_i <- model_matrix(as.numeric(x[i,]))
  W_app <- (i-1)/i*W_app + N_i/i


  # J-opt approximate
  fTM <- t(f) %*% N_0_j_inv
  zk_j_opt_app <- (fTM %*% W_app %*% t(fTM) + -sum(diag(N_0_j_inv%*% W_app)))[1,1]

  if (zk_j_opt_app>=c_k_j_opt_app){
    nk_j_opt_app <- nk_j_opt_app+1
    N_0_j_inv <- nk_j_opt_app/(nk_j_opt_app-1)*(N_0_j_inv-1/(nk_j_opt_app-1+fTM%*%f)[1,1]*t(fTM)%*%fTM)
    x_0_j <- rbind(x_0_j,x[i,,drop=FALSE])
    idx <- c(idx,i)
  }


  b_k_j_opt_app <- min(c(1/f_k_j_opt_app,beta_0*(i-1)^gamma))

  # Update Ck
  c_k_j_opt_app <- c_k_j_opt_app +b_k_j_opt_app/(i)^q*( (zk_j_opt_app>=c_k_j_opt_app)-my_alpha )

  # Update hk
  h_k_j_opt_app <- h_j_opt_app/(i^gamma)

  # Update fk
  f_k_j_opt_app <- f_k_j_opt_app + ((abs(zk_j_opt_app-c_k_j_opt_app)<=h_k_j_opt_app)/(2*h_k_j_opt_app)-f_k_j_opt_app)/(i^q)

}
result <- list(x_j = x_0_j, idx = idx)
return(result)
}
