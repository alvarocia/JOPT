# Title: Efficiency Comparison of J-optimal and D-optimal Subdata Selection
#
# Description:
# This script compares the efficiency of two subdata selection methods:
# J-optimal and D-optimal. It evaluates their performance based on
# predefined efficiency criteria and generates comparative plots
# to visualize the results.
# Theoretical J-optimal is included for comparison.

#' Run Efficiency Comparison Example
#'
#' This function provides an example of running an efficiency comparison
#' using the package's capabilities. The function contains the example
#' code previously executed at the top level.
#'
#' @import Matrix
#' @import latex2exp
#'
#' @examples
#' # To run the example:
#' run_efficiency_comparison_example()
#'
#' @export
run_efficiency_comparison_example <- function() {
# library(Matrix)
# library(latex2exp)

my_seed <- 1234
set.seed(my_seed)

n <- 1e5 # Total sample size (streaming data)

x <- rnorm(n)
x2 <- x^2


# Theoretical expression of matrix W. For comparison
# Non-centered moments for standard normal distribution:


mux <- 0
mux2 <- 1
mux3 <- 0
mux4 <- 3
W <- Matrix(c(1, mux, mux2, mux ,mux2,mux3, mux2, mux3,mux4), nrow = 3, ncol = 3, byrow = TRUE)


df_crit <- data.frame()


#####################################
# STEP 1

# Parameters of the algorithm
my_alpha <- 1e-1 # Subsampling proportion
m <- 3 # Number of model parameters
k0 <- 5*m # Initial subsample size
q <- 5/8
gamma <- 1/10
eps1 <- 0


#####################################
# STEP 2

# Initial subsample
x_0_th <- x[1:k0]
x_0_d <- x_0_th
x_0_j <- x_0_th


xx <- mean(x_0_j)
xx2 <- mean(x_0_j^2)
xx3 <- mean(x_0_j^3)
xx4 <- mean(x_0_j^4)
N_0 <- Matrix(c(1, xx, xx2, xx,xx2,xx3,xx2,xx3,xx4), nrow = 3, ncol = 3, byrow = TRUE)


N_0_th_inv <- solve(N_0)
N_0_d_inv <- N_0_th_inv
N_0_j_inv <- N_0_th_inv

# Approximated expression for W. Updated in each step of the algorithm
xx_app <- xx
xx2_app <- xx2
xx3_app <- xx3
xx4_app <- xx4
W_app <- Matrix(c(1, xx_app, xx2_app, xx_app,xx2_app,xx3_app,xx2_app,xx3_app,xx4_app), nrow = 3, ncol = 3, byrow = TRUE)


# Criteria (maximize)
crit_0_th <- -sum(diag(N_0_th_inv %*% W)) # Theoretical J-opt. For comparison
crit_0_d <- -log(det(N_0_d_inv)) # D-opt
crit_0_j <- -sum(diag(N_0_j_inv %*% W_app)) # Approximated J-opt. PROPOSED

# Directional derivatives
dd_j_opt_th <- numeric(k0)
dd_j_opt_app <- numeric(k0)
dd_d_opt <- numeric(k0)

for (i in 1:k0){
  xx <- x_0_th[i]
  f <- Matrix(c(1,xx,xx^2),nrow=3,ncol=1, byrow=TRUE)
  fTM <- t(f) %*% N_0_th_inv
  dd_d_opt[i] <- - m + fTM %*% f
  dd_j_opt_th[i] <- + fTM %*% W %*% t(fTM) + crit_0_th
  dd_j_opt_app[i] <- + fTM %*% W_app %*% t(fTM) + crit_0_j
}

nk_j_opt_th <- k0
nk_j_opt_app <- k0
nk_d_opt <- k0

k0plus <- ceiling((1-my_alpha/2)*k0)
k0minus <- max(c(floor((1-3*my_alpha/2)*k0),1))

dd_j_opt_th_sorted <- sort(dd_j_opt_th)
dd_j_opt_app_sorted <- sort(dd_j_opt_app)
dd_d_opt_sorted <- sort(dd_d_opt)

# Initialization

# Quantile estimation: c_k_0 hat
c_k_j_opt_th <- dd_j_opt_th_sorted[ceiling((1-my_alpha)*k0)]
c_k_j_opt_app <- dd_j_opt_app_sorted[ceiling((1-my_alpha)*k0)]
c_k_d_opt <- dd_d_opt_sorted[ceiling((1-my_alpha)*k0)]

beta_0 <- k0/(k0plus-k0minus)

h_j_opt_th <- dd_j_opt_th_sorted[k0plus] - dd_j_opt_th_sorted[k0minus]
h_j_opt_app <- dd_j_opt_app_sorted[k0plus] - dd_j_opt_app_sorted[k0minus]
h_d_opt <- dd_d_opt_sorted[k0plus] - dd_d_opt_sorted[k0minus]

h_k_j_opt_th <- h_j_opt_th/(k0^gamma)
h_k_j_opt_app <- h_j_opt_app/(k0^gamma)
h_k_d_opt <- h_d_opt/(k0^gamma)

f_k_j_opt_th <- sum(abs(dd_j_opt_th_sorted-c_k_j_opt_th)<=h_k_j_opt_th)/(2*k0*h_k_j_opt_th)
f_k_j_opt_app <- sum(abs(dd_j_opt_app_sorted-c_k_j_opt_app)<=h_k_j_opt_app)/(2*k0*h_k_j_opt_app)
f_k_d_opt <- sum(abs(dd_d_opt_sorted-c_k_d_opt)<=h_k_d_opt)/(2*k0*h_k_d_opt)


#####################################
# STEP 3
# Sequential selection: from k_0+1 to n

pb = txtProgressBar(min = k0+1, max = n, initial = k0+1)

for (i in (k0+1):n){

  setTxtProgressBar(pb,i)

  f <- Matrix(c(1,x[i],x[i]^2),nrow=3,ncol=1, byrow=TRUE)
  N_i <- f%*%t(f)
  W_app <- (i-1)/i*W_app + N_i/i


  # Theoretical J-opt
  fTM <- t(f) %*% N_0_th_inv
  zk_j_opt <- (fTM %*% W %*% t(fTM)  -sum(diag(N_0_th_inv %*% W)))[1,1]

  if (zk_j_opt>=c_k_j_opt_th){
    nk_j_opt_th <- nk_j_opt_th+1
    N_0_th_inv <- nk_j_opt_th/(nk_j_opt_th-1)*(N_0_th_inv-1/(nk_j_opt_th-1+fTM%*%f)[1,1]*t(fTM)%*%fTM)
    x_0_th <- c(x_0_th,x[i])
  }

  # Proposed J-opt approximate
  fTM <- t(f) %*% N_0_j_inv
  zk_j_opt_app <- (fTM %*% W_app %*% t(fTM) + -sum(diag(N_0_j_inv%*% W_app)))[1,1]

  if (zk_j_opt_app>=c_k_j_opt_app){
    nk_j_opt_app <- nk_j_opt_app+1
    N_0_j_inv <- nk_j_opt_app/(nk_j_opt_app-1)*(N_0_j_inv-1/(nk_j_opt_app-1+fTM%*%f)[1,1]*t(fTM)%*%fTM)
    x_0_j <- c(x_0_j,x[i])
  }

  # D-opt
  fTM <- t(f) %*% N_0_d_inv
  zk_d_opt <- (-m + fTM %*% f )[1,1]

  if (zk_d_opt>=c_k_d_opt){
    nk_d_opt <- nk_d_opt+1
    N_0_d_inv <- nk_d_opt/(nk_d_opt-1)*(N_0_d_inv-1/(nk_d_opt-1+fTM%*%f)[1,1]*t(fTM)%*%fTM)
    x_0_d <- c(x_0_d,x[i])
  }


  b_k_j_opt <- min(c(1/f_k_j_opt_th,beta_0*(i-1)^gamma))
  b_k_j_opt_app <- min(c(1/f_k_j_opt_app,beta_0*(i-1)^gamma))
  b_k_d_opt <- min(c(1/f_k_d_opt,beta_0*(i-1)^gamma))

  # Update Ck
  c_k_j_opt_th <- c_k_j_opt_th +b_k_j_opt/(i)^q*( (zk_j_opt>=c_k_j_opt_th)-my_alpha )
  c_k_j_opt_app <- c_k_j_opt_app +b_k_j_opt_app/(i)^q*( (zk_j_opt_app>=c_k_j_opt_app)-my_alpha )
  c_k_d_opt <- c_k_d_opt +b_k_d_opt/(i)^q*( (zk_d_opt>=c_k_d_opt)-my_alpha )

  # Update hk
  h_k_j_opt_th <- h_j_opt_th/(i^gamma)
  h_k_j_opt_app <- h_j_opt_app/(i^gamma)
  h_k_d_opt <- h_d_opt/(i^gamma)

  # Update fk
  f_k_j_opt_th <- f_k_j_opt_th + ((abs(zk_j_opt-c_k_j_opt_th)<=h_k_j_opt_th)/(2*h_k_j_opt_th)-f_k_j_opt_th)/(i^q)
  f_k_j_opt_app <- f_k_j_opt_app + ((abs(zk_j_opt_app-c_k_j_opt_app)<=h_k_j_opt_app)/(2*h_k_j_opt_app)-f_k_j_opt_app)/(i^q)
  f_k_d_opt <- f_k_d_opt + ((abs(zk_d_opt-c_k_d_opt)<=h_k_d_opt)/(2*h_k_d_opt)-f_k_d_opt)/(i^q)


  df_crit <- rbind(df_crit,data.frame(i=i,j_opt=sum(diag(N_0_th_inv %*% W)),
                                      d_opt=sum(diag(N_0_d_inv %*% W)),
                                      j_opt_app=sum(diag(N_0_j_inv %*% W))))

}



#####################################
# PLOT

# dir.create("PLOT", showWarnings = FALSE)
# pdf(file = "PLOT/histogram.pdf", width = 8, height = 4)

par(mfrow=c(1,2),mai=c(.7,.7,.2,.2))

plot(density(x_0_d,bw=1/100,kernel="epanechnikov"),ylab="",lwd=1.7,xlim=c(-4,4), main="",xlab="")
title(xlab = latex2exp::TeX(r'($x$)'), mgp = c(2, 1, 0))
title(ylab = latex2exp::TeX(r'($h^*_m$)'), mgp=c(2,1,0))

plot(density(x_0_j,bw=1/100,kernel="epanechnikov"),ylab="",lwd=1.7,xlim=c(-4,4), main="",xlab="")
title(xlab = latex2exp::TeX(r'($x$)'), mgp = c(2, 1, 0))
title(ylab = latex2exp::TeX(r'($h^*_m$)'), mgp=c(2,1,0))

# dev.off()





# pdf(file = "PLOT/efficiency.pdf", width = 7, height = 5)
par(mfrow=c(1,1))
plot(log10(df_crit$i), df_crit$j_opt, col="blue", pch="",xlab="", ylab="",ylim=c(0.2,1.1) )
title(xlab=latex2exp::TeX(r'($log_{10}k$)'), mgp = c(2, 1, 0))
title(ylab = latex2exp::TeX(r'($eff_{\Phi_J}$)'), mgp=c(2,1,0))
abline(h=1, col = adjustcolor("black", alpha = 0.4), lty="dotted")
abline(h=0.8635933945394123, col = adjustcolor("black", alpha = 0.4), lty="dotted")
lines(log10(df_crit$i), 1.800994153417123/df_crit$j_opt, col="chartreuse3", lty=1,lwd=1.5)
lines(log10(df_crit$i), 1.800994153417123/df_crit$j_opt_app, col="coral2",lty=4, ylab="J-Eff",lwd=1.5)
lines(log10(df_crit$i), 1.800994153417123/df_crit$d_opt, col="blue", lty=3,lwd=1.5)
legend("bottomright",legend=c("D-OPT",latex2exp::TeX(r'($J_k-OPT$)'),"J-OPT"), col=c("blue","coral2","chartreuse3"),
       lty=c(3,4,1), lwd=c(1.5,1.5,1.5), ncol=1)
# dev.off()

}
