##install.packages("expm")
require(expm)
# library(ape)
# mu_G <- vcv(rtree(2))
alpha_G <- 0.2
sigma_G <- diag(0.1,2)
mu_G <- matrix(c(1,0.8,0.8,1),nrow = 2)
delta_t <- 0.0001
G_0 <-  diag(rep(1,2))
N <- 10000 #no of points

AI_sim <- function (alpha_G,sigma_G, mu_G, G_0, delta_t, N) {
  ans <- array(NA, dim = c(2,2,N))
  ans[,,1] <- G_0
  
  for(i in 2:N) {
    W_t <- rnorm(3,0,1)
    W_t <- matrix(c(W_t[1],W_t[3]/sqrt(2),W_t[3]/sqrt(2),W_t[2]), nrow = 2)
    sqrt_G <- sqrtm(ans[,,i-1])
    ans[,,i] <- sqrt_G %*% expm(alpha_G * solve(sqrt_G) %*% mu_G %*% 
                                 solve(sqrt_G) * delta_t +
                                 sigma_G  %*% (sqrt(delta_t) * W_t)) %*% sqrt_G
    
  }
  return(ans)
}
X <- AI_sim(alpha_G,sigma_G,mu_G,G_0,delta_t,N)
x1 <- X[1,1,1:10000]
x2 <- X[1,2,1:10000]
plot(x1)
plot(x2)
plot(X)

