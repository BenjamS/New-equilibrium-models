n <- 3
Cmat <- matrix(runif(n * n), n, n)
Ci <- rowSums(Cmat)
R <- Ci * (1 + runif(n))
lam <- rep(1, n)
#lam <- c(1, 1.2, 1)
basal <- exp(rnorm(n, 1, 1.5))
#basal <- c(1.277015, 5.748106, 2.156822)
U <- 1 + exp(rnorm(n, 1, 1))
Ry <- R / U
#---
h <- as.vector(diag(lam) %*% diag(1 / R) %*% Ci)
A <- diag(h) %*% diag(1 / Ci) %*% Cmat
I <- diag(rep(1, n))
lBetaA <- diag(A %*% t(log(A)))
v <- (I - diag(h)) %*% log(R) - log(basal) + diag(h) %*% lam - lBetaA
lw <- solve(A - I) %*% v
w <- as.vector(exp(lw))
#w
# lw
# h
# basal
hist(lw)
sum(h > 1)
# hist(h)
# hist(log(basal))
# hist(lBetaA)
y <- Ry / w
Q <- y * U
hist(y)
X <- Cmat %*% diag(1 / w)
Xused <- colSums(X)
Q - Xused
