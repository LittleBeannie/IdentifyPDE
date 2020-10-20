# ------------------------------#
#        data generation        #
# ------------------------------#
N = 50
M = 50
t = seq(0.1, 0.5, length.out = N)
x1 = runif(M,0.1, 0.8) #seq(0.1, 0.5, length.out = M)
x2 = runif(M,0.1, 0.9) #seq(0.2, 0.9, length.out = M)
x = cbind(x1, x2)
v = c(-1,-1) 
u.true = matrix(0, nrow = M, ncol = N) #each row is a time series
for (i in 1:M) {
  for (n in 1:N) {
    u.true[i,n] =  sum( (x[i,] + t[n]*v)^2 )
  }
}
u.obs = u.true + matrix( rnorm( M*N, mean = 0, sd = 0.01), nrow = M, ncol = N)
# ------------------------------#
# derive the derivative w.r.t t #
# ------------------------------#
kernel.fun1 = function(x) {
  output = NULL
  for (i in 1:length(x)) {
    output = c(output, max(0, 3/4*(1 - x[i]^2)))
  }
  return(output)
}

derivative.t = matrix(0, nrow = M, ncol = N)
for (i in 1:M) {
  x.fix.index = i
  for (n in 1:N) {
    t.taylor = t[n]
    X.tep.col = t - t.taylor
    X.tep = cbind(X.tep.col^0, X.tep.col^1, X.tep.col^2, X.tep.col^3, X.tep.col^4)
    y.tep = u.obs[x.fix.index,]
    h.N = N^(-1/7)
    weight.matrix = diag( kernel.fun1(X.tep.col/h.N) )
    b.t = solve(t(X.tep) %*% weight.matrix %*% X.tep) %*% t(X.tep) %*% weight.matrix %*% y.tep
    derivative.t[i,n] = b.t[2]
  }
}
# ------------------------------#
# derive the derivative w.r.t x #
# ------------------------------#
kernel.spt = function(x, x.taylor, bandwidth){
  n = dim(x)[1]
  d = dim(x)[2]
  output.K = matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      element1 = x[i,] - x.taylor
      element2 = x[j,] - x.taylor
      output.K[i,j] = max(0, 3/4*(1-element1^2))/bandwidth * max(0, 3/4*(1-element2^2))/bandwidth
    }
  }
  return(output.K)
}

local.poly.base = function(x, x.taylor){
  n = dim(x)[1]
  output.L = matrix(0, n, n)
  dd = x - x.taylor
  output.L = cbind( rep(1, n), dd, dd[,1]^2, dd[,1]*dd[,2], dd[,2]^2)
  return(output.L)
}

derivative.x.order1 = list(d.x1 = matrix(0, nrow = M, ncol = N), 
                           d.x2 = matrix(0, nrow = M, ncol = N))
derivative.x.order2 = list(d.x11 = matrix(0, nrow = M, ncol = N), 
                           d.x12 = matrix(0, nrow = M, ncol = N),
                           d.x22 = matrix(0, nrow = M, ncol = N))
lambda.spt = 1e-2
for (n in 1:N) {
  t.fix.index = n
  for (i in 1:M) {
    x.taylor.index = i
    y.spt = u.obs[ , t.fix.index]
    L.x = local.poly.base(x, x.taylor = x[x.taylor.index,])
    K.x = kernel.spt(x, x.taylor = x[x.taylor.index,], 1e-3)
    c.x = solve( t(L.x) %*% L.x) %*% t(L.x) %*% K.x %*% solve(K.x + lambda.spt*diag(1,M,M)) %*% y.spt
    derivative.x.order1$d.x1[i,n] = c.x[2]
    derivative.x.order1$d.x2[i,n] = c.x[3]
    derivative.x.order2$d.x11[i,n] = c.x[4]
    derivative.x.order2$d.x12[i,n] = c.x[5]*2
    derivative.x.order2$d.x22[i,n] = c.x[6]
  }
}

# ------------------------------#
#     solve the LASSO problem   #
# ------------------------------#
y.lasso = as.vector(t(derivative.t))
X.lasso.1 = rep(1, M*N)
X.lasso.u = as.vector(t(u.obs))
X.lasso.d1.x1 = as.vector(t(derivative.x.order1$d.x1))
X.lasso.d1.x2 = as.vector(t(derivative.x.order1$d.x2))
X.lasso.d2.x11 = as.vector(t(derivative.x.order2$d.x11))
X.lasso.d2.x12 = as.vector(t(derivative.x.order2$d.x12))
X.lasso.d2.x22 = as.vector(t(derivative.x.order2$d.x22))
X.lasso = cbind(X.lasso.1, X.lasso.u, 
                X.lasso.d1.x1, X.lasso.d1.x2,
                X.lasso.d2.x11, X.lasso.d2.x12, X.lasso.d2.x22)
library(glmnet)
fit = glmnet(X.lasso, y.lasso, nlambda = 10, intercept = F, standardize = F)
plot(fit, xvar="lambda", label=TRUE)
coeff = fit$beta[, length(fit$lambda)]
error[sim] = max(abs(v - coeff[c(3,4)]))
print(c(sim, error[sim]))

