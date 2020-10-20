error = rep(0,1000)
for (sim in 1:1000) {
  
  N = 50
  M = ceiling(N^(9/7))
  t = seq(0.1, 0.5, length.out = N)
  x1 = rnorm(M, mean = 0.4, sd = 0.05)
  x2 = rnorm(M, mean = 0.4, sd = 0.05)
  x = cbind(x1, x2)
  v = c(-1,-1) 
  u.true = matrix(0, nrow = M, ncol = N) #each row is a time series
  for (i in 1:M) {
    for (n in 1:N) {
      u.true[i,n] =  sqrt(sum( (x[i,] + t[n]*v)^2 ) )
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
  
  # ---------------------------------------------#
  # derive the derivative w.r.t x  (first order) #
  # ---------------------------------------------#
  kron.n = function(x,n){
    output = 1
    for (i in 1:n) {
      output = kronecker(output, x)
    }
    return(output)
  }
  kernel.fun2 = function(x, bandwidth) {
    n = dim(x)[1]
    d = dim(x)[2]
    output = NULL
    for (i in 1:n) {
      #p1 = exp(-x[i,1]^2/2/bandwidth^2)
      #p2 = exp(-x[i,2]^2/2/bandwidth^2)
      p1 = max(0, 3/4*(1-x[i,1]^2))/bandwidth
      p2 = max(0, 3/4*(1-x[i,2]^2))/bandwidth
      output = c(output, p1*p2)
    }
    return(output)
  }
  
  derivative.x.order1 = list(d.x1 = matrix(0, nrow = M, ncol = N), 
                             d.x2 = matrix(0, nrow = M, ncol = N))
  
  for (n in 1:N) {
    t.fix.index = n
    for (i in 1:M) {
      x.taylor = x[i,]
      dd = x - x.taylor
      X.spt.col1 = dd
      X.spt.col2 = cbind( dd[,1]^2, 2*dd[,1]*dd[,2], dd[,2]^2)
      X.spt.col3 = cbind( dd[,1]^3, 3*dd[,1]^2*dd[,2], 3*dd[,1]*dd[,2]^2, dd[,2]^3)
      X.spt.col4 = cbind( dd[,1]^4, 4*dd[,1]^3*dd[,2], 4*dd[,1]^2*dd[,2]^2, 4*dd[,1]*dd[,2]^3, dd[,2]^4 )
      
      X.spt = cbind(rep(1,M), X.spt.col1, X.spt.col2, X.spt.col3, X.spt.col4)
      y.spt = u.obs[, t.fix.index]
      h.M = M^(-1/10)
      weight.spt = diag(kernel.fun2(x - x.taylor, bandwidth = h.M))
      c.t = solve(t(X.spt) %*% weight.spt %*% X.spt) %*% t(X.spt) %*% weight.spt %*% y.spt
      derivative.x.order1$d.x1[i,n] = c.t[2]
      derivative.x.order1$d.x2[i,n] = c.t[3]
    }
  }
  
  # ---------------------------------------------#
  # derive the derivative w.r.t x (second order) #
  # ---------------------------------------------#
  derivative.x.order2 = list(d.x11 = matrix(0, nrow = M, ncol = N), 
                             d.x12 = matrix(0, nrow = M, ncol = N),
                             d.x22 = matrix(0, nrow = M, ncol = N))
  for (n in 1:N) {
    t.fix.index = n
    for (i in 1:M) {
      x.taylor = x[i,]
      dd = x - x.taylor
      X.spt.col1 = dd
      X.spt.col2 = cbind( dd[,1]^2, 2*dd[,1]*dd[,2], dd[,2]^2)
      X.spt.col3 = cbind( dd[,1]^3, 3*dd[,1]^2*dd[,2], 3*dd[,1]*dd[,2]^2, dd[,2]^3)
      X.spt.col4 = cbind( dd[,1]^4, 4*dd[,1]^3*dd[,2], 4*dd[,1]^2*dd[,2]^2, 4*dd[,1]*dd[,2]^3, dd[,2]^4 )
      X.spt.col5 = cbind( dd[,1]^5, 5*dd[,1]^4*dd[,2], 5*dd[,1]^3*dd[,2]^2, 5*dd[,1]^2*dd[,2]^3, 5*dd[,1]*dd[,2]^4, dd[,2]^5 )
      
      X.spt = cbind(rep(1,M), X.spt.col1, X.spt.col2, X.spt.col3, X.spt.col4, X.spt.col5)
      y.spt = u.obs[, t.fix.index]
      h.M = M^(-1/10)
      weight.spt = diag(kernel.fun2(x - x.taylor, bandwidth = h.M))
      c.t = solve(t(X.spt) %*% weight.spt %*% X.spt) %*% t(X.spt) %*% weight.spt %*% y.spt
      derivative.x.order2$d.x11[i,n] = c.t[4] 
      derivative.x.order2$d.x12[i,n] = c.t[5] 
      derivative.x.order2$d.x22[i,n] = c.t[6] 
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
}






















