# ------------------------------#
#        data generation        #
# ------------------------------#
N = 80
M = ceiling(N^(8/7))
t = seq(0.1, 0.2, length.out = N)
x1 = seq(0.1, 0.3, length.out = M)
x2 = seq(0.2, 0.3, length.out = M)
x = NULL
for (i in 1:M) {
  for (j in 1:M) {
    x = rbind(x, c(x1[i], x2[j]))  # x1[1], x2[1]
                                   # x1[1], x2[2]
  }
}
v = c(1,1) 
u.true = array(0, dim = c(M, M, N))
for (i in 1:M) {     # u.true[,,1] = x1[1]x2[1] x1[1]x2[2] x1[1]x2[3] ...
  for (j in 1:M) {   #               x1[2]x2[1] x1[2]x2[2] x1[2]x2[3] ...
    for (n in 1:N) {
      u.true[i,j,n] =  sin(sum( (c(x1[i], x2[j]) + t[n]*v)^2 ))
    }
  }
}
u.obs = u.true + array( rnorm( M*M*N, mean = 0, sd = 0.01), dim = c(M, M, N))
# ------------------------------#
# derive the derivative w.r.t t #
# ------------------------------#
kernel.tep = function(x) {
  output = NULL
  for (i in 1:length(x)) {
    output = c(output, max(0, 3/4*(1 - x[i]^2)))
  }
  return(output)
}

derivative.t = array(0, dim = c(M, M, N))
for (i in 1:M) {
  x1.fix.index = i
  for (j in 1:M) {
    x2.fix.index = j
    y.tep = u.obs[x1.fix.index, x2.fix.index,]
    for (n in 1:N) {
      t.taylor = t[n]
      X.tep.col = t - t.taylor
      X.tep = cbind(X.tep.col^0, X.tep.col^1, X.tep.col^2, X.tep.col^3, X.tep.col^4)
      h.N = N^(-1/7)
      weight.matrix = diag( kernel.tep(X.tep.col/h.N) )
      b.t = solve(t(X.tep) %*% weight.matrix %*% X.tep) %*% t(X.tep) %*% weight.matrix %*% y.tep
      derivative.t[i,j,n] = b.t[2]
    }
  }
  
}
# --------------------------------------------#
# derive the derivative w.r.t x (first order) #
# --------------------------------------------#
kernel.spt = function(x, bandwidth){
  n = dim(x)[1]
  d = dim(x)[2]
  output = NULL
  for (i in 1:n) {
      element1 = x[i,1] 
      element2 = x[i,2]
      output = c(output, max(0, 3/4*(1-element1^2))/bandwidth * max(0, 3/4*(1-element2^2))/bandwidth) 
  }
  return(output)
}
derivative.x1 = array(0, dim = c(M, M, N))
derivative.x2 = array(0, dim = c(M, M, N))
for (n in 1:N) {
  t.fix.index = n
  y.spt = as.vector( t(u.obs[ , ,t.fix.index]) ) # first fix x1, then range x2
  for (i in 1:M) {
    for (j in 1:M) {
      x.taylor = c(x1[i], x2[j])
      dd = x - x.taylor
      X.spt.col1 = dd
      X.spt.col2 = cbind( dd[,1]^2, dd[,1]*dd[,2], dd[,2]^2)
      X.spt.col3 = cbind( dd[,1]^3, dd[,1]^2*dd[,2], dd[,1]*dd[,2]^2, dd[,2]^3)
      X.spt.col4 = cbind( dd[,1]^4, dd[,1]^3*dd[,2], dd[,1]^2*dd[,2]^2, dd[,1]*dd[,2]^3, dd[,2]^4 )
      X.spt = cbind(rep(1,M), X.spt.col1, X.spt.col2, X.spt.col3, X.spt.col4)
      
      h.M = M^(-1/10)
      weight.spt = diag(kernel.spt(dd, bandwidth = h.M))
      c.x = solve(t(X.spt) %*% weight.spt %*% X.spt) %*% t(X.spt) %*% weight.spt %*% y.spt
      derivative.x1[i,j,n] = c.x[2]
      derivative.x2[i,j,n] = c.x[3]
    }
  }
}

# ---------------------------------------------#
# derive the derivative w.r.t x (second order) #
# ---------------------------------------------#
derivative.x1x1 = array(0, dim = c(M, M, N))
derivative.x1x2 = array(0, dim = c(M, M, N))
derivative.x2x2 = array(0, dim = c(M, M, N))

for (n in 1:N) {
  t.fix.index = n
  y.spt = as.vector( t(u.obs[ , ,t.fix.index]) ) # first fix x1, then range x2
  for (i in 1:M) {
    for (j in 1:M) {
      x.taylor = c(x1[i], x2[j])
      dd = x - x.taylor
      X.spt.col1 = dd
      X.spt.col2 = cbind( dd[,1]^2, dd[,1]*dd[,2], dd[,2]^2)
      X.spt.col3 = cbind( dd[,1]^3, dd[,1]^2*dd[,2], dd[,1]*dd[,2]^2, dd[,2]^3)
      X.spt.col4 = cbind( dd[,1]^4, dd[,1]^3*dd[,2], dd[,1]^2*dd[,2]^2, dd[,1]*dd[,2]^3, dd[,2]^4 )
      X.spt.col5 = cbind( dd[,1]^5, dd[,1]^4*dd[,2], dd[,1]^3*dd[,2]^2, dd[,1]^2*dd[,2]^3, dd[,1]*dd[,2]^4, dd[,2]^5 )
      X.spt = cbind(rep(1,M), X.spt.col1, X.spt.col2, X.spt.col3, X.spt.col4, X.spt.col5)
      
      h.M = M^(-1/10)
      weight.spt = diag(kernel.spt(dd, bandwidth = h.M))
      c.x = solve(t(X.spt) %*% weight.spt %*% X.spt) %*% t(X.spt) %*% weight.spt %*% y.spt
      derivative.x1x1[i,j,n] = c.x[4]
      derivative.x1x2[i,j,n] = c.x[5]
      derivative.x2x2[i,j,n] = c.x[6]
    }
  }
}



# ------------------------------#
#     solve the LASSO problem   #
# ------------------------------#
# x1[1] x2[1] t[1]
# x1[2] x2[1] t[1]
# x1[3] x2[1] t[1]
# x1[4] x2[1] t[1]
# ...
# x1[M] x2[1] t[1]
# x1[1] x2[2] t[1]
# ...
# x1[M] x2[M] t[1]
# x1[1] x2[1] t[2]
# ...
y.lasso = as.vector(derivative.t)
X.lasso.1 = rep(1, M*M*N)
X.lasso.u = as.vector(u.obs)
X.lasso.d1.x1 = as.vector(derivative.x1)
X.lasso.d1.x2 = as.vector(derivative.x2)
X.lasso.d2.x11 = as.vector(derivative.x1x1)
X.lasso.d2.x12 = as.vector(derivative.x1x2)
X.lasso.d2.x22 = as.vector(derivative.x2x2)
X.lasso = cbind(X.lasso.1, X.lasso.u, 
                X.lasso.d1.x1, X.lasso.d1.x2,
                X.lasso.d2.x11, X.lasso.d2.x12, X.lasso.d2.x22,
                X.lasso.u*X.lasso.d1.x1, X.lasso.u*X.lasso.d1.x2, X.lasso.u*X.lasso.d2.x11, X.lasso.u*X.lasso.d2.x12, X.lasso.u*X.lasso.d2.x22,
                X.lasso.d1.x1*X.lasso.d1.x2, X.lasso.d1.x1*X.lasso.d2.x11, X.lasso.d1.x1*X.lasso.d2.x12, X.lasso.d1.x1*X.lasso.d2.x22,
                X.lasso.d1.x2*X.lasso.d2.x11, X.lasso.d1.x2*X.lasso.d2.x12, X.lasso.d1.x2*X.lasso.d2.x22,
                X.lasso.d2.x11*X.lasso.d2.x12, X.lasso.d2.x11*X.lasso.d2.x22,
                X.lasso.d2.x12*X.lasso.d2.x22)
library(glmnet)
fit = glmnet(X.lasso, y.lasso, nlambda = 10, intercept = F, standardize = F)
plot(fit, xvar="lambda", label=TRUE)
cv_output <- cv.glmnet(X.lasso, y.lasso, alpha = 1, lambda = 10^seq(1, -2, by = -.1))
coef(fit, s = cv_output$lambda.min)
