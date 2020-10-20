N = 20
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
v = c(1,0) 
u.true = array(0, dim = c(M, M, N))
for (i in 1:M) {     # u.true[,,1] = x1[1]x2[1] x1[1]x2[2] x1[1]x2[3] ...
  for (j in 1:M) {   #               x1[2]x2[1] x1[2]x2[2] x1[2]x2[3] ...
    for (n in 1:N) {
      u.true[i,j,n] =  sin(sum( (c(x1[i], x2[j]) + t[n]*v)^2 ))
    }
  }
}
u.obs = u.true + array( rnorm( M*M*N, mean = 0, sd = 0.02), dim = c(M, M, N))

obs.plot <- list(dim1 = x1, dim2 = x2, dim3 = u.obs[,,1])
fig.obs <- plot_ly(x = obs.plot$dim1, y = obs.plot$dim2, z = obs.plot$dim3, showscale = FALSE) %>% add_surface()
fig.obs <- fig %>% add_surface(z = u.obs[,,10], opacity = 0.98)
fig.obs <- fig %>% add_surface(z = u.obs[,,20], opacity = 0.98)
fig.obs

true.plot <- list(dim1 = x1, dim2 = x2, dim3 = u.true[,,1])
fig.true <- plot_ly(x = true.plot$dim1, y = true.plot$dim2, z = true.plot$dim3, showscale = FALSE) %>% add_surface()
fig.true <- fig %>% add_surface(z = u.true[,,10], opacity = 0.98)
fig.true <- fig %>% add_surface(z = u.true[,,20], opacity = 0.98)
fig.true