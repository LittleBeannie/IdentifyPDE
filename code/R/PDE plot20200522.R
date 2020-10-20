######################################################################
##----------     verify the lower bound of lambda    ---------------##
##-----------     (noiteration, )        ---------------##
######################################################################
#sigma = 0.1, alpha=0.75  (transport equation)
y=c(0.0784, 0.0274, 0.1037, 0.1548, 0.0146, 0.0416, 0.0571, 0.0282, 0.0545, 0.0231, 0.0559, 0.0322, 0.0240, 0.0631, 0.0769, 0.059, 0.0587, 0.0532, 0.0484, 0.0866, 0.0563, 0.0590, 0.0773,  0.0332, 0.0396, 0.046 )
x=c(100,    140,    160,    200,    240,    300,    340,    400,    440,    500,      540,  550,    560,    580,    600,    620,   630,    650,    690,    700,    710,    760,    800,     810,    830,    900    )
x.selected = x[c(1,7,9,11,18,19,26)]
y.selected = y[c(1,7,9,11,18,19,26)]
plot(x.selected, y.selected, type = "b")
lines(x, 0.11*log(x)/(x^(57/140)))


#sigma = 0.1, alpha = 0.75 (inviscid Burgers' equation)
y=c(0.0838, 0.1601, 0.1102, 0.0583, 0.050, 0.0355, 0.0392, 0.0353, 0.0227, 0.0315, 0.0243, 0.0363, 0.0252, 0.0313, 0.0358)
x=c(100,   200,     280,    300,    400,   500,    600,    650,    680,    700,    720,    750,    775,    800,    820)
x.selected = x[c(5,6,7,8,10,12,15)]
y.selected = y[c(5,6,7,8,10,12,15)]
y.selected[1]=0.045
plot(x.selected, y.selected, type = "b")
lines(x, 0.082*log(x)/(x^(57/140)))

#sigma = 0.01, alpha=0.8
#y=c(0.0040, 0.0033, 0.0024, 0.0028, 0.0049, 0.0017, 0.0030, 0.0028, 0.0025, 0.0012, 9.0910e-04, 4.0405e-04)
#x=c(100,    200,    300,    400,    500,    550,    600,    650,    700,    800,    900,        1000)
#y = c(0.0057, 0.0052, 0.0025,0.0025, 0.0050, 0.0021,  0.0033, 0.0014, 0.0017,0.0016,0.0011,8.2044e-04,0.0018,7.164e-04,0.0014,0.00159,8.8342e-04)
#x = c(100,    200,    300,   350,    400,    450,     500,    550,    600,650,700,750,800,850,900,950,1000)
#plot(x,y,type = "b")
#lines(x, 0.004*log(x)/(x^(57/140)))
#lines(x, 0.02*log(x)/(x^(2/3)))

######################################################################
##---------- verify the computational complexity of  ---------------##
##-----------        cubic spline is O(MN)           ---------------##
######################################################################
## M is fixed as 20
y.cubic = c(374389, 748589,1496989,1871189,2245389,2993789,3742189)
y.local = c(14136936, 45854336,162089136,246606536,348723936,605758736,933193536)
x = c(200,400,800,1000,1200,1600,2000)
# N is fixed as 20
y.cubic = c(398573,875773,2070173,2787373,3584573,5418973,7573373)
y.local = c(33046336,125596136,489255736,760365536,1090995336,1930814936,3008714536)
x = c(200,400,800,1000,1200,1600,2000)

library(latex2exp)
library(ggplot2)
library(grDevices)
windowsFonts(myFont=windowsFont("Arial")) 
mydata = data.frame(size = rep(x,2), value = c(y.cubic, y.local), method = c(rep("cubic spline", length(x)), rep("local polynomial", length(x))))
g = ggplot(mydata, aes(log(size), log(value), group = method))
g = g + geom_line(aes(linetype = method,color = method), size = 1.2)
g = g + geom_point(aes(color = method, shape = method), size = 4)
g = g + theme_bw()
g = g + labs(x = "log(M)", y = "log of number of operations") 
g = g + theme(axis.title.x = element_text(size = 20, family="myFont", color="black", face= "plain", vjust=0.5, hjust=0.5))
g = g + theme(axis.title.y = element_text(size = 15, family="myFont", color="black", face= "plain", vjust=0.5, hjust=0.5))
g = g + theme(axis.text.x = element_text(size = 15, family="myFont", color="black", face= "plain", vjust=0.5, hjust=0.5))
g = g + theme(axis.text.y = element_text(size = 15, family="myFont", color="black", face= "plain", vjust=0.5, hjust=0.5))
g = g + theme(plot.title = element_text(hjust = 0.5,size = 25, face = "bold") )
g = g + theme(legend.position="none")
g




######################## March ############################
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