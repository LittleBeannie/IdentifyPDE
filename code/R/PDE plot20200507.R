######################################################################
##---------- verify the computational complexity of  ---------------##
##-----------        cubic spline is O(MN)           ---------------##
######################################################################
## N is fixed as 200
y = c(1976723, 16490247, 288505197,1040932354,2614647113,9773435629,25359871947,100195682179,271882080413,597979065447,1148846635279,2007644793713,3270333538947,5045672868379,7455222787013,102327504232947,489376844340079,1501363243120013, 7371187216671079, 22933879424902948)
x = c(10,20,50,75,100,150,200,300,400,500,600,700,800,900,1000,2000,3000,4000,6000,8000)
# M is fixed as 200
x =c(10,20,50,100,200,300,400,500,600,700,1000,2000,4000,6000,8000)
y = c(9323837557,  10167839367, 12699844797, 16919853847,25359871947, 33799890047, 42239908147, 50679926247, 59119944347,67559962447, 92880016747, 177280197747,346080559747,514880921747, 683681283747     )
library(latex2exp)
library(ggplot2)
library(grDevices)
windowsFonts(myFont=windowsFont("Arial")) 
mydata = data.frame(N = x, complexity = y)
g = ggplot(mydata, aes(log(N), log(complexity)))
g = g + geom_line()
g = g + geom_point()
g = g + theme_bw()
g = g + labs(x = "log(N)", y = "log of number of operations") 
g = g + theme(axis.title.x = element_text(size = 20, family="myFont", color="black", face= "plain", vjust=0.5, hjust=0.5))
g = g + theme(axis.title.y = element_text(size = 15, family="myFont", color="black", face= "plain", vjust=0.5, hjust=0.5))
g = g + theme(axis.text.x = element_text(size = 15, family="myFont", color="black", face= "plain", vjust=0.5, hjust=0.5))
g = g + theme(axis.text.y = element_text(size = 15, family="myFont", color="black", face= "plain", vjust=0.5, hjust=0.5))
g = g + theme(plot.title = element_text(hjust = 0.5,size = 25, face = "bold") )
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