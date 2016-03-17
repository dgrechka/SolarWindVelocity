#This test script does approximation of sample data 
# It loads sample data
# approximates it
# draws observations and approximation

obs <- read.csv('SampleData/observations.csv', header = T, sep = ';')
pred <- read.csv('SampleData/prediction.csv', header = T, sep = ';')
colnames(obs) <- c('t','v')
colnames(pred) <- c('t','v')

mv <- 300

filtered_pred <- pred[pred$v != mv,]

joined <- rbind(obs,filtered_pred)

nodes <- joined$t
nodes <- nodes[c(c(T),logical(9))] #take every 10th reference point as node point
if(nodes[length(nodes)]!=joined$t[length(joined$t)])
  nodes <- c(nodes,joined$t[length(joined$t)]) #the last ref point must always present

source("2OrderPolynomialApprox.R")


plot(obs$t,obs$v,
     xlab = "Time (hours since 03.02.2016  18:00:00)",
     ylab = "Velocity (km/s)",
     main = "Solar Wind near Earth\n2nd order piecewise polynomial approximation of observation data",
     xlim=c(0,760),
     ylim=c(250,850),
     sub="Black points - observations; Blue points - predictions; Red curve - trend")

res1_poly <- getApproxF(nodes, joined$t, joined$v)

points(filtered_pred$t,filtered_pred$v,col="blue")
curve(res1_poly,add = T, from = 0, to = 760,col = "red",n=2000,lwd=2)