obs <- read.csv('SampleData/observations.csv', header = T, sep = ';')
colnames(obs) <- c('t','v')

step <- 25
left <-0
right <- 650

nodes <- seq(from=left, to=right, by=step)

poly_f <-function(p,step,t) {
  i <- t %/% step
  return (p[i*3+1]*t*t+p[i*3+2]*t+p[i*3+3]);
}

getP <- function(initSlope,X,Y) {
  res <- numeric((length(X)-1)*3)
  Y1 <- numeric(length(Y))
  Y1[1] <- initSlope
  for(i in 1:(length(X)-1)) {
    x_i <- X[i]
    x_i1 <- X[i+1]
    a <- matrix(
      c(x_i1*x_i1,x_i1,1,
      x_i*x_i,x_i,1,
      2*x_i,1,0), nrow = 3, byrow = T)
    r_part <- c(Y[i+1],Y[i],Y1[i])
    s <- solve(a,r_part)
    Y1[i+1]=2*s[1]*x_i1+s[2]
    res[(i-1)*3+1]=s[1]
    res[(i-1)*3+2]=s[2]
    res[(i-1)*3+3]=s[3]
  }
  return(res)
}

poly_fitted <- function(initSlope,X,Y,step) {
  p1 <- getP(initSlope,X,Y)
  
  return(function(x) {
    poly_f(p1,step,x)
    }
  );
}

testInitSlope <- 0

toMinimize <- function(x) {
  sqErr <- 0
  cur_x_fitted <- poly_fitted(testInitSlope,nodes,x,step)
  for(j in 1:nrow(obs)){
    dif1 <- cur_x_fitted(obs$t[j])-obs$v[j]
    if(! is.na(dif1))
      sqErr <- sqErr + dif1*dif1
  }
  return(sqErr)
}

init <- rep(400, length(nodes))

opt_res1 <- optim(init,toMinimize, , method = "BFGS",control = list(maxit = 500))

plot(obs$t,obs$v,
     xlab = "Time (hours since 03.02.2016  18:00:00)",
     ylab="Velocity (km/s)",
     main = "Solar Wind near Earth\n2nd order piecewise polynomial approximation of observation data")

res1_poly <- poly_fitted(testInitSlope,nodes,opt_res1$par,step)
curve(res1_poly,add = T, from = 0, to = 650)

#testX = c(0,5,10,15,20)
#testY = c(10,5,5,12,14)

#p1 <-getP(-3,testX,testY)



#poly_f1 <- function(x) {
#  return(poly_f(p1,5,x))
#}
#plot(testX,testY,ylim=c(0,20))
#curve(poly_f1,add = T, from = 0, to = 20)
