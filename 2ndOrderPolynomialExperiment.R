obs <- read.csv('SampleData/observations.csv', header = T, sep = ';')
colnames(obs) <- c('t','v')
#plot(obs$t,obs$v)

step <- 25

nodes <- seq(from=0, to=650, by=step)

p <- numeric((length(nodes)-1)*3)

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

testX = c(0,5,10,15,20)
testY = c(10,5,5,12,14)

p1 <-getP(-3,testX,testY)



poly_f1 <- function(x) {
  return(poly_f(p1,5,x))
}
plot(testX,testY,ylim=c(0,20))
curve(poly_f1,add = T, from = 0, to = 20)
