
poly_f <- function(p,x,t) {
  i <- findInterval(t,x)-1
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

p_to_polyF <- function(p,nodes) {
  return(function(x) {
    poly_f(p,nodes,x)
  }
  );
}

getApproxF <- function(nodes,x,y) {
  testInitSlope <- 0
  
  init <- rep(mean(y), length(nodes))
  
  toMinimize <- function(x1) {
  sqErr <- 0
  p1 <- getP(testInitSlope,nodes,x1)
  cur_fitted <- p_to_polyF(p1,nodes)
  for(j in 1:length(x)){
    dif1 <- cur_fitted(x[j])-y[j]
    if(! is.na(dif1)) {
      sqErr <- sqErr + dif1*dif1 #difference with measured data
    }
  }
  return(sqErr)
  }
  
  opt_res1 <- optim(init,toMinimize, method = "BFGS",control = list(maxit = 500))
  
  res1_p <- getP(testInitSlope,nodes,opt_res1$par)
  res1_poly <- p_to_polyF(res1_p,nodes)
  return(res1_poly)
}

