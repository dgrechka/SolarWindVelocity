#This test script visually shows how 2ndOrder poly interpolation works

source('2OrderPolynomialApprox.R')


testX = c(0,7,10,18,20)
testY = c(10,5,5,12,14)

p1 <-getP(-3,testX,testY)
poly_p1 <- p_to_polyF(p1,testX)


plot(testX,testY,ylim=c(0,20))
curve(poly_p1,add = T, from = 0, to = 20)