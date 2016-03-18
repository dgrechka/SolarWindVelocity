ch_data <- read.csv('SampleData/CHBase_0193_R_CH_CSA_2015_1h.csv')

#manually filtering out outliers
ch_data2 <- ch_data[-c(65,83,264,270,750,1423,1576,2495,2374,2659, 2734,2857,2866,2885,2895,2901,2922,2952,2968,2976,3725,6873,7185,7305,7347,7808,8247,8531,8623,8650,10503,11498,11594),]

source("2OrderPolynomialApprox.R")

ch_data3 <- head(ch_data2,3315) # from 30 up to 2140 hours = 2110 hours; about 3 months

X <- ch_data3[,1]
Y <- ch_data3[,2]

nodes <- seq(from=29,to=2141, by=6) #node is every 6 hours

plot(ch_data3)

res1_poly <- getApproxF(nodes, X, Y)

curve(res1_poly,add = T, from = 30, to = 2142,col = "red",n=2000,lwd=2)

x_synth <- seq(from=29,to=2140)
y_synth <- res1_poly(x_synth)

out <- data.frame(t=x_synth,CHarea=y_synth)
write.csv(file="1.ChAreaApproximated.csv",out,row.names = F)