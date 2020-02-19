setwd("Documents/MTP/Part_1/After_MidSem/Code/gregory")




x = c(-2,-1,-0.3,-0.2)
f = 1/(x^2)
d = -2/(x^3)

x = seq(-5,5,1)
f = c(10,5,-3,-15,-45,-30,15,25,40,66,70)
d = c(-1,-5,-6,-5,-10,30,15,5,4,3,2)

x = c(-10,-5,2,6,16,30)
f = c(10,-3,-15,15,40,70)
d = c(-1,-6,-10,1,4,8)


# plot(x,f,lwd = 2)


a = min(x)
b = max(x)

# png(filename="Cubic_Spline_Interpolation.png")
# plot(spline(x, f, n = 10000*length(x),method ="natural"),lwd=2, cex = 0.3,xlab = "x", ylab = "y", main = "Cubic Spline(Natural) ")
# points(x,f,lwd=6,col = 'red')
# dev.off()

####### HERMITE INTERPOLATION #####

hermite = splinefunH(x, f, d)
ls(envir = environment(hermite))
xs = seq(min(x),max(x),1/10000)

png(filename="Hermite_Spline_Interpolation.png",height = 800,width = 950)
plot(xs,hermite(xs),lwd=2, cex = 0.3,xlab = "x", ylab = "y", main = "Hermite Spline ")
points(x,f,lwd=6,col = 'red')
dev.off()




###### CONVEXITY PRESERVING RATIONAL INTERPOLATION #######

x_s = list()
y_s = list()

for ( i in seq(1,length(x) - 1) ){
	#print (i)
	x_1 = x[i]
	x_2 = x[i+1]
	f_1 = f[i]
	f_2 = f[i+1]
	d_1 = d[i]
	d_2 = d[i+1]
	h_i = x_2 - x_1

	xs = seq(x_1,x_2,0.0001)

	delta_i = (f_2 - f_1) / (x_2 - x_1)
	z = (xs - x_1) / h_i


	r_i = 1 + (d_2 - delta_i)/(delta_i - d_1) + (delta_i - d_1)/(d_2 - delta_i)
	p <-  f_2 * (z^3)  + ( (r_i * f_2) - (h_i * d_2) )*( (z^2) * (1 - z) ) +  ( (r_i * f_1) + (h_i * d_1) )*( (z) * ( (1 - z)^2 ) ) + f_1*( (1 - z)^3 )
	q <- 1 + (r_i - 3)*(z * (1 - z))

	######### learn append to a vector #########
	x_s[[i]] = xs
	#print (x_s[[i]][1:50])
	y_s[[i]] = p/q


} 

####yy[duplicated(yy) == TRUE]

xx = unlist(x_s)
yy = unlist(y_s)

png(filename="Convexity_Preserving_Rational_Interpolation.png",height = 800,width = 950)
plot(xx,yy,lwd = 0.5,cex = 0.5,xlab = "x", ylab = "y",ylim = c(min(f)-2,max(f)+2), main = "Convexity_Preserving_Rational_Interpolation")
points(x,f,lwd=6,col = "red")
dev.off()






##### MONOTONICITY PRESERVING RATIONAL INTERPOLATION #######



x_s = list()
y_s = list()

for ( i in seq(1,length(x) - 1) ){
	#print (i)
	x_1 = x[i]
	x_2 = x[i+1]
	f_1 = f[i]
	f_2 = f[i+1]
	d_1 = d[i]
	d_2 = d[i+1]
	h_i = x_2 - x_1

	xs = seq(x_1,x_2,0.0001)

	delta_i = (f_2 - f_1) / (x_2 - x_1)
	z = (xs - x_1) / h_i


	r_i = 1 + ( (d_1 + d_2) / (delta_i) )
	p <-  f_2 * (z^3)  + ( (r_i * f_2) - (h_i * d_2) )*( (z^2) * (1 - z) ) +  ( (r_i * f_1) + (h_i * d_1) )*( (z) * ( (1 - z)^2 ) ) + f_1*( (1 - z)^3 )
	q <- 1 + (r_i - 3)*(z * (1 - z))

	######### learn append to a vector #########
	x_s[[i]] = xs
	#print (x_s[[i]][1:50])
	y_s[[i]] = p/q


} 

####yy[duplicated(yy) == TRUE]

xx = unlist(x_s)
yy = unlist(y_s)

png(filename="Monotonicity_Preserving_Rational_Interpolation.png",height = 800,width = 950)
plot(xx,yy,lwd = 0.5,cex = 0.5,xlab = "x", ylab = "y",main = "Monotonicity_Preserving_Rational_Interpolation")
points(x,f,lwd=6,col = "red")
dev.off()






##### APPROXIMATIONS FOR THE DERIVATIVE PARAMETERS #####

h = vector()
delta = vector()
for (i in seq(1, length(x) - 1)){
	x_1 = x[i]
	x_2 = x[i+1]
	f_1 = f[i]
	f_2 = f[i+1]
	h_i = x_2 - x_1
	delta_i = (f_2 - f_1) / (x_2 - x_1)
	h[i] = h_i
	delta[i] = delta_i
}

#### Three point Difference Approximation #### 
d = vector()
for (i in seq(1, length(x))){
	if(i == 1){
		d[i] = (1 + h[i]/h[i+1])*delta[i] - ( (h[i]/h[i+1]) * ( (f[3] - f[1])/(x[3] - x[1]) ) )
	} else if(i == length(x)){
		d[i] = (1 + (h[i-1]/h[i-2]) ) - (h[i-1]/h[i-2]) * ( (f[i] - f[i-2])/(x[i] - x[i-2]) ) 
	} else{
		d[i] = (h[i]*delta[i-1] + h[i-1]*delta[i]) / (h[i-1] + h[i])
	}
}



###### CONVEXITY PRESERVING RATIONAL INTERPOLATION #######

x_s = list()
y_s = list()

for ( i in seq(1,length(x) - 1) ){
	#print (i)
	x_1 = x[i]
	x_2 = x[i+1]
	f_1 = f[i]
	f_2 = f[i+1]
	d_1 = d[i]
	d_2 = d[i+1]
	h_i = x_2 - x_1

	xs = seq(x_1,x_2,0.00001)

	delta_i = (f_2 - f_1) / (x_2 - x_1)
	z = (xs - x_1) / h_i


	r_i = 1 + (d_2 - delta_i)/(delta_i - d_1) + (delta_i - d_1)/(d_2 - delta_i)
	p <-  f_2 * (z^3)  + ( (r_i * f_2) - (h_i * d_2) )*( (z^2) * (1 - z) ) +  ( (r_i * f_1) + (h_i * d_1) )*( (z) * ( (1 - z)^2 ) ) + f_1*( (1 - z)^3 )
	q <- 1 + (r_i - 3)*(z * (1 - z))

	######### learn append to a vector #########
	x_s[[i]] = xs
	#print (x_s[[i]][1:50])
	y_s[[i]] = p/q


} 

####yy[duplicated(yy) == TRUE]

xx = unlist(x_s)
yy = unlist(y_s)

png(filename="Convexity_Preserving_Rational_Interpolation_with_three_point_derivatives.png",height = 800,width = 950)
plot(xx,yy,lwd = 0.5,cex = 0.5,xlab = "x", ylab = "y",ylim = c(min(f) - 3 , max(f) + 3 ) , main = "Convexity_Preserving_Rational_Interpolation_with_three_point_derivatives" )
points(x,f,lwd=6,col = "red")
dev.off()






##### GEOMETRIC MEAN APPROXIMATION #####

d = vector()
for (i in seq(1, length(x))){
	if(i == 1){
		d[i] = delta[i]^(1 + (h[i]/h[i+1]) ) *   ((f[3] - f[1])/(x[3] - x[1]))^(-h[i]/h[i+1]) 
	} else if(i == length(x)){
		d[i] = delta[i-1]^(1 + (h[i-1]/h[i-2]) ) * ( (f[i] - f[i-2])/(x[i] - x[i-2]) )^(-h[i-1]/h[i-2]) 
	} else{
		d[i] = delta[i-1]^h[i] +  delta[i]^(h[i-1] / (h[i-1] + h[i]) )
	}
}



###### CONVEXITY PRESERVING RATIONAL INTERPOLATION #######

x_s = list()
y_s = list()

for ( i in seq(1,length(x) - 1) ){
	#print (i)
	x_1 = x[i]
	x_2 = x[i+1]
	f_1 = f[i]
	f_2 = f[i+1]
	d_1 = d[i]
	d_2 = d[i+1]
	h_i = x_2 - x_1

	xs = seq(x_1,x_2,0.00001)

	delta_i = (f_2 - f_1) / (x_2 - x_1)
	z = (xs - x_1) / h_i


	r_i = 1 + (d_2 - delta_i)/(delta_i - d_1) + (delta_i - d_1)/(d_2 - delta_i)
	p <-  f_2 * (z^3)  + ( (r_i * f_2) - (h_i * d_2) )*( (z^2) * (1 - z) ) +  ( (r_i * f_1) + (h_i * d_1) )*( (z) * ( (1 - z)^2 ) ) + f_1*( (1 - z)^3 )
	q <- 1 + (r_i - 3)*(z * (1 - z))

	######### learn append to a vector #########
	x_s[[i]] = xs
	#print (x_s[[i]][1:50])
	y_s[[i]] = p/q


} 

####yy[duplicated(yy) == TRUE]

xx = unlist(x_s)
yy = unlist(y_s)

png(filename="Convexity_Preserving_Rational_Interpolation_with_geometric_derivatives.png",height = 800,width = 950)
plot(xx,yy,lwd = 0.5,cex = 0.5,xlab = "x", ylab = "y", ,main = "Convexity_Preserving_Rational_Interpolation_with_geometric_derivatives")
points(x,f,lwd=6,col = "red")
dev.off()

