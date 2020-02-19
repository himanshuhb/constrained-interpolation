setwd("Documents/Academics/MTP/Part_2/Before_MidSem/Code/Range")

x = c(0,1,2,3,4)
f = c(0.4,0.9,0.9,0.4,0.9)
l = c(0.32,0.82,0.82,0.32,0.82)
u_given = c(0.42,0.92,0.92,0.42,0.92)

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

klen = length(x)

for ( i in seq(1,length(x) - 2) ){
### Bessel Derivatives
	d[[i+1]] = ( h[[i+1]]*delta[[i]] + h[[i]]*delta[[i+1]] ) / (h[[i+1]] + h[[i]])
}
d[[1]] = 2*delta[[1]] - d[[2]]
d[[klen]] = 2*delta[[klen-1]] - d[[klen-1]]

der = unlist(d)

##### CUBIC SPLINE #######
cubicspline = splinefunH(x, f, unlist(d))

xx = seq(x[1],x[klen],0.0001)
cubicspline(xx) -> y
plot(xx,y,lwd=0.5,cex=0.5,ylim = c(0.3,1.0),xlab = "x",ylab = "y",main = "Cubic Spline")

for ( i in seq(1,length(x) - 1) ){
	x_1 = x[i]
	x_2 = x[i+1]
	l_1 = l[i]
	l_2 = l[i+1]
	u_1 = u_given[i]
	u_2 = u_given[i+1]
	lines(c(x_1,x_2),c(l_1,l_2),lwd=1.25,col = "red")
	lines(c(x_1,x_2),c(u_1,u_2),lwd=1.25,col = "red")
}

png(filename="cubic_with_range_0.2.png",height = 800,width=800)
plot(xx,y,lwd=0.5,cex=0.5,ylim = c(0.3,1.0),xlab = "x",ylab = "y",main = "Cubic Spline")

for ( i in seq(1,length(x) - 1) ){
	x_1 = x[i]
	x_2 = x[i+1]
	l_1 = l[i]
	l_2 = l[i+1]
	u_1 = u_given[i]
	u_2 = u_given[i+1]
	lines(c(x_1,x_2),c(l_1,l_2),lwd=1.25,col = "red")
	lines(c(x_1,x_2),c(u_1,u_2),lwd=1.25,col = "red")
}
dev.off()




x_s = list()
y_s = list()

for ( i in seq(1,length(x) - 1) ){
	#print (i)
	x_1 = x[i]
	x_2 = x[i+1]
	f_1 = f[i]
	f_2 = f[i+1]
	h_i = x_2 - x_1
	delta_i = (f_2 - f_1) / (x_2 - x_1)
	d_1 = der[i]
	d_2 = der[i+1]
	l_1 = l[i]
	l_2 = l[i+1]
	u_1 = u_given[i]
	u_2 = u_given[i+1]

	xs = seq(x_1,x_2,0.0001)
	#z = (xs - x_1) / h_i

	u = (xs - x_1) / h_i
	v = (x_2 - xs) / h_i

	Aa = ( delta_i - d_1 ) / (f_1 - l_1)
	Bb = ( delta_i - d_1 ) / (f_1 - u_1)
	Cc = ( d_2 - delta_i ) / (f_2 - l_2)
	Dd = ( d_2 - delta_i ) / (f_2 - u_2)

	Zz = max(Aa,Bb,Cc,Dd)

	Kk = h_i * Zz

	alpha = -3.5 + Kk


	#r_i = 1 + (d_2 - delta_i)/(delta_i - d_1) + (delta_i - d_1)/(d_2 - delta_i)
	
	#p <-  f_2 * (z^3)  + ( (r_i * f_2) - (h_i * d_2) )*( (z^2) * (1 - z) ) +  ( (r_i * f_1) + (h_i * d_1) )*( (z) * ( (1 - z)^2 ) ) + f_1*( (1 - z)^3 )
	#q <- 1 + (r_i - 3)*(z * (1 - z))

	p = ((d_1 - delta_i)*v) + ((delta_i - d_2)*u)
	q = 1 + alpha*u*v

	y_s[[i]] = (f_1 * v) + (f_2 * u) + ((h_i*u*v)*(p/q))

	######### learn append to a vector #########
	x_s[[i]] = xs
	#print (x_s[[i]][1:50])


} 

####yy[duplicated(yy) == TRUE]

xx = unlist(x_s)
yy = unlist(y_s)

plot(xx,yy,lwd = 0.5,cex=0.25,ylim = c(0.3,1.0),xlab = "x",ylab = "y",main = "Range restricted intepolation")
for ( i in seq(1,length(x) - 1) ){
	x_1 = x[i]
	x_2 = x[i+1]
	l_1 = l[i]
	l_2 = l[i+1]
	u_1 = u_given[i]
	u_2 = u_given[i+1]
	lines(c(x_1,x_2),c(l_1,l_2),lwd=1.25,col = "red")
	lines(c(x_1,x_2),c(u_1,u_2),lwd=1.25, col = "red")
}


png(filename="range_restricted_0.2.png",height = 800,width=800)

plot(xx,yy,lwd = 0.5,cex=0.5,ylim = c(0.3,1.0),xlab = "x",ylab = "y",main = "Range restricted intepolation")
for ( i in seq(1,length(x) - 1) ){
	x_1 = x[i]
	x_2 = x[i+1]
	l_1 = l[i]
	l_2 = l[i+1]
	u_1 = u_given[i]
	u_2 = u_given[i+1]
	lines(c(x_1,x_2),c(l_1,l_2),lwd=1.25,col = "red")
	lines(c(x_1,x_2),c(u_1,u_2),lwd=1.25, col = "red")
}

dev.off()




