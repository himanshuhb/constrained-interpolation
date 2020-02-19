setwd("D:\\MTP/Code")

a <- c(0,1,10)
b <- c(1,3,6)

x1 <- a[1]
x2 <- b[1]
y1 <- a[2]
y2 <- b[2]
d1 <- a[3]
d2 <- b[3]

h <- x2 - x1
delta <- 2
alpha <- 5
beta <- 3

require(polynom)

p <- function(x){
	( ( d1 + d2 - 2*(delta) ) / (h^2) ) * (x - x1)^3   +    ( ( -2*d1 - d2 + 3*delta ) / h ) * (x - x1)^2    +    d1(x - x1)    +   y1
}


polynomial(rev(c((( d1 + d2 - 2*(delta) ) / (h^2) ),( ( -2*d1 - d2 + 3*delta ) / h ),d1,y1)))  -> p

change.origin(p,-x1) -> p

xs <- seq(min(x1,x2),max(x1,x2),0.00001)

png(filename="Hermite_Cubic_ODE.png")
plot(smooth.spline(xs,predict(p,xs)),lwd=1,cex=0.5, xlab = "x", ylab = "y", main = "Hermite cubic spline")
points(c(x1,y1),c(x2,y2),col = "black",lwd = 8)
dev.off()



## monotonicity region
ellipse <- function(x,y){
	x^2 + y^2 + x*y - 6*x - 6*y     # + 9 omitted
}



## for alpha >= beta

K = 1.1 * (alpha / 3) 
r = 0.8 * ( (4 - beta) / (4*K - beta) )
beta_r = beta*( (1 - r) / (1 - r*K))

for ( i in seq(1,3,0.1) ) {
	
	if (beta_r > 3) {
		if ( ellipse(i,beta_r) < (-9) ) {
			alpha_r = i
			break
		} 
	}
	else {
		alpha_r = i
		break
	}
	
}


xnew = x1 + r*(x2 - x1)
ynew = y1 + r*(x2 - x1)*K*delta
dnew = alpha_r*delta*( (1 - r*K) / (1 - r) )




delta_l <- (ynew - y1) / (xnew - x1)
delta_r <- (y2 - ynew) / (x2 - xnew) 
h_l <- xnew - x1
h_r <- x2 - xnew


polynomial(rev(c((( d1 + dnew - 2*(delta_l) ) / (h_l^2) ),( ( -2*d1 - dnew + 3*delta_l ) / h_l ),d1,y1)))  -> m
change.origin(m,-x1) -> m
polynomial(rev(c((( dnew + d2 - 2*(delta_r) ) / (h_r^2) ),( ( -2*dnew - d2 + 3*delta_r ) / h_r ),dnew,ynew)))  -> n
change.origin(n,-xnew) -> n



xs_l <- seq(min(x1,xnew),max(x1,xnew),0.0001)
xs_r <- seq(min(xnew,x2),max(xnew,x2),0.0001)

all_x <- c(x1,x2,xnew)
all_y <- c(y1,y2,ynew)
all_point <- cbind(all_x,all_y)

png(filename="ODE_Algo.png")
plot(smooth.spline(c(xs_l,xs_r),c(predict(m,xs_l),predict(n,xs_r))),lwd=1,cex=0.5,xlab = "x", ylab = "y", main = "Plot after adding knots")
points(all_point,col = "black",lwd = 8)
dev.off()
#par(new =TRUE)
# points(smooth.spline(xs_r,predict(n,xs_r)),lwd=2,col="green")





cbind(c(xs_l,xs_r),c(predict(m,xs_l),predict(n,xs_r))) -> data 
plot(data)




### testing
predict(m,xnew)
predict(n,xnew)