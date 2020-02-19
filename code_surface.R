# x = seq(0.01,30,1)
# y = seq(0.01,30,1)

# z = exp(-x - y) + 0.01

# require(mapplots)


# byx = 1
# byy = 1

# xlim <- c(min(x),max(x))
# ylim <- c(min(y),max(y))

# grd <- make.grid(x , y ,z , byx, byy, xlim, ylim)


# require(fields)

  
# grid.l<-list( abcissa= seq( 0.01,30,1), ordinate= seq( 0.01,30,1)) 
# xg<-make.surface.grid( grid.l)
# # xg is a 300X2 matrix that has all pairs of X and Y grid values 
# z<- exp(- xg[,1] - xg[,2]) + 0.01  
# # now fold z in the matrix format needed for persp 
# out.p<-as.surface( xg, z) 
# persp( out.p) 


f <- function(x,y) exp(-x - y) + 0.01
n <- seq(0.01,30,1)
z <- outer(n,n,f)

## MESH PLOT
persp(x,y,z,theta=30, phi=30, zlim = c(-0.5,1.2),expand=0.6,ticktype='detailed')

## SURFACE PLOT
persp(x,y,z,theta=30, phi=30,zlim = c(-0.5,1.2), expand=0.6,col='lightblue', shade=0.75, ltheta=120,ticktype='detailed')


# require(plot_ly)
# plot_ly (z=z,type = 'surface' )




klen = length(n)
d_x = matrix(0,30,30)

for (iter in seq(1,klen,1) ) {

h = vector()
delta  = vector()

for (i in seq(1,length(n)-1)) {
	h[i] = n[i+1] - n[i]
	delta[i] = ( z[iter,i+1] - z[iter,i] ) /h[i]
}



# for ( i in seq(1,length(x) - 2) ){
# ### Bessel Derivatives
# 	d[[i+1]] = ( h[[i+1]]*delta[[i]] + h[[i]]*delta[[i+1]] ) / (h[[i+1]] + h[[i]])
# }
# d[[1]] = 2*delta[[1]] - d[[2]]
# d[[klen]] = 2*delta[[klen-1]] - d[[klen-1]]

d_1 = delta[1] + (delta[1] - delta[2]) * (h[1]/(h[1] + h[2]))
if (delta[1] == 0 || sign(delta[1]) != sign(d_1)){
	# print("Hello")
	d_x[iter,1] = 0
} else {
	d_x[iter,1] = d_1
}


for ( i in seq(1,length(n) - 2) ){
### Method 2 Preserving the monotonic shape of data (Section 6)
	if(delta[i] == 0 || delta[i+1] == 0){
		d_x[iter,i+1] = 0
	} else{
		d_x[iter,i+1] = ( h[i+1]*delta[i] + h[i]*delta[i+1] ) / (h[i+1] + h[i])
}
}

d_n = delta[klen-1] + (delta[klen-1] - delta[klen - 2]) * (h[klen-1]/(h[klen-2] + h[klen-1]))
if (delta[klen-1] == 0 || sign(delta[klen-1]) != sign(d_n)){
	d_x[iter,klen] = 0
} else {
	d_x[iter,klen] = d_n
}


}

klen = length(n)


# h_mat = matrix(0,length(n)-1,length(n)-1)
# delta_mat = matrix(0,length(n)-1,length(n)-1)


# for (i in seq(1,klen-1,1) ) {
# 	for (j in seq(1,klen-1,1)){
# 		h_mat[i,j] = n[j+1] - n[j]
# 		delta_mat[i,j] = ( z[i,j+1] - z[i,j] ) / h_mat[i,j]

# }
# }


d_y = matrix(0,30,30)

for (iter in seq(1,klen,1) ) {

h = vector()
delta  = vector()

for (i in seq(1,length(n)-1)) {
	h[i] = n[i+1] - n[i]
	delta[i] = ( z[i+1,iter] - z[i,iter] ) /h[i]
}



# for ( i in seq(1,length(x) - 2) ){
# ### Bessel Derivatives
# 	d[[i+1]] = ( h[[i+1]]*delta[[i]] + h[[i]]*delta[[i+1]] ) / (h[[i+1]] + h[[i]])
# }
# d[[1]] = 2*delta[[1]] - d[[2]]
# d[[klen]] = 2*delta[[klen-1]] - d[[klen-1]]

d_1 = delta[1] + (delta[1] - delta[2]) * (h[1]/(h[1] + h[2]))
if (delta[1] == 0 || sign(delta[1]) != sign(d_1)){
	# print("Hello")
	d_y[1,iter] = 0
} else {
	d_y[1,iter] = d_1
}


for ( i in seq(1,length(n) - 2) ){
### Method 2 Preserving the monotonic shape of data (Section 6)
	if(delta[i] == 0 || delta[i+1] == 0){
		d_y[i+1,iter] = 0
	} else{
  }
}

d_n = delta[klen-1] + (delta[klen-1] - delta[klen - 2]) * (h[klen-1]/(h[klen-2] + h[klen-1]))
if (delta[klen-1] == 0 ||  sign(delta[klen-1]) != sign(d_n)){
	d_y[klen,iter] = 0
} else {
	d_y[klen,iter] = d_n
} 

}

 

f_xy <- function(x,y) exp(-x - y)
n_xy <- seq(0.01,30,1)
d_xy <- outer(n_xy,n_xy,f_xy)


# 

# f_x <- function(x,y) -exp(-x - y)
# n_x <- seq(0.01,30,1)
# d_xx <- outer(n_x,n_x,f_x)

tuple <- function(x,y) {
  as.matrix(expand.grid(x=x, y=y))
}



for (i in seq(1,length(n)-1,1)){
	for (j in seq(1,length(n)-1,1)){

		h[i] = 1
		hh[j] = 1
		f_mat = matrix(0,4,4)
		a_i = matrix(0,1,4)
		a_j = matrix(0,1,4)

		q_i = r_ij * (1 - theta)^3 + u_ij * theta * (1 - theta)^2 + v_ij * theta^2 * (1 - theta) + w_ij * theta^3
		a0 = ( r_ij * (1 - theta)^3 + u_ij * theta * (1 - theta)^2 ) / q_i
		a1 = ( v_ij * theta^2 * (1 - theta) + w_ij * theta^3 ) / q_i
		a2 = ( h[i] * r_ij * theta * (1 - theta)^2 ) / q_i
		a3 = ( -h[i] * w_ij * theta^2 * (1 - theta) ) / q_i

		qq_i = rr_ij * (1 - theta_y)^3 + uu_ij * theta_y * (1 - theta_y)^2 + vv_ij * theta_y^2 * (1 - theta_y) + ww_ij * theta_y^3
		aa0 = ( rr_ij * (1 - theta_y)^3 + uu_ij * theta_y * (1 - theta_y)^2 ) / qq_i
		aa1 = ( vv_ij * theta_y^2 * (1 - theta_y) + ww_ij * theta_y^3 ) / qq_i
		aa2 = ( hh[j] * rr_ij * theta_y * (1 - theta_y)^2 ) / qq_i
		aa3 = ( -hh[j] * ww_ij * theta_y^2 * (1 - theta_y) ) / qq_i

		tl = tuple(i,j) ## top left
		tr = tuple(i,j+1) ## top right
		bl = tuple(i+1,j) ## bottom left
		br = tuple(i+1,j+1) ## bottom right

		f_mat[1,1] = z[tl]
		f_mat[1,2] = z[tr]
		f_mat[2,1] = z[bl]
		f_mat[2,2] = z[br]

		f_mat[1,3] = d_y[tl]
		f_mat[1,4] = d_y[tr]
		f_mat[2,3] = d_y[bl]
		f_mat[2,4] = d_y[br]

		f_mat[3,1] = d_x[tl]
		f_mat[3,2] = d_x[tr]
		f_mat[4,1] = d_x[bl]
		f_mat[4,2] = d_x[br]

		f_mat[3,3] = d_xy[tl]
		f_mat[3,4] = d_xy[tr]
		f_mat[4,3] = d_xy[bl]
		f_mat[4,4] = d_xy[br]		

	}
}