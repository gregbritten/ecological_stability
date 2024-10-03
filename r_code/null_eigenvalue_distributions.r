library(MASS) #package for 2D kernel density estimation

##############################################################################
## PARAMETERS ################################################################
##############################################################################
ps <- c(2,3,5,10,20,50,100,200,500)        #different choices for matrix dimension

circ    <- seq(0,2*pi,length=200)          #points to draw unit circle
circ_xy <- t(rbind( 0+sin(t), 0+cos(t)))   #coordinates of unit circle

###############################################################################
## DISTRIBUTION OF EIGENVALUES FOR A SINGLE RANDOM MATRIX #####################
###############################################################################

##--Gaussian random matrix--############
par(mfrow=c(3,3),mar=c(2,2,2,2),oma=c(3,3,3,3))
for(i in 1:length(ps)){
	p <- ps[i]                                           
	eig <- eigen(matrix(rnorm(p*p, 0, 1),ncol=p))      #draw p*p standard normal realizations
	plot(Re(eig$values)/sqrt(p), Im(eig$values)/sqrt(p),ylim=c(-1.5,1.5),xlim=c(-1.5,1.5)) #scatterplot of real vs. imaginary eigenvalues
	lines(circ_xy,type='l',lwd=2,col='blue')           #plot unit circle
}
mtext(expression('Real Eigenvalues (Re('*lambda['i']*')/'*sqrt(p)*')'),outer=TRUE,side=1,line=1)
mtext(expression('Imaginary Eigenvalues (Im( '*lambda['i']*')/'*sqrt(p)*')'),outer=TRUE,side=2,line=1)

##--Uniform random matrix--############
par(mfrow=c(3,3),mar=c(2,2,2,2),oma=c(3,3,3,3))
for(i in 1:length(ps)){
	p <- ps[i]
	eig <- eigen(matrix(runif(p*p, -2, 2),ncol=p))     #draw p*p uniform realizations
	plot(Re(eig$values)/sqrt(p), Im(eig$values)/sqrt(p),ylim=c(-1.5,1.5),xlim=c(-1.5,1.5))
	lines(circ_xy,type='l',lwd=2,col='blue')
}
mtext(expression('Real Eigenvalues (Re('*lambda['i']*')/'*sqrt(p)*')'),outer=TRUE,side=1,line=1)
mtext(expression('Imaginary Eigenvalues (Im( '*lambda['i']*')/'*sqrt(p)*')'),outer=TRUE,side=2,line=1)


###############################################################################
## DISTRIBUTION OF EIGENVALUES FOR A SINGLE RANDOM MATRIX #####################
###############################################################################
niter <- 50 
ps <- c(2,3,5,10,20,50,100,200,500)        #different choices for matrix dimension

for(j in 1:length(ps)){
	p <- ps[j]
	eigs_re = eigs_im <- matrix(NA,ncol=p,nrow=niter)

		for(i in 1:niter){
			eig <- eigen(matrix(rnorm(p*p, 0, 1),ncol=p))
			eigs_re[i,] <- Re(eig$values)
			eigs_im[i,] <- Im(eig$values)
		}

	plot(c(eigs_re)/sqrt(p),c(eigs_im)/sqrt(p),ylim=c(-1.5,1.5),xlim=c(-1.5,1.5))
	lines(circ_xy,type='l',lwd=2,col='blue')
}


##############################################################################
## KERNEL DENSITY ############################################################
##############################################################################
par(mfrow=c(2,2))
kde <- kde2d(x=c(eigs_re)/sqrt(p), y=c(eigs_im)/sqrt(p))

filled.contour(kde)






