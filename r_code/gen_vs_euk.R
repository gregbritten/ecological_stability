library(rstan)
options(mc.cores=parallel::detectCores())
################################################################
## COMPILE MODEL ###############################################
################################################################
mod <- stan_model('~/dropbox/working/ecological_stability/r_code/MVAR1.stan')

################################################################
## READ DATA ###################################################
################################################################
euk <- read.csv('~/dropbox/data/martin_platero_etal_2018/csvs/eukaryotic_OTUs.csv')
bac <- read.csv('~/dropbox/data/martin_platero_etal_2018/csvs/bacterial_OTUs.csv')
env <- read.csv('~/dropbox/data/martin_platero_etal_2018/csvs/environment.csv')

###########################################
##--AGGREGATE BY GENUS--###################
###########################################
taxa <- c('Phyla','Class','Order','Family','Genus')

TAX_bac <- list()
TAX_euk <- list()
for(i in 3:7){
	nms_bac <- unique(bac[,i])
		nms_bac[nchar(as.character(nms_bac))>3]
	nms_euk <- unique(euk[,i])
		nms_euk[nchar(as.character(nms_euk))>3]
	
	TAX_bac[[i-2]] <- data.frame()
	TAX_euk[[i-2]] <- data.frame()
	for(j in 1:length(nms_bac)){
		xtmp <- log(apply(as.data.frame(bac[bac[,i]==nms_bac[j],9:ncol(bac)]),2,sum)+1)
		TAX_bac[[i-2]] <- rbind(TAX_bac[[i-2]],xtmp)  
	}
	for(j in 1:length(nms_euk)){
		xtmp <- log(apply(as.data.frame(euk[euk[,i]==nms_euk[j],9:ncol(euk)]),2,sum)+1)
		TAX_euk[[i-2]] <- rbind(TAX_euk[[i-2]],xtmp)
	}
	
	TAX_bac[[i-2]] <- TAX_bac[[i-2]][rowSums(TAX_bac[[i-2]])>500,]
	TAX_euk[[i-2]] <- TAX_euk[[i-2]][rowSums(TAX_euk[[i-2]])>500,]
}


par(mfrow=c(3,2),mar=c(2,5,2,2))
for(i in 1:5){
	matplot(scale(t(TAX_bac[[i]]),scale=FALSE),type='l',lty=1,ylim=c(-8,8),ylab='')
	mtext(taxa[i])
}





par(mfrow=c(3,2),mar=c(2,5,2,2))
for(i in 1:5){
	matplot(scale(t(TAX_euk[[i]]),scale=FALSE),type='l',lty=1,ylim=c(-8,8))
}


dat_bac=dat_euk <- list()
for(i in 1:5){
	dat_bac[[i]] <- list(T=ncol(TAX_bac[[i]]),
	#					 p=nrow(TAX_bac[[i]]),
						 Y=TAX_bac[[i]])
	dat_euk[[i]] <- list(T=ncol(TAX_euk[[i]]),
						 p=nrow(TAX_euk[[i]]),
						 Y=TAX_euk[[i]])
}


###################################################################
## FIT MODELS #####################################################
###################################################################
MCMC_bac=MCMC_euk <- list()
POST_bac=POST_euk <- list()


mcmc <- sampling(mod, data=dat_bac[[2]], open_progress=TRUE)
post <- extract(mcmc)

for(i in 5:5){
	MCMC_bac[[i]] <- sampling(mod, data=dat_bac[[i]], open_progress=TRUE)
	MCMC_euk[[i]] <- sampling(mod, data=dat_euk[[i]], open_progress=TRUE)

	POST_bac[[i]] <- extract(MCMC_bac[[i]])
	POST_euk[[i]] <- extract(MCMC_euk[[i]])
}

OPT_bac=OPT_euk <- list()

opt <- o(mod, data=dat_bac[[i]])

##--SAVE/LOAD--#############
#save(MCMC_bac, file='d:/dropbox/working/ecological_stability/fits/MCMC_bac.rdata')
#save(MCMC_euk, file='d:/dropbox/working/ecological_stability/fits/MCMC_euk.rdata')
load('d:/dropbox/working/ecological_stability/fits/MCMC_bac.rdata')
load('d:/dropbox/working/ecological_stability/fits/MCMC_euk.rdata')

####################################################################
## ANALYZE MODELS ##################################################
####################################################################
niter <- 4000

eig_gen_bac=eig_gen_euk <- array(NA, dim=c(niter,2,47))
eig_fam_bac=eig_fam_euk <- array(NA, dim=c(niter,2,33))
#vec_bac=vec_euk         <- array(NA, dim=c(niter,10,47))

post_fam_bac <- extract(MCMC_bac[[3]])
post_gen_bac <- extract(MCMC_bac[[4]])

for(i in 1:niter){
	eig_gen_bac[i,1,] <- sort(Re(eigen(post_gen_bac$PHI[i,,])$values),decreasing=TRUE)
	eig_gen_bac[i,2,] <- sort(Im(eigen(post_gen_bac$PHI[i,,])$values),decreasing=TRUE)
	#eig_gen_euk[i,1,] <- sort(Re(eigen(post_gen_euk$PHI[i,,])$values),decreasing=TRUE)
	#eig_gen_euk[i,2,] <- sort(Im(eigen(post_gen_euk$PHI[i,,])$values),decreasing=TRUE)
	
	eig_fam_bac[i,1,] <- sort(Re(eigen(post_fam_bac$PHI[i,,])$values),decreasing=TRUE)
	eig_fam_bac[i,2,] <- sort(Im(eigen(post_fam_bac$PHI[i,,])$values),decreasing=TRUE)
	#eig_fam_euk[i,1,] <- sort(Re(eigen(post_fam_euk$PHI[i,,])$values),decreasing=TRUE)
	#eig_fam_euk[i,2,] <- sort(Im(eigen(post_fam_euk$PHI[i,,])$values),decreasing=TRUE)

	#vec_bac[i,,] <- Re(eigen(post_bac$PHI[i_samp,,])$vectors)
	#vec_euk[i,,] <- Re(eigen(post_euk$PHI[i_samp,,])$vectors)
}



brks <- seq(0.9,1.2,0.0025)
par(mfrow=c(3,1))
hist(eig_gen_bac[,1,1],freq=FALSE,ylim=c(0,120),col=adjustcolor('blue',alpha.f=0.3),main='')
	lines(density(eig_bac[,1,1]),col=adjustcolor('blue',alpha.f=0.5))
hist(eig_gen_euk[,1,1],freq=FALSE,add=TRUE,breaks=brks,col=adjustcolor('red',alpha.f=0.3))
	lines(density(eig_euk[,1,1]),col=adjustcolor('red',alpha.f=0.5))
hist(eig_fam_bac[,1,1],freq=FALSE,ylim=c(0,120),breaks=brks,col=adjustcolor('blue',alpha.f=0.3),main='')
	lines(density(eig_bac[,1,1]),col=adjustcolor('blue',alpha.f=0.5))
hist(eig_fam_euk[,1,1],freq=FALSE,add=TRUE,breaks=brks,col=adjustcolor('red',alpha.f=0.3))
	lines(density(eig_fam_euk[,1,1]),col=adjustcolor('red',alpha.f=0.5))

##--Raindrop plot across eigenvalues--#######################
##--genera--#########
par(mfrow=c(2,2))
#hist(eig_bac[,1,1],freq=FALSE,ylim=c(0,30),col=adjustcolor('blue',alpha.f=0.3),xlim=c(0,1.2))
hist(-999,freq=FALSE,ylim=c(0,35),col=adjustcolor('blue',alpha.f=0.3),xlim=c(0,1.2),main='',xlab='')
for(i in 1:10){
	lines(density(eig_fam_euk[,1,i]),col=i)
}
hist(-999,freq=FALSE,ylim=c(0,35),col=adjustcolor('blue',alpha.f=0.3),xlim=c(0,1.2),main='',xlab='')
for(i in 1:10){
	lines(density(eig_fam_euk[,1,i]),col=i)
}
##############################################
##--families--################################
##############################################
par(mfrow=c(3,2),mar=c(2,2,2,2),oma=c(2,2,2,2))

matplot(scale(t(TAX_bac[[4]]),scale=FALSE),type='l',lty=1,ylim=c(-5,5),ylab='',bty='n'); mtext(taxa[4],cex=1.2)
	mtext(side=2,line=2.5,'Normalized Abundance',cex=0.7); mtext(side=1,'Days',line=2.5)
matplot(scale(t(TAX_bac[[5]]),scale=FALSE),type='l',lty=1,ylim=c(-5,5),ylab='',bty='n'); mtext(taxa[5],cex=1.2)
	mtext(side=1,'Days',line=2.5)
#hist(eig_bac[,1,1],freq=FALSE,ylim=c(0,30),col=adjustcolor('blue',alpha.f=0.3),xlim=c(0,1.2))
hist(-999,freq=FALSE,ylim=c(0,20),col=adjustcolor('blue',alpha.f=0.3),xlim=c(0,1.2),main='',xlab=''); 
	mtext(side=1,'Eigenvalue',line=2.5)
	for(i in 3:20){
		lines(density(eig_fam_bac[,1,i]))
	}
	mtext(side=2,line=2.5,'Probability Density',cex=0.7)

hist(-999,freq=FALSE,ylim=c(0,20),col=adjustcolor('blue',alpha.f=0.3),xlim=c(0,1.2),main='',xlab='')
	mtext(side=1,'Eigenvalue',line=2.5)
for(i in 3:20){
	lines(density(eig_gen_bac[,1,i]))
}

t <- seq(0,8,0.01)

PHI <- apply(post_fam_bac$PHI,c(2,3),mean)
eig <- eigen(PHI)

plot(t,1*exp(-Re(eig$values[1])*t)*Re(eig$vectors[1,1]),type='l',ylim=c(-0.3,0.3),bty='n',xlab='',ylab=''); mtext(side=1,'Days',line=2.5)
	for(i in 3:20){
		lines(t,exp(-Re(eig$values[i])*t)*Re(eig$vectors[1,i]))
	}
	mtext(side=2,line=2.5,expression('X - X'[0]),cex=0.7)

PHI <- apply(post_gen_bac$PHI,c(2,3),mean)
eig <- eigen(PHI)

plot(t,1*exp(-Re(eig$values[1])*t)*Re(eig$vectors[1,1]),type='l',ylim=c(-0.3,0.3),bty='n',xlab='',ylab=''); mtext(side=1,'Days',line=2.5)
for(i in 3:20){
	lines(t,exp(-Re(eig$values[i])*t)*Re(eig$vectors[1,i]))
}



#####################################################################
## EIGENMODES #######################################################
#####################################################################

par(mfrow=c(1,2))
PHI <- apply(post_fam_bac$PHI,c(2,3),mean)
eig <- eigen(PHI)

plot(t,1*exp(-Re(eig$values[1])*t)*Re(eig$vectors[1,1]),type='l',ylim=c(-0.6,0.6))
for(i in 2:20){
	lines(t,exp(-Re(eig$values[i])*t)*Re(eig$vectors[1,i]))
}

PHI <- apply(post_gen_bac$PHI,c(2,3),mean)
eig <- eigen(PHI)

plot(t,1*exp(-Re(eig$values[1])*t)*Re(eig$vectors[1,1]),type='l',ylim=c(-0.6,0.6))
for(i in 2:20){
	lines(t,exp(-Re(eig$values[i])*t)*Re(eig$vectors[1,i]))
}



