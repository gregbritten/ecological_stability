library(rstan)
################################################################
## COMPILE MODEL ###############################################
################################################################
options(mc.cores=parallel::detectCores())
mod <- stan_model('d:/dropbox/working/ecological_stability/r_code/MVAR1.stan')

################################################################
## READ DATA ###################################################
################################################################
euk <- read.csv('d:/dropbox/data/martin_platero_etal_2018/csvs/eukaryotic_OTUs.csv')
bac <- read.csv('d:/dropbox/data/martin_platero_etal_2018/csvs/bacterial_OTUs.csv')
env <- read.csv('d:/dropbox/data/martin_platero_etal_2018/csvs/environment.csv')

##--AGGREGATE BY GENUS--##
gen_bac <- unique(bac[,7])   
GEN_bac <- data.frame()      
for(i in 1:length(gen_bac)){
  xtmp     <- apply(as.data.frame(bac[bac[,7]==gen_bac[i],9:ncol(bac)]),2,sum) 
  GEN_bac  <- rbind(GEN_bac,xtmp)                                                #attach as rows to the empty data frame
}

gen_euk <- unique(euk[,7])
GEN_euk <- data.frame()
for(i in 1:length(gen_euk)){
  xtmp <- apply(as.data.frame(euk[euk[,7]==gen_euk[i],9:ncol(euk)]),2,sum)
  GEN_euk <- rbind(GEN_euk,xtmp)
}

dat_gen_bac <- list(T=ncol(GEN_bac),
                    p=nrow(GEN_bac),
                    Y=log(GEN_bac))
dat_gen_euk <- list(T=ncol(GEN_euk),
                    p=nrow(GEN_euk),
                    Y=log(GEN_euk))

par(mfrow=c(2,1),mar=c(2,5,2,2))
matplot(as.matrix(t(dat_gen_bac$Y[1:10,])),type='l')
matplot(as.matrix(t(dat_gen_euk$Y[1:10,])),type='l')

mcmc_bac <- sampling(mod, data=dat_gen_bac)

