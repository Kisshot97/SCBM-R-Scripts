################################################################################
# Module: Make Plot Data
# Description: This script is part of the SCBM framework and implements
# core functionality related to make plot data.
################################################################################

ferr=function(v){
  d=list()
  noise=0.5
  a=1
  nn0=c(1,2,3)
  pp0=c(1,2,3,4)
  npnp=numeric()
  nn=rep(nn0,length(pp0))
  pp=rep(pp0,length(nn0))
  npnp0=rbind(nn,pp)[,-c(4:11)]
  ddad0=c(1:12)
  ddad=rep(ddad0,each=ncol(npnp0))
  npnp=matrix(rep(npnp0, length(ddad0)), nrow = nrow(npnp0), byrow = FALSE)
  npnp=rbind(ddad,npnp)
  propmse=prop_psmse=prop_STRAmse=TOmse= BCMmse= BCM1mse=VTmse=CFmse=numeric()
  propbias=prop_psbias=prop_STRAbias=TObias= BCMbias= BCM1bias=VTbias=CFbias=numeric()
  for (np in 1:ncol(npnp)) {
    ddad=as.numeric(npnp[1,np])
    nn=as.numeric(npnp[2,np])
    pp=as.numeric(npnp[3,np])
    nx=c(200,400,1000)
    px=c(25,50,100,150)
    #GL

# --- Load CSV Data ---
    tau1=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,3]
    #GL ps

# --- Load CSV Data ---
    tau2=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,4]
    #GL STRA

# --- Load CSV Data ---
    tau3=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,5]
    #to

# --- Load CSV Data ---
    tau4=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,6]
    #BCM

# --- Load CSV Data ---
    tau5=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,7]
    #BCM1

# --- Load CSV Data ---
    tau6=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,8]
    #CF

# --- Load CSV Data ---
    tau7=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,9]
    #vt

# --- Load CSV Data ---
    tau8=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,10]
    #true tau

# --- Load CSV Data ---
    tau0=read.csv(sprintf(paste0(getwd(), "/out/scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))[,2]
    
    
    #scaled MSE
    propmse=c(propmse,mean((tau1-tau0)^2))
    prop_psmse=c(prop_psmse,mean((tau2-tau0)^2))
    prop_STRAmse=c(prop_STRAmse,mean((tau3-tau0)^2))
    TOmse=c(TOmse,mean((tau4-tau0)^2))
    BCMmse=c(BCMmse,mean((tau5-tau0)^2))
    BCM1mse=c(BCM1mse,mean((tau6-tau0)^2))
    CFmse=c(CFmse,mean((tau7-tau0)^2))
    VTmse=c(VTmse,mean((tau8-tau0)^2))

    #prop_singl_nonadaparametermse=c(prop_singl_nonadaparametermse,mean((tau7-tau0)^2))
    #noglmse=c(noglmse,mean((tau8-tau0)^2))
    #BIAS
    propbias=c(propbias,mean((tau1-tau0)))
    prop_psbias=c(prop_psbias,mean((tau2-tau0)))
    prop_STRAbias=c(prop_STRAbias,mean((tau3-tau0)))
    TObias=c(TObias,mean((tau4-tau0)))
    BCMbias=c(BCMbias,mean((tau5-tau0)))
    BCM1bias=c(BCM1bias,mean((tau6-tau0)))
    CFbias=c(CFbias,mean((tau7-tau0)))
    VTbias=c(VTbias,mean((tau8-tau0)))

  }
  
  erlist=list(propmse,prop_psmse,prop_STRAmse,TOmse,BCMmse,BCM1mse,CFmse,VTmse)
  biaslist=list(propbias,prop_psbias,prop_STRAbias,TObias,BCMbias,BCM1bias,CFbias,VTbias)

  re=list()
  re[[1]]=biaslist
  re[[2]]=erlist
# --- Return Final Output ---
  return(re)
}

#===============MAIN============================

# ====================
# Load Required Libraries
# ====================
library(doParallel)

# ====================
# Load Required Libraries
# ====================
library(foreach)
Sys.sleep(0)
cl <- makeCluster(detectCores()-1) # create a cluster with n-1 cores
registerDoParallel(cl)
result <-foreach(i = 1:100, .combine = 'rbind') %dopar% {
  tryCatch(ferr(i), error = function(e) NA)
}
stopCluster(cl)
for (i in 1:8) {
  aaaf1=do.call(rbind, lapply(result[,1], "[", i))
  aaaf2=do.call(rbind, lapply(result[,2], "[", i))

# --- Save Results to CSV ---
  write.csv(t(sapply(aaaf1, unlist)),sprintf(paste0(getwd(), "/pdat/bias1_res%s.csv"),i))

# --- Save Results to CSV ---
  write.csv(t(sapply(aaaf2, unlist)),sprintf(paste0(getwd(), "/pdat/mse1_res%s.csv"),i))
}