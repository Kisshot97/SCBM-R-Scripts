ferr=function(v,sim,deg){
  d=list()
  a=1
  param_pairs <- data.frame(
    nn = c(1, 2, 3, 3),
    pp = c(1, 2, 3, 4)
  )
  
  grid <- merge(
    param_pairs,
    expand.grid(ddad = 1:12, KEEP.OUT.ATTRS = FALSE),
    by = NULL,
    all = TRUE
  )       
  npnp=t(grid)
  propmse=prop_psmse=prop_STRAmse=TOmse= BCMmse= BCM1mse=CBARTmse=VTmse=CFmse=numeric()
  propbias=prop_psbias=prop_STRAbias=TObias= BCMbias= BCM1bias=CBARTbias=VTbias=CFbias=numeric()
  for (np in 1:ncol(npnp)) {
    ddad=as.numeric(npnp[3,np])
    nn=as.numeric(npnp[1,np])
    pp=as.numeric(npnp[2,np])
    nx=c(200,400,1000)
    px=c(25,50,100,150)
    #GL
    tau1=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,3]
    #GL ps
    tau2=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,4]
    #GL STRA
    tau3=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,5]
    #to
    tau4=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,6]
    #BCM
    tau5=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,7]
    #BCM1
    tau6=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,8]
    #BART
    tau7=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,9]
    #CF
    tau8=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,10]
    #vt
    tau9=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,11]
    #true tau
    tau0=read.csv(sprintf("D:/personal/work/test326/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg))[,2]
    #scaled MSE
    propmse=c(propmse,mean((tau1-tau0)^2))
    prop_psmse=c(prop_psmse,mean((tau2-tau0)^2))
    prop_STRAmse=c(prop_STRAmse,mean((tau3-tau0)^2))
    TOmse=c(TOmse,mean((tau4-tau0)^2))
    BCMmse=c(BCMmse,mean((tau5-tau0)^2))
    BCM1mse=c(BCM1mse,mean((tau6-tau0)^2))
    CBARTmse=c(CBARTmse,mean((tau7-tau0)^2))
    CFmse=c(CFmse,mean((tau8-tau0)^2))
    VTmse=c(VTmse,mean((tau9-tau0)^2))
    #BIAS
    propbias=c(propbias,mean((tau1-tau0)))
    prop_psbias=c(prop_psbias,mean((tau2-tau0)))
    prop_STRAbias=c(prop_STRAbias,mean((tau3-tau0)))
    TObias=c(TObias,mean((tau4-tau0)))
    BCMbias=c(BCMbias,mean((tau5-tau0)))
    BCM1bias=c(BCM1bias,mean((tau6-tau0)))
    CBARTbias=c(CBARTbias,mean((tau7-tau0)))
    CFbias=c(CFbias,mean((tau8-tau0)))
    VTbias=c(VTbias,mean((tau9-tau0)))
    
  }
  
  erlist=list(propmse,prop_psmse,prop_STRAmse,TOmse,BCMmse,BCM1mse,CBARTmse,CFmse,VTmse)
  biaslist=list(propbias,prop_psbias,prop_STRAbias,TObias,BCMbias,BCM1bias,CBARTbias,CFbias,VTbias)
  
  re=list()
  re[[1]]=biaslist
  re[[2]]=erlist
  return(re)
}

#===============MAIN============================
library(doParallel)
library(foreach)
for (sim in 1) {
  for (deg in 1:3) {
    cl <- makeCluster(detectCores()-1) # create a cluster with n-1 cores
    registerDoParallel(cl)
    result <-foreach(i = 1:100, .combine = 'rbind') %dopar% {
      tryCatch(ferr(i,sim,deg), error = function(e) NA)
    }
    stopCluster(cl)
    for (i in 1:9) {
      aaaf1=do.call(rbind, lapply(result[,1], "[", i))
      aaaf2=do.call(rbind, lapply(result[,2], "[", i))
      af1=do.call(rbind,aaaf1)
      af2=do.call(rbind,aaaf2)
      write.csv(t(af1), sprintf("D:/personal/work/test326/out/pdat/bias%s_deg%s_var%s.csv", sim, deg, i))
      write.csv(t(af2), sprintf("D:/personal/work/test326/out/pdat/mse%s_deg%s_var%s.csv", sim, deg, i))
    }
  }
  
}