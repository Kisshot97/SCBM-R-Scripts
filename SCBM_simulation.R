################################################################################
# Module: Simulation
# Description: This script is part of the SCBM framework and implements
# core functionality related to simulation.
################################################################################


# ====================
# Load Required Libraries
# ====================
library(parallel)
set.seed(10086)
om=function(f){ ############sim func#################
 

# --------------------------------------------------
# Function: pscore
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
  pscore <- function(data,exp){

# ====================
# Load Required Libraries
# ====================
    library(grf)
    exper=c(rep('rct',6),rep('obv',6))
    data.trans <- data
    expt=exper[exp]
    y <- data.trans[,1]
    g <- data.trans[,2]
    
    if(is.null(exp)){
# --- Return Final Output ---
      return(data)
    }
    
    if(expt == "rct"){
      pi <- rep(0.5,nrow(data))
    }
    if(expt == "obv"){
      # Linearized probit regression PS:
      m.out1 <- predict(regression_forest(data.trans[,c(-1,-2)],data.trans[,2]))
      dis <- m.out1$predictions
      dis[dis > 0.90] <- 0.90
      dis[dis < 0.10] <- 0.10
      pi <- dis 
    }
# --- Return Final Output ---
    return(pi)
  }
  ##############definition######################
  v=f[3]
  nn=f[1]
  pp=f[2]
  ddad=f[4]
  genen=1
  noise=0.5
  ############==========GEN==========###########

# ====================
# Load Required Libraries
# ====================
  library(grf)

# --------------------------------------------------
# Function: transOUT
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
  transOUT<- function(Y,Z,X){
# --- Construct a Data Frame ---
    datq=data.frame(Z,X)
    mdl=regression_forest(X=X,Y=Z)
# --- Construct a Data Frame ---
    XX=data.frame(X)
    Zpre=predict(mdl,newdata=XX)
    Zpre[which(Zpre>0.90)]=0.90
    Zpre[which(Zpre<0.10)]=0.10
    Zpe=Z*Y/Zpre+(1-Z)*(-Y)/(1-Zpre)
# --- Return Final Output ---
    return(Zpe)
  }

# --------------------------------------------------
# Function: transOUT.rct
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
  transOUT.rct<- function(Y,Z,X){
# --- Construct a Data Frame ---
    dat=data.frame(Z,X)
# --- Construct a Data Frame ---
    XX=data.frame(X)
    Zpre=0.5
    Zpe=Z*Y/Zpre+(1-Z)*(-Y)/(1-Zpre)
# --- Return Final Output ---
    return(Zpe)
  }
  pmpv=function(X){
    for (i in 1:length(X)) {
      print(paste0("[[",i,"]]"))
      print(mean(X[[i]]))
      print(var(X[[i]]))
    }
  }

# --------------------------------------------------
# Function: fnum_tau
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
  fnum_tau<- function(X){
    # mx=c(3,2,3,2,3,2,3,2,3,2,3,2)
    # tx=c(5,6,7,8,9,10,5,6,7,8,9,10)
    
    f1 <- rep(0,nrow(X))
    
    f2 <- -3*X[,1]-3*X[,3]
    
    f3 <- f2
    
    f4 <- 0
    
    f5 <- 0
    
    f6 <- X[,1]+2*X[,3]+X[,5]+2*X[,7]+X[,9]#new
    
    f7 <- 5*as.numeric(X[,1]>1)*as.numeric(X[,7] > 0) + 5*as.numeric(X[,5] > 1)-5#new
    
    f8 <- X[,1]*X[,5] + 2*X[,3]*X[,7] + 2*X[,2]*X[,4] -2 #new
    
    f9 <- 2*sin(pi*X[,1]*X[,7]) + 2*(X[,3]+0.5)^2 + X[,9] -2 #new
    
    f10 <- 1/(1+exp(-X[,1]))+2*X[,3]^2+2*X[,9]-2#new
    
    f11 <- 0
    
    return (list(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11))
  }
  
  
  
  fnum_mu<- fnum_tau
  

# --------------------------------------------------
# Function: Xgen
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
  Xgen<- function(n,p){
    di=1:p
    n=n
    ji <- seq(from=1,to=length(di),by=2)
    ou <- setdiff(di,ji)
    X=matrix(numeric(n*length(di)),nrow = n)
    for (i in 1:length(di)) {
      if(i%in%ji){
        X[,i]=rnorm(n)
      }
      else{
        X[,i]=rbinom(n,1,0.5)
      }
    }
# --- Return Final Output ---
    return(X)
  }

# --------------------------------------------------
# Function: datagen
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
  datagen<- function(mu,tau,type,X,noise=1){
    n=nrow(X)
    if (type==1) {
      t=rbinom(n,1,0.5)
    }else{
      t=rbinom(n, 1, 0.3 + 0.2 * (X[, 1] > 1)+ 0.2 * (X[, 3] > 1))
      
    }
    
    y=mu+(t-0.5)*tau+rnorm(n,0,noise)
    
    
    if (type==1) {#response
# --- Construct a Data Frame ---
      dat=data.frame(transOUT.rct(y,t,X),t,tau,X,y)
    }else{
# --- Construct a Data Frame ---
      dat=data.frame(transOUT(y,t,X),t,tau,X,y)
    }
    
# --- Return Final Output ---
    return(dat)
  }
 
  mx=c(3,2,3,2,3,2,3,2,3,2,3,2)
  tx=c(5,6,7,8,9,10,5,6,7,8,9,10)
  nx=c(200,400,1000)
  px=c(25,50,100,150)
  typex=c(rep(1,6),rep(0,6))
  
  if (genen==1) {
    #=======train data=====================
    X=Xgen(n=nx[nn],p=px[pp])
    fn_tau=fnum_tau(X)
    fn_mu=fnum_mu(X)
    mu0=fn_mu[[mx[ddad]]]
    tau0=fn_tau[[tx[ddad]]]
    
    datalist1=datagen(mu0,tau0,typex[ddad],X,noise)

# --- Save Results to CSV ---
    write.csv(datalist1,sprintf(file.path(getwd(), "scenario%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise),row.names = FALSE)
    #=======test data=====================
    X=Xgen(n=nx[nn],p=px[pp])
    fn_tau=fnum_tau(X)
    fn_mu=fnum_mu(X)
    mu0=fn_mu[[mx[ddad]]]
    tau0=fn_tau[[tx[ddad]]]
    
    datalist2=datagen(mu0,tau0,typex[ddad],X,noise)

# --- Save Results to CSV ---
    write.csv(datalist2,sprintf(file.path(getwd(), "scenario%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise),row.names = FALSE)
    #############read################
    
  }else{

# --- Load CSV Data ---
    datalist1=read.csv(sprintf(file.path(getwd(), "scenario%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))

# --- Load CSV Data ---
    datalist2=read.csv(sprintf(file.path(getwd(), "scenario%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))
  }
  #======================================
  traindat=datalist1
  testdat=datalist2
  trainX=traindat[,c(-1,-2,-3,-ncol(traindat))]
  trainY=traindat[,ncol(traindat)]
  trainG=traindat[,2]
  train=cbind(trainY,trainG,trainX)
  testX= testdat[,c(-1,-2,-3,-ncol(traindat))]
  testY= testdat[,ncol(traindat)]
  testG= testdat[,2]
  test=cbind( testY, testG, testX)

  ######################################
# --- Construct a Data Frame ---
  res.save <- data.frame(matrix(0,nrow(train),9))
  colnames(res.save) <- c("true","GL","GLPS","GLSTRA","TO","BCM","BCM1","CF","VT")
  ps=pscore(train,ddad)
  stra=stratify(ps,trainG,min.per.arm =100)[[1]]
  stra1=stratify(ps,trainG)[[1]]
  res.save$true<-testdat$tau
  res.save.train <-res.save
  
  # TODO: Continue with model fitting and result saving (e.g., causal MARS, GLPS, causal forests)
  #################MARS GL GLPS TO###################################  
  source(file.path(getwd(), "fun_GL_PS.R"))
  traindat=cbind(1,datalist1)
  testdat=cbind(1,datalist2)
  tau=testdat$tau
  tau_train=traindat$tau
  hte=marslasso_GLPSTO(traindat,testdat,ps,deg=2,boost=40,per_resamp=0.5)
  hte_test=hte[[1]]
  hte_test_ps=hte[[2]]
  hte_test_TO=hte[[3]]
  res.save$GL=hte_test
  res.save$GLPS=hte_test_ps
  res.save$TO=hte_test_TO
  print("1")

  #################STRA###################################
  source(file.path(getwd(), "fun_GL_PS_stratify.R"))
  traindat=cbind(1,datalist1)
  testdat=cbind(1,datalist2)
  hte_stra=marslasso_stra(traindat,testdat,ps,deg=2,boost=40,per_resamp=0.5)
  res.save$GLSTRA=hte_stra
  print("1")
  
  # #################BCM###################################
  if (ddad<7) {
    fit_cm <- bagged.causalMARS(trainX,trainG,trainY)
    fit_cm1 <- fit_cm 
  }else{
  fit_cm <- bagged.causalMARS(trainX,trainG,trainY,propensity=TRUE,stratum = stra)
  fit_cm1 <- bagged.causalMARS(trainX,trainG,trainY,propensity=TRUE,stratum = stra1)
  }
  pred_bcm <- predict(fit_cm, newx =testX)
  res.save$BCM<- pred_bcm
  pred_bcm1 <- predict(fit_cm1, newx =testX)
  res.save$BCM1<- pred_bcm1
  
  ###################4.VT#############################
  
  teX1=trainX[traindat$t==1,]
  teY1=trainY[traindat$t==1]
  teX0=trainX[traindat$t==0,]
  teY0=trainY[traindat$t==0]
  f1=randomForest(teX1,teY1)
  f0=randomForest(teX0,teY0)
  Yf1=predict(f1,trainX)
  Yf0=predict(f0,trainX)
  hte=Yf1-Yf0
# --- Construct a Data Frame ---
  datvt=data.frame(trainX,hte)
  mdlvt=rpart(hte~.,datvt)
  output1=predict(mdlvt,newdata=testX)
  res.save$VT <- output1
  print("3")
  ###################################################
  #3.causal forest
  c.forest<- causal_forest(trainX,trainY,trainG,W.hat = ps,tune.parameters = "all")
  c.pred <- predict(c.forest, testX)$prediction
  res.save$CF<-c.pred
  ######################save########################

# --- Save Results to CSV ---
  write.csv(res.save,sprintf(file.path(getwd(), "scenario_%s_n%s_p%s_%s_noise%s.csv"),ddad,nx[nn],px[pp],v,noise))
# --- Return Final Output ---
  return(res.save)
}







#####################
########simulation start############
#####################

lcore=detectCores()-2
cl <- makeCluster(lcore)
clusterEvalQ(cl, {

# ====================
# Load Required Libraries
# ====================
  library(rpart)

# ====================
# Load Required Libraries
# ====================
  library(Rcpp)

# ====================
# Load Required Libraries
# ====================
  library(grf)

# ====================
# Load Required Libraries
# ====================
  library(earth)

# ====================
# Load Required Libraries
# ====================
  library(BART)

# ====================
# Load Required Libraries
# ====================
  library(bartCause)

# ====================
# Load Required Libraries
# ====================
  library(randomForest)

# ====================
# Load Required Libraries
# ====================
  library(causalLearning)
})
#==========================
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
res=list()
n=1

for (i in 1:ncol(npnp)) {
    for (v in c(1:100)) {
    ddad=npnp[1,i]
    nn=npnp[2,i]
    pp=npnp[3,i]
    res[[n]]=c(nn,pp,v,ddad)
    n=n+1
    }
}

t=n=1
for (qwe in seq(1,length(res),25)) {
  rres=list()
  for (i in 1:25) {
    rres[[i]]=res[[n]]
    n=n+1}
  print(rres[[1]])
  t=1
  #============================
  repeat {
    t=t+1
    rply=try(parSapply(cl,rres,om),silent = FALSE)
    print(t)
    if(class(rply)[1]!= "try-error") break
    if (t>5) break
  }
  
}
stopCluster(cl)