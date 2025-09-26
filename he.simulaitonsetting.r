
om=function(invar){ ############sim func#################
  set.seed(invar)
  param_pairs <- data.frame(
    nn = c(1, 2, 3, 3),
    pp = c(1, 2, 3, 4)
  )
  grid <- merge(
    param_pairs,
    expand.grid(ddad = 1:12, deg_pa=1:3 ,KEEP.OUT.ATTRS = FALSE),
    by = NULL,
    all = TRUE
  )   
  pscore <- function(data,exp){
    library(grf)
    exper=c(rep('rct',6),rep('obv',6))
    data.trans <- data
    expt=exper[exp]
    y <- data.trans[,1]
    g <- data.trans[,2]
    
    if(is.null(exp)){
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
    return(pi)
  }
  ##############definition######################
  f <- unlist(grid[invar, ], use.names = FALSE)
  ##############definition######################
  v=invar
  nn=f[1]
  pp=f[2]
  ddad=f[3]
  deg_pa=f[4]
  noise=0.5
  bs_pa=100
  re_pa=1
  ############==========GEN==========###########

  transOUT<- function(Y,Z,X){
    datq=data.frame(Z,X)
    mdl=regression_forest(X=X,Y=Z)
    XX=data.frame(X)
    Zpre=predict(mdl,newdata=XX)
    Zpre[which(Zpre>0.90)]=0.90
    Zpre[which(Zpre<0.10)]=0.10
    Zpe=Z*Y/Zpre+(1-Z)*(-Y)/(1-Zpre)
    return(Zpe)
  }
  pscore <- function(data,exp){
    exper=c(rep('rct',6),rep('obv',6))
    data.trans <- data
    expt=exper[exp]
    y <- data.trans[,1]
    g <- data.trans[,2]
    
    if(is.null(exp)){
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
    return(pi)
  }
  transOUT.rct<- function(Y,Z,X){
    dat=data.frame(Z,X)
    XX=data.frame(X)
    Zpre=0.5
    Zpe=Z*Y/Zpre+(1-Z)*(-Y)/(1-Zpre)
    return(Zpe)
  }
  pmpv=function(X){
    for (i in 1:length(X)) {
      print(paste0("[[",i,"]]"))
      print(mean(X[[i]]))
      print(var(X[[i]]))
    }
  }
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
    return(X)
  }
  datagen<- function(mu,tau,type,X,noise=1){
    n=nrow(X)
    if (type==1) {
      t=rbinom(n,1,0.5)
    }else{
      t=rbinom(n, 1, 0.3 + 0.2 * (X[, 1] > 1)+ 0.2 * (X[, 3] > 1))
      
    }
    
    y=mu+(t-0.5)*tau+rnorm(n,0,noise)
    
    
    if (type==1) {#response
      dat=data.frame(transOUT.rct(y,t,X),t,tau,X,y)
    }else{
      dat=data.frame(transOUT(y,t,X),t,tau,X,y)
    }
    
    return(dat)
  }
  
  mx=c(3,2,3,2,3,2,3,2,3,2,3,2)
  tx=c(5,6,7,8,9,10,5,6,7,8,9,10)
  nx=c(200,400,1000)
  px=c(25,50,100,150)
  typex=c(rep(1,6),rep(0,6))
  
  
  #=======train data=====================
  X=Xgen(n=nx[nn],p=px[pp])
  fn_tau=fnum_tau(X)
  fn_mu=fnum_mu(X)
  mu0=fn_mu[[mx[ddad]]]
  tau0=fn_tau[[tx[ddad]]]
  
  datalist1=datagen(mu0,tau0,typex[ddad],X,noise)
  addr0=getwd()
  addr_save <- sprintf("/train/scenario%s_n%s_p%s_%s_noise%s.csv",ddad,nx[nn],px[pp],v,noise)
  input_path1 <-paste0(addr0, addr_save )
  write.csv(datalist1,input_path1)
  #write.csv(datalist1,sprintf("E:/SCBM/dat.24.10.24/train/scenario%s_n%s_p%s_%s_noise%s.csv",ddad,nx[nn],px[pp],v,noise),row.names = FALSE)
  #=======test data=====================
  X=Xgen(n=nx[nn],p=px[pp])
  fn_tau=fnum_tau(X)
  fn_mu=fnum_mu(X)
  mu0=fn_mu[[mx[ddad]]]
  tau0=fn_tau[[tx[ddad]]]
  
  datalist2=datagen(mu0,tau0,typex[ddad],X,noise)
  addr0=getwd()
  addr_save <- sprintf("/test/scenario%s_n%s_p%s_%s_noise%s.csv",ddad,nx[nn],px[pp],v,noise)
  input_path2 <-paste0(addr0, addr_save )
  write.csv(datalist2,input_path2)
  #write.csv(datalist2,sprintf("E:/SCBM/dat.24.10.24/test/scenario%s_n%s_p%s_%s_noise%s.csv",ddad,nx[nn],px[pp],v,noise),row.names = FALSE)

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
  res.save <- data.frame(matrix(0,nrow(train),10))
  colnames(res.save) <- c("true","GL","GLPS","GLSTRA","TO","BCM","BCM1","CBART","CF","VT")
  ps=pscore(train,ddad)
  stra=stratify(ps,trainG,min.per.arm =100)[[1]]
  stra1=stratify(ps,trainG)[[1]]
  res.save$true<-testdat$tau
  res.save.train <-res.save
  # TODO: Continue with model fitting and result saving (e.g., causal MARS, GLPS, causal forests)
  #################MARS GL GLPS TO################################### 
  addr_source <- sprintf("/fun_GL_PS.R")
  source_path <-paste0(addr0, addr_source )
  source(source_path)
  traindat=cbind(1,datalist1)
  testdat=cbind(1,datalist2)
  tau=testdat$tau
  tau_train=traindat$tau
  hte=marslasso_GLPSTO(traindat,testdat,ps,deg=deg_pa,boost=bs_pa,per_resamp=re_pa)
  hte_test=hte[[1]]
  hte_test_ps=hte[[2]]
  hte_test_TO=hte[[3]]
  res.save$GL=hte_test
  res.save$GLPS=hte_test_ps
  res.save$TO=hte_test_TO
  print("1")
  #################STRA###################################
  addr_source <- sprintf("/fun_GL_PS_stratify.R")
  source_path <-paste0(addr0, addr_source )
  source(source_path)
  traindat=cbind(1,datalist1)
  testdat=cbind(1,datalist2)
  hte_stra=marslasso_stra(traindat,testdat,ps,deg=deg_pa,boost=bs_pa,per_resamp=re_pa)
  res.save$GLSTRA=hte_stra
  print("1")
  ##################BCM###################################
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
  ##################VT#############################
  teX1=trainX[traindat$t==1,]
  teY1=trainY[traindat$t==1]
  teX0=trainX[traindat$t==0,]
  teY0=trainY[traindat$t==0]
  f1=randomForest(teX1,teY1)
  f0=randomForest(teX0,teY0)
  Yf1=predict(f1,trainX)
  Yf0=predict(f0,trainX)
  hte=Yf1-Yf0
  datvt=data.frame(trainX,hte)
  mdlvt=rpart(hte~.,datvt)
  output1=predict(mdlvt,newdata=testX)
  res.save$VT <- output1
  print("3")
  ##################CF#################################
  #3.causal forest
  c.forest<- causal_forest(trainX,trainY,trainG,W.hat = ps,tune.parameters = "all")
  c.pred <- predict(c.forest, testX)$prediction
  res.save$CF<-c.pred
  ##################BARTc#################################
  bc=bartc(trainY, trainG, trainX,keepTrees = TRUE)
  hte.bart= predict(bc,type="ite",newdata =testX)
  hte.bart=apply(hte.bart,2,mean)
  res.save$CBART <- hte.bart
  ######################save########################
  addr0=getwd()
  addr_save <- sprintf("/out/OUT_v2_%s_n%s_p%s_%s_deg%s.csv",ddad,nx[nn],px[pp],v,deg_pa)
  outpath <-paste0(addr0, addr_save )
  write.csv(res.save,outpath)
  return(NULL)
}







####################################################################################
########simulation start######################################################
####################################################################################
## --- auto-install + load -------------------------------------------------
pkgs <- c("rpart","Rcpp","grf","earth","BART",
          "bartCause","randomForest","causalLearning",
          "doSNOW","foreach","parallel")          # add/remove as needed

need <- pkgs[!pkgs %in% rownames(installed.packages())]   # not yet present?
if (length(need)) install.packages(need, dependencies = TRUE)

## ------------------------------------------------------------------------
install.packages("remotes")
remotes::install_github("saberpowers/causalLearning")
library(doSNOW)
library(doParallel)
library(rpart)
library(Rcpp)
library(grf)
library(earth)
library(BART)
library(bartCause)
library(randomForest)
library(causalLearning)
# brings in foreach, %dopar%
lcore <- max(1, detectCores() - 10)      # keep a couple of cores free
cl    <- makeCluster(lcore)
registerDoSNOW(cl)
# ── R packages each worker will need ─────────────────────────────────────────
pkgs <- c("rpart","Rcpp","grf","earth","BART",
          "bartCause","randomForest","causalLearning")

clusterExport(cl, "pkgs")       
clusterEvalQ(cl, lapply(pkgs, require, character.only = TRUE))

n_tasks <- nrow(grid)          # 4 800 in your case

foreach(r = 1:100,
          .packages      = pkgs ) %dopar% {
          aft = om(r) 
          aft=NULL# any error here stops everything
          }

close(pb)
stopCluster(cl)

