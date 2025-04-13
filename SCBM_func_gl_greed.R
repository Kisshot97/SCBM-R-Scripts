################################################################################
# Module: Func Gl Greed
# Description: This script is part of the SCBM framework and implements
# core functionality related to func gl greed.
################################################################################


# --------------------------------------------------
# Function: marslasso_GLPSTO
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
marslasso_GLPSTO <- function(traindat,testdat,ps,deg,boost,per_resamp){

# --------------------------------------------------
# Function: earthlassog
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
  earthlassog <- function(datain1,TO,deg,maxn,samp){

# ====================
# Load Required Libraries
# ====================
    library("earth")

# ====================
# Load Required Libraries
# ====================
    library(glmnet)
    SAMDAT=sample(nrow(datain1),samp*nrow(datain1))
    g_dat=datain1[SAMDAT,]
    g_oobdat=datain1
    TTO=TO[SAMDAT]
    oobTO=TO
    X=g_dat[,c(-1,-2,-3,-4,-ncol(g_dat))]
    oobX=g_oobdat[,c(-1,-2,-3,-4,-ncol(g_oobdat))]
    ytt=g_dat[,ncol(g_dat)]
    oobytt=g_oobdat[,ncol(g_oobdat)]
# --- Construct a Data Frame ---
    dat0=data.frame(Y=TTO,X)# sampled data
    
    
    rawtauar=numeric()
    rss3=numeric()
    for (d in deg) {
      q0=earth(Y~.,data=dat0,degree=d,thresh=0.001,nk=maxn,pmethod="none")
      bx=q0$bx[,-1]
      assign(paste0("bx",d),bx)
      dir=q0$dirs[-1,]
      assign(paste0("dir",d),dir)
      cut=q0$cuts[-1,]
      assign(paste0("cut",d),cut)
      coe=q0$coefficients[-1,]
      assign(paste0("coe",d),coe)
      XX_test=numeric()
      py1=predict(q0,newdata = oobX)
      rss=(py1-oobTO)^2
      #=====================================
      rss3[d]=sum(rss)
      assign(paste0("q0",d),q0)
    }
    df=which(rss3==min(na.omit(rss3)))[1]
    bx=get(paste0("bx",df))
    dir=get(paste0("dir",df))
    cut=get(paste0("cut",df))
    coe=get(paste0("coe",df))
    q0=get(paste0("q0",df))
    bcd=list(bx,cut,dir,coe,q0)
# --- Return Final Output ---
    return(bcd)
  }

# ====================
# Load Required Libraries
# ====================
  library(glmnet)

# ====================
# Load Required Libraries
# ====================
  library("grpreg")

# ====================
# Load Required Libraries
# ====================
  library(parallel)

# ====================
# Load Required Libraries
# ====================
  library(grf)
  #==============pre-process=====
  datfortrain=traindat
  datfortest=testdat
  #==============================
  oobdat=datfortrain
  #########################
  group=oobdat[,3]###??????
  g1 <- sort(unique(group))
  gn <- length(group)
  g<-numeric(gn)
  g[which(group==g1[1])]=0
  g[which(group==g1[2])]=1
  gt<-numeric(gn)
  gt[which(group==g1[1])]=1
  gt[which(group==g1[2])]=0
  #########################
  oobX=oobdat[,c(-1,-2,-3,-4,-ncol(oobdat))]
  oobytt=oobdat[,ncol(oobdat)]
  oobto=oobdat[,2]
  pscore=g*1/ps+(g-1)*1/(1-ps)
  Zpe=g*oobytt/ps-(1-g)*(-oobytt)/(1-ps)
  y_ps=Zpe
    #============data process===============
  nboots=boost
  maxn=ncol(oobX)
  indat=datfortrain
  bcd=list()
  for (i in 1:nboots) {
    bcd[[i]]=earthlassog(datain=indat,TO=Zpe,deg=deg,maxn = maxn,samp=per_resamp)
  }
  cut0=numeric()
  dir0=numeric()
  coe0=numeric()
  q0=list()
  for (i in 1:nboots) {
    cut0=rbind(cut0,bcd[[i]][[2]])
    dir0=rbind(dir0,bcd[[i]][[3]])
# --- Construct a Data Frame ---
    coe0=rbind(coe0,data.frame(bcd[[i]][[4]]))
    q0[[i]]=bcd[[i]][[5]]
  }
  
 #============group lasso==============#

  XX_test=numeric()
  for (i in 1:nrow(oobX)) {
    xin=matrix(unlist(rep(oobX[i,],nrow(cut0))),nrow = nrow(cut0),byrow = TRUE)-cut0
    dirin=xin*dir0
    a02=which(dir0==2)##   ???
    dirin[which(dirin==0)]=1
    dirin[dirin<0]=0
    dirin[a02]=xin[a02]
# --- Construct a Data Frame ---
    bs=data.frame(apply(dirin,1,prod))#*coe0#*coe0#连???
    bs=unlist(bs)
    XX_test=rbind(XX_test,bs)
  }
  
  ###########################
  xxg=cbind(XX_test*g,XX_test*gt)
  groupx=c(c(1:ncol(XX_test)),c(1:ncol(XX_test)))
  assign(paste0("cv.glmdl",1),cv.grpreg(xxg,oobytt,group=groupx,nfolds=10))
 # assign(paste0("cv.glmdl",2),cv.grpreg(xxg,y_ps,group=groupx,nfolds=10))
 # assign(paste0("cv.glmdl",3),cv.grpreg(xxg,y_ps1,group=groupx,nfolds=10))
  fittau0=NA
  

  #========================================#
  ############test#######################
  #========================================#
  dat_test=datfortest
  dtest=dat_test[,c(-1,-4)]
  dtestx=dtest[,c(-1,-2,-ncol(dtest))]
  gtest=dat_test[,3]
  gtest0 <- ifelse(gtest == 1, 0, 1)
  XX_test1=numeric()
 
  for (i in 1:nrow(dtestx)) {
    xin=matrix(unlist(rep(dtestx[i,],nrow(cut0))),nrow = nrow(cut0),byrow = TRUE)-cut0
    dirin=xin*dir0
    a02=which(dir0==2)##   ???
    dirin[which(dirin==0)]=1
    dirin[dirin<0]=0
    dirin[a02]=xin[a02]
# --- Construct a Data Frame ---
    bs=data.frame(apply(dirin,1,prod))
    bs=unlist(bs)
    XX_test1=rbind(XX_test1,bs)
  }
  px1t=cbind(XX_test1*gtest,XX_test1*gtest0)
  yb1t=predict(cv.glmdl1,px1t)
  rres=yb1t
# --- Return Final Output ---
  return(rres)
  


###################################earthlasso##############################

##################################
####################################
####################################



}