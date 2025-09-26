marslasso_stra <- function(traindat,testdat,ps,VIMP=0,deg=c(1:3),boost=50,per_resamp=0.3){
  library(glmnet)
  library("grpreg")
  library(parallel)
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
  y_ps=oobytt*pscore
  Zpe=g*oobytt/ps+(1-g)*(-oobytt)/(1-ps)
  #============data process===============
  nboots=boost
  maxn=ncol(oobX)
  indat=datfortrain
  bcd=list()
  for (i in 1:nboots) {
    bcd[[i]]=earthlasso(indat,Zpe,deg=deg,maxn = maxn,samp=per_resamp)
  }
  cut0=numeric()
  dir0=numeric()
  coe0=numeric()
  q0=list()
  for (i in 1:nboots) {
    cut0=rbind(cut0,bcd[[i]][[2]])
    dir0=rbind(dir0,bcd[[i]][[3]])
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
    bs=data.frame(apply(dirin,1,pprod))#*coe0#*coe0#连???
    bs=unlist(bs)
    XX_test=rbind(XX_test,bs)
  }
####################stratification################################### 
  n_stra <- 5  # 
  min_stra <- 1  # 
  batch_size <- ceiling(nrow(XX_test) / n_stra)
  od_XX_test <- XX_test[order(ps),]
  od_y <- oobytt[order(ps)]
  od_g <- g[order(ps)]
  XX_test_group <- cut(seq_along(ps), breaks = 4, labels = FALSE)
  min_01=min(sum(g==0),sum(g==1))
  batches=list()
  g_batch=list()
  y_batch=list()
  for (i in 1:max(XX_test_group)) {
    batches[[i]]=od_XX_test[which(XX_test_group==i),]
    g_batch[[i]]=od_g[which(XX_test_group==i)]
    y_batch[[i]]=od_y[which(XX_test_group==i)]
  }
  if (min_01 < min_stra) stop("Error: min_stra over")
  
  ###########################
  for (i in 1:max(XX_test_group)) {  
    groupx=c(c(1:ncol(batches[[i]])),c(1:ncol(batches[[i]])))
    
    xxg=cbind(batches[[i]]*g_batch[[i]],batches[[i]]*(1-g_batch[[i]]))
    assign(paste0("cv.glmdl",i),cv.grpreg(xxg,y_batch[[i]],group=groupx,nfolds=10))
  }
  
  

  #========================================#
  ############test#######################
  #========================================#
  dat_test=datfortest
  dtest=dat_test[,c(-1,-4)]
  dtestx=dtest[,c(-1,-2,-ncol(dtest))]
  XX_test1=numeric()

  for (i in 1:nrow(dtestx)) {
    xin=matrix(unlist(rep(dtestx[i,],nrow(cut0))),nrow = nrow(cut0),byrow = TRUE)-cut0
    dirin=xin*dir0
    a02=which(dir0==2)##   ???
    dirin[which(dirin==0)]=1
    dirin[dirin<0]=0
    dirin[a02]=xin[a02]
    bs=data.frame(apply(dirin,1,pprod))#*coe0#*coe0#连???
    bs=unlist(bs)
    XX_test1=rbind(XX_test1,bs)
  }
  px1t=cbind(XX_test1,XX_test1*0)
  px2t=cbind(XX_test1*0,XX_test1)
  output=data.frame(matrix(ncol =max(XX_test_group), nrow =nrow(datfortest)))
  ybft=matrix(nrow = nrow(XX_test1),ncol = max(XX_test_group))
  for (i in 1:max(XX_test_group)) {  
    yb1t=predict(get(paste0("cv.glmdl",i)),px1t)
    yb2t=predict(get(paste0("cv.glmdl",i)),px2t)
    ybft[,i]=yb1t-yb2t
    #assign(paste0("ybft",i),ybft)
  }
  y_stra=apply(ybft, 1, mean)

  
#====================================#
#############VIMP#####################
#====================================#
  VIMP_PLOT=list()
if (VIMP==1){
  cl <- makeCluster(14)
  p=c(1:ncol(X))#ncol(X)
  bcm=list()
  bcm=parLapply(cl,p,fitfun,oobX,oobytt,cut0,dir0,cv.glmdl1,g,gt)
  stopCluster(cl)
  
  #=====================================
  
  
  fitm=numeric()
  for (i in 1:ncol(X)) {
    
    xxg2t=cbind(bcm[[i]]*1,bcm[[i]]*0)
    xxg2tt=cbind(bcm[[i]]*0,bcm[[i]]*1)
    aa1=predict(cv.glmdl1,X=xxg2t)
    aa0=predict(cv.glmdl1,X=xxg2tt)
    aatau=aa1-aa0
    
    fitm=rbind(fitm,aatau)
  }
 
  fitminu=matrix(rep(fittau0,ncol(X)),ncol = length(fittau0),byrow = TRUE)
  errm=apply((fitm-fitminu)^2, 1,sum)

  reserr=normalize(errm,method= "range",range = c(0,100))
  
  names(reserr)=c(sprintf("x%s",1:ncol(X)))
  VIMP_PLOT=barplot(sort(reserr[which(reserr>0)],decreasing=T), main = "variable importance")

  
}else{
  reserr=NULL
}
  rres=y_stra
  return(rres)
  
}
## prod  连???######################
pprod <- function(vec){
  mul= prod(vec)
  return(mul)
}
## 2値か多値かを判断する関数#######
Measure <- function(vec){
  flug=numeric()
  uni <- unique(vec)
  if (length(uni)==2){ flug <- 0}
  if (length(uni)>2) { flug <- 1}
  return(flug)
}
#============================
fitfun <- function(p,oobX,oobytt,cut0,dir0,cv.glmdl1,g,gt){
  ## prod  连???######################
  pprod <- function(vec){
    mul= prod(vec)
    return(mul)
  }
  #=================================
  ec0=which(dir0[,p]%in%0==FALSE)
  XX_test2=numeric()
  for (i in 1:nrow(oobX)) {
    xin=matrix(unlist(rep(oobX[i,],nrow(cut0))),nrow = nrow(cut0),byrow = TRUE)-cut0
    dirin=xin*dir0
    a02=which(dir0==2)##   ???
    dirin[which(dirin==0)]=1
    dirin[dirin<0]=0
    dirin[a02]=xin[a02]
    dirin[ec0,p]=0
    bs=apply(dirin,1,pprod)#连???
    XX_test2=rbind(XX_test2,bs)
  }
  #xxg2=cbind(XX_test2*g,XX_test2*gt)
  #fit1=predict.grpreg(cv.glmdl1,X=xxg2)
  #res=(oobytt-fit1)^2
  return(XX_test2)
}
###################################earthlasso##############################
earthlasso <- function(datain,TO,btype=1,deg,maxn,samp){
  pprod <- function(vec){#liancheng
    mul= prod(vec)
    return(mul)
  }
  
  ############start####################
  
  library("earth")
  library(glmnet)
  DAT=sample(nrow(datain),samp*length(datain),replace = TRUE)
  data1=datain[DAT,]
  data2=datain#[-DAT,]
  TTO=TO[DAT]
  oobTO=TO
  dat=data1
  oobdat=data2
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
  X=dat[,c(-1,-2,-3,-4,-ncol(dat))]
  oobX=oobdat[,c(-1,-2,-3,-4,-ncol(oobdat))]
  ytt=dat[,ncol(dat)]
  oobytt=oobdat[,ncol(oobdat)]
  Z=dat[,3]
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
  return(bcd)
}
##################################
####################################
####################################











