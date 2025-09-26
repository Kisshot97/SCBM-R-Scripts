library(BART)
library(glmnet)
library("grpreg")
library(parallel)
library(grf)
library(foreach)
library(doParallel)
set.seed(99)
for (pc in 1) {
  

#############data input################
data(ACTG175)
ACTG=ACTG175
new_data <- ACTG[, c("pidnum","arms","age", "wtkg", "hemo", "homo", "drugs",  "karnof" , "preanti", "race", "gender", "symptom", "cd40", "cd80", "cd420")]
new_data <- subset(new_data, arms == 0 | arms == 1)
group=new_data[,2]###
g1 <- sort(unique(group))
gn <- length(group)
g<-numeric(gn)
g[which(group==g1[1])]=0
g[which(group==g1[2])]=1
#############data process############
X=new_data[,c(-1,-2,-ncol(new_data))]
cd40=new_data$cd40
cd420=new_data$cd420
quan=function(x) {
  ul=quantile(x, c(0.05,0.95))
  low=which(x < ul[1])
  up=which(x > ul[2])
  x[low]=ul[1]
  x[up]=ul[2]
  return(x)}
cd40=quan(cd40)
cd420=quan(cd420)
ytt=(cd420-cd40)/cd40
datfortrain=data.frame(ytt,g,X)
ps=0.5
Zpe=g*ytt/ps+(1-g)*(-ytt)/(1-ps)
y_ps=Zpe
psFROEST=regression_forest(X,g)
nboots=100
maxn=ncol(X)
per_resamp=1
indat=datfortrain
bcd=list()
deg=2

###################loop####################

vimp_ma=matrix(0,10,12)
for (r in 1:10) {
  vimp=numeric(ncol(X))
  cl <- makeCluster(25)  
  registerDoParallel(cl)
  
  # Parallel execution with foreach
  bcd <- foreach(i = 1:nboots) %dopar% {
    earthlasso(indat, Zpe, deg = deg, maxn = maxn, samp = per_resamp)
  }
  
  # Stop the cluster
stopCluster(cl)
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
for (i in 1:nrow(X)) {
  xin=matrix(unlist(rep(X[i,],nrow(cut0))),nrow = nrow(cut0),byrow = TRUE)-cut0
  dirin=xin*dir0
  a02=which(dir0==2)##
  dirin[which(dirin==0)]=1
  dirin[dirin<0]=0
  dirin[a02]=xin[a02]
  bs=data.frame(apply(dirin,1,prod))#*coe0#*coe0#
  bs=unlist(bs)
  XX_test=rbind(XX_test,bs)
}
xxg=cbind(XX_test*g,XX_test*(1-g))
groupx=c(c(1:ncol(XX_test)),c(1:ncol(XX_test)))
assign(paste0("cv.glmdl",1),cv.grpreg(xxg,ytt,group=groupx,nfolds=10))
#assign(paste0("cv.glmdl",2),cv.grpreg(xxg,y_ps,group=groupx,nfolds=10))
#assign(paste0("cv.glmdl",3),cv.grpreg(xxg,y_ps1,group=groupx,nfolds=10))
fittau0_err=predict(cv.glmdl1,X=xxg)
fittau_err=abs(ytt-fittau0_err)

fitfun <- function(p,oobX,cut0,dir0){
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
    bs=apply(dirin,1,prod)#连???
    XX_test2=rbind(XX_test2,bs)
  }
  return(XX_test2)
}
#########################vimp#################################

cl <- makeCluster(24)
p=c(1:ncol(X))#ncol(X)
bcm=list()
bcm=parLapply(cl,p,fitfun,X,cut0,dir0)
stopCluster(cl)
fitm=numeric()
for (i in 1:ncol(X)) {
  xxg2t=cbind(bcm[[i]]*g,bcm[[i]]*(1-g))
  #xxg2tt=cbind(bcm[[i]]*0,bcm[[i]]*1)
  aa1=predict(cv.glmdl1,X=xxg2t)
  #aa0=predict(cv.glmdl1,X=xxg2tt)
  #aatau=aa1-aa0
  aatau_err=abs(ytt-aa1)
  fitm=rbind(fitm,aatau_err)
}
fitminu=matrix(rep(fittau_err,ncol(X)),ncol = length(fittau_err),byrow = TRUE)
errm=apply(abs(fitm-fitminu), 1,sum)
reserr=errm
names(reserr)=colnames(X)
vimp=rbind(vimp,reserr)
print(paste0("i=",i))
vimp=vimp[-1,]
vimp_ma[r,] <- vimp
}
vimp_means=apply(vimp_ma, 2, mean)
names(vimp_means)=colnames(X)
vimp_means=BBmisc::normalize(vimp_means,method= "range",range = c(0,100))
sorted_indices <- order(vimp_means, decreasing = TRUE)
vimp_sorted <- vimp_means[sorted_indices]


#VIMP_PLOT=barplot(sort(reserr[which(reserr>0)],decreasing=T), main = "variable importance")
#######################################PD###################################################

library(ggplot2)
library(dplyr)
library(tidyr)

data <- data.frame(
  Variable = names(vimp_sorted),
  Importance = vimp_sorted
)
# 按 Importance 排序并重新设置因子顺序
data <- data[order(-data$Importance), ]
data$Variable <- factor(data$Variable, levels = data$Variable)
# 使用 ggplot 作图
VIMP_PLOT=ggplot(data, aes(x = Variable, y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +  # 条形图
  scale_fill_gradient(low = "#87A922", high = "#0B60B0") +  # 定义过渡色
  labs(
    title = "Variable Importance",
    x = "Variables",
    y = "Importance",
    fill = "Importance"
  ) +
  # 最小化主题
  theme(
    legend.position = "none",
    text = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
  )

# Optionally save the grid plot
ggsave(
  filename = paste0("D:/personal/work/test326/ACTGout/pic/", "He_VIMPLOT", ".eps"), # 修改扩展名为 .eps
  plot = VIMP_PLOT, # 图像对象
  width = 5, # 图像宽度
  height = 3, # 图像高度
  device = cairo_ps # 使用 Cairo PostScript 设备
)


}

############################################
############################################
############################################
earthlasso <- function(datain,TO,btype=1,deg,maxn,samp){

  ############start####################
  library("earth")
  library(glmnet)
  DAT=sample(nrow(datain),round(samp*nrow(datain)),replace = TRUE)
  dat=datain[DAT,]
  oobdat=datain#[-DAT,]
  TTO=TO[DAT]
  oobTO=TO
  X=dat[,c(-1,-2)]
  oobX=oobdat[,c(-1,-2)]
  ytt=dat[,1]
  oobytt=oobdat[,1]
  Z=dat[,2]
  group=oobdat[,2]###??????
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









