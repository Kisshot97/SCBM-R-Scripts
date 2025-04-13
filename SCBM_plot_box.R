################################################################################
# Module: Plot Box
# Description: This script is part of the SCBM framework and implements
# core functionality related to plot box.
################################################################################


# ====================
# Load Required Libraries
# ====================
library(reshape2)

# ====================
# Load Required Libraries
# ====================
library(Rmisc)

# ====================
# Load Required Libraries
# ====================
library(ggplot2)

# ====================
# Load Required Libraries
# ====================
library(grid)

# ====================
# Load Required Libraries
# ====================
library(cowplot)

# ====================
# Load Required Libraries
# ====================
library(gridExtra)
for (tp in c("bias","mse")) {
  

pl.type=tp#########bias or mse#########

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
datalist=bartlist=list()
nop=8
for (i in 1:nop) {

# --- Load CSV Data ---
  datalist[[i]]=read.csv(sprintf(paste0(getwd(), "/pdat/%s1_res%s.csv"),pl.type,i))[-1]######input######
}
dat=list()
qp <-as.numeric(c(1:(12*ncol(npnp0)))) 
for (np in qp){#####extract first 1000 samples
  ddad=as.numeric(npnp[1,np])
  nn=as.numeric(npnp[2,np])
  pp=as.numeric(npnp[3,np])
  nx=c(200,400,1000)
  px=c(25,50,100,150)
  aa=list()
  for (k in 1:nop) {
    aa[[k]]=na.omit(na.omit(datalist[[k]][np])[1:100,])
  }
# --- Construct a Data Frame ---
  aa=as.data.frame(as.list(aa))
  colnames(aa)=c("prop","pr_ps","stra","TO","BCM","BCM0","CF","VT")[1:nop]
  dat[[np]]=aa
}
# dat1=lapply(dat,function(df) {
# quan=function(x) {
#   ul=quantile(x, c(0,1))
#     low=which(x < ul[1])
#     up=which(x > ul[2])
#     x[low]=ul[1]
#     x[up]=ul[2]
#     #x=x[-c(low,up)]
# --- Return Final Output ---
#     return(x)}
# df1=numeric()
# for (i in 1:ncol(df)) {
#   df1=cbind(df1,quan(df[,i]))
# }
# colnames(df1)=colnames(df)
# --- Construct a Data Frame ---
# df1=as.data.frame(df1)
# --- Return Final Output ---
#   return(df1)
# })
dat1=dat
mdat=melt(dat1)
mdat=na.omit(mdat)
mdat[,2]=abs(mdat[,2])
mdat=mdat[-which(mdat$variable=="TO"),]####delete#####nops
pic1=list()


for (i in 1:12) {
ni=4*(i-1)+1
mf1=mdat[mdat[,3]==c(ni,ni+1,ni+2,ni+3),]
colnames(mf1)=c("method","MSE","p")
lable_name=list("n=200,p=25","n=400,p=50", "n=400,p=100", "n=1000,p=150")
lable_name_labeller=function(variable,value){
# --- Return Final Output ---
  return(lable_name[((as.numeric(value)-1)%%4+1)])
}
 name_mf1=unique(mf1[,1]) # model picking!!!
 mf1=mf1[is.element(mf1[,1],name_mf1),] # model picking!!!
 ymax=quantile(na.omit(mf1[,2]),c(0,0.95))# y scale!!
 
######pic ##################
pic1[[i]]=ggplot(mf1) +
  geom_boxplot(aes( y = MSE, x = method,fill=method)) +
  scale_y_continuous(trans = "sqrt", limits = ymax) +
  ylab(pl.type) +
  #scale_x_discrete(labels = c("n=100,p=50","n=200,p=100", "n=500,p=200", "n=500,p=400")) +
  xlab("scenario") +
  theme(legend.position = "none", text = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x  = element_text(size = 8), axis.text.y  = element_text(size = 15), plot.title = element_text(size = 20), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  ggtitle(sprintf("Simulation%s", i)) +
  scale_fill_manual(values = c("#0B60B0","#B7C9F2", "#D96098", "#1FE5BD", "#EA5455", "#CDF2CA", "#FFDEB4", "#57C5B6","#87A922","#E69F00", "#EE99C2","#FFF6E9","#12372A")) +
  facet_wrap(.~factor(p),labeller=lable_name_labeller,strip.position = "bottom")
}


# Create a DataFrame
# --- Construct a Data Frame ---
data <- data.frame(
  Xdata = rnorm(nop), Ydata = rnorm(nop),
  Methods =c("PTO","CM","CF","VT","CBART","SCBM","prop_ps","TO","TO.S","LR","LL","Tmars","Tmars.ps")[1:nop])
data$Methods= factor(data$Methods, levels =c("PTO","CM","CF","VT","CBART","SCBM","prop_ps","TO","TO.S","LR","LL","Tmars","Tmars.ps")[1:nop])
#ptomse,cmmse,cfmse,bartmse,cbartmse,propmse,prop_psmse,TOmse,TO.Smse,LRmse,LLmse,Tmarsmse,Tmars.psmse
# Create a Scatter Plot
gplot <- ggplot(data, aes(Xdata, Ydata, fill = Methods)) +
  geom_boxplot(size = rep(1,nop)) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 14)
  ) +guides(fill = guide_legend(nrow = 1))+
  scale_fill_manual(
    values = c("#0B60B0","#B7C9F2", "#D96098", "#1FE5BD", "#EA5455", "#CDF2CA", "#FFDEB4", "#57C5B6","#87A922","#E69F00", "#EE99C2","#FFF6E9","#12372A")
  )

# Draw Only Legend without plot
# Grab legend from gplot
legend <- get_legend(gplot)


# ====================
# Load Required Libraries
# ====================
library(gridExtra)

for (i in c(1,7)) {
  # 将多个图表组合成一个页面
  ppout = grid.arrange(pic1[[i]], pic1[[i+1]], pic1[[i+2]], pic1[[i+3]], pic1[[i+4]], pic1[[i+5]], ncol = 3)
  
  # 保存为 EPS 文件

# --- Save Plot to File ---
  ggsave(
    filename = paste0(paste0(getwd(), "/pic/"), pl.type, "0", i, ".eps"),
    plot = ppout,
    width = 15,
    height = 10,
    device = cairo_ps  # 使用 Cairo PostScript 设备生成 EPS 文件
  )
}


}